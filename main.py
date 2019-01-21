from main_parser import open_parse_chemkin_input
from json import loads
import numpy as np
from scipy.integrate import odeint
from math import exp
import matplotlib.pyplot as plt
import os


class KineticMechanism:

    def __init__(self, name, chemkin_path, initial_conditions_path):
        self.name = name
        self.elements, self.species, self.reactions = self.parse_input(chemkin_path)
        self.__initial_conditions = self.get_initial_conditions_from(initial_conditions_path)
        self.max_time = 1000

    def parse_input(self, path_to_file):
        return open_parse_chemkin_input(path_to_file)

    def get_initial_conditions_from(self, path_to_file):
        try:
            with open(path_to_file, 'r') as input_conditions:
                return loads(input_conditions.read())
        except OSError as e:
            print(e)

    def make_calculations(self):
        matrices = ReactionMatrices(self.species, self.reactions, self.__initial_conditions['temperature'])
        system = ODESystem(matrices, self.__initial_conditions, self.max_time)
        return system.solve()

    def set_max_time(self, time):
        self.max_time = time

    def build_plot(self):
        data, t = self.make_calculations()
        data = data.T
        for i, specie in enumerate(self.species):
            if np.amax(data[i]) > 0.01:
                plt.plot(t, data[i], label=specie)
        plt.legend(loc='best')
        self.decorate_plot()
        plt.show()
        # self.save_figure()

    def decorate_plot(self):
        plt.grid(True, 'both', alpha=0.5)
        plt.xlabel('time, s')
        plt.xscale('log')
        plt.yticks(np.linspace(0, 1, 11))
        plt.ylabel('Molar fraction')
        plt.ylim(0, 1)
        plt.title(f'T={self.__initial_conditions["temperature"]}K')

    def save_figure(self):
        plot_path = os.path.join('output/plots', f'{self.name}_T={self.__initial_conditions["temperature"]}K.png')
        plt.savefig(plot_path, format='png')


class ODESystem:

    def __init__(self, matrices, initial_conditions, max_time=1000):
        self.matrices = matrices
        self.initial_conditions = initial_conditions
        self.max_time = max_time

    def solve(self):
        C0 = self.get_initial_concentrations()
        t = np.linspace(0, self.max_time, 10000000)
        dC_dt = self.get_model()
        solution = odeint(dC_dt, C0, t)
        return solution, t

    def get_initial_concentrations(self):
        return np.array([value for value in self.initial_conditions['concentrations'].values()])

    def get_model(self):
        build_equation = self.build_equation

        def model(C, t=0):
            return np.array([build_equation(C, i) for i in range(len(C))])
        return model

    def build_equation(self, C, i):
        equation = 0
        for j in range(self.matrices.reactants_stoich.shape[1]):
            C_mul_reac, C_mul_prod = 1, 1
            for k in range(self.matrices.reactants_stoich.shape[0]):
                C_mul_reac *= C[k]**self.matrices.reactants_stoich[k][j]
                C_mul_prod *= C[k]**self.matrices.products_stoich[k][j]
            stoich = (self.matrices.products_stoich[i][j] - self.matrices.reactants_stoich[i][j])
            equation += stoich * (self.matrices.k_for[j] * C_mul_reac - self.matrices.k_rev[j] * C_mul_prod)

        return equation


class ReactionMatrices:

    def __init__(self, species, reactions, temperature):
        self.__reactions = reactions
        self.__species = species
        self.__temperature = temperature
        self.k_for, self.k_rev = self.get_K_vectors()
        self.reactants_stoich, self.products_stoich = self.get_stoich()

    @staticmethod
    def calculate_K(parameters):
        return parameters['A']*parameters['T']**parameters['beta']*exp(-parameters['E']/parameters['T'])

    def get_K_vectors(self):
        k_for = np.empty(len(self.__reactions)-1)
        k_rev = np.empty(len(self.__reactions)-1)
        for i in range(1, len(self.__reactions)):
            for_params = self.__reactions[i]['FOR_ARR_COEFFS']
            rev_params = self.__reactions[i]['REV_ARR_COEFFS']
            for_params.update({"T": self.__temperature})
            rev_params.update({"T": self.__temperature})
            k_for[i-1] = self.calculate_K(for_params)
            k_rev[i-1] = self.calculate_K(rev_params)
        return k_for, k_rev

    def get_stoich(self):
        shape = (len(self.__reactions) - 1, len(self.__species))
        reactants_matrix = np.zeros(shape)
        products_matrix = np.zeros(shape)
        self.fill_matrices(reactants_matrix, products_matrix)
        return reactants_matrix.T.astype('int'), products_matrix.T.astype('int')

    def fill_matrices(self, reactants_matrix, products_matrix):
        for i in range(1, len(self.__reactions)):
            reactants = self.__reactions[i]['reactants']
            products = self.__reactions[i]['products']
            reactants_matrix[i - 1] = self.map_species_list_to(reactants)
            products_matrix[i - 1] = self.map_species_list_to(products)

    def map_species_list_to(self, list):
        result = np.zeros(len(self.__species))
        for index, specie in enumerate(self.__species):
            for item in list:
                if specie in item:
                    result[index] += item[1]
        return result


class UserDialog:

    def __init__(self):
        self.__available_files = {
            "evans": ("KINETIC_MECHANISM_EVANS.DAT", "INITIAL_EVANS.json"),
            "jakimovski": ("KINETIC_MECHANISM_JAKIMOVSKI.DAT", "INITIAL_JAKIMOVSKI.json"),
            "jakimovski_short": ("KINETIC_MECHANISM_JAK_SHORT.DAT", "INITIAL_JAK_SHORT.json")
            }
        self.__finished = False

    def ask_for_mech(self):
        print('\nAvailable kinetic mechanisms:')
        for i, key in enumerate(self.__available_files):
            print(f'    [{i}] {key}')
        return int(input(f'\nSelect mechanism[0-{i}]: '))

    @staticmethod
    def ask_time():
        return float(input('Enter max time in seconds: '))

    def get_paths(self):
        mech_index = self.ask_for_mech()
        mech_name = list(self.__available_files.keys())[mech_index]
        file_names_list = list(self.__available_files.values())
        return mech_name, \
            os.path.join('input', file_names_list[mech_index][0]), \
            os.path.join('input', file_names_list[mech_index][1])

    def start(self):
        while not self.__finished:
            mechanism_name, mechanism_path, initial_conditions_path = self.get_paths()
            self.confirm_edit(initial_conditions_path)

            mechanism = KineticMechanism(mechanism_name, mechanism_path, initial_conditions_path)
            mechanism.set_max_time(self.ask_time())

            mechanism.build_plot()

            self.confirm_exit()

    def confirm_exit(self):
        answer = input('\nProceed to other mechanism[y/N]: ') or 'N'
        if answer == 'N':
            self.__finished = True

    @staticmethod
    def confirm_edit(initial_conditions_path):
        with open(initial_conditions_path, 'r') as ini_conditions:
            print(f'\n{ini_conditions.read()}')
        answer = input('Edit initial conditions JSON[y/N]: ') or 'N'
        if answer == 'y':
            os.system(f'vim {initial_conditions_path}')


if __name__ == '__main__':

    dialog = UserDialog()
    dialog.start()

