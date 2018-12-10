from json import loads
from parsers import *
from errors import *


def elements_from_json():
    """
    parses json with elements' atomic weight data
    :return: dictionary of element-weight pairs
    """
    weights_json = open("elementsWeights.json", "r")
    weights_data = weights_json.read()
    return loads(weights_data)


def chemkin_parser(input_file, acceptable_elements_dict):
    """ 
    Parses chemkin input file line by line
    """
    elements_dict = {}
    species_dict = {}
    reactions = []
    spec = ''
    reaction_index = -1
    flags = {}
    flags_names = ['elem', 'spec', 'ther', 'reac']
    for flag_name in flags_names:
        flags['is_in_' + flag_name] = False
        flags['was_in_' + flag_name] = False
    line_number = 0
    for line in input_file:
        try:
            line = line.upper()
            line_number += 1
            if line.startswith("!") or line.strip() == "" or line.startswith("END"):
                continue
            elif inelem(line, flags['is_in_elem']):
                flags['is_in_elem'] = True
                flags['was_in_elem'] = True
                parse_elements_line(line, elements_dict, acceptable_elements_dict)
            elif inspec(line, flags['is_in_spec']):
                if flags['is_in_elem']:
                    flags['is_in_elem'] = False
                elif not flags['was_in_elem']:
                    missingblockerror("ELEMENTS")
                flags['is_in_spec'] = True
                flags['was_in_spec'] = True
                parse_spec_line(line, species_dict, elements_dict)
            elif line.startswith("THER"):
                flags['is_in_spec'] = False
                flags['is_in_ther'] = True
                flags['was_in_ther'] = True
                if "ALL" not in line:
                    open_parse_themo_data('input/thermo.dat', species_dict)  # TODO file path is hardcoded
            elif flags['is_in_ther'] and not line.startswith("REAC"):
                spec = parse_therm_line(line, spec, species_dict)
            elif line.startswith("REAC"):
                if not flags['was_in_ther']:
                    if flags['was_in_spec']:
                        open_parse_themo_data('input/thermo.dat', species_dict)
                        flags['is_in_spec'] = False
                    else:
                        missingblockerror("SPECIES")
                flags['is_in_ther'] = False
                flags['is_in_reac'] = True
                flags['was_in_reac'] = True
            elif flags['is_in_reac']:
                reaction_index = parse_reac_line(line, line_number, reaction_index, reactions, species_dict)
        except:
            inblockerror(input_file.name, line_number)
    return


def inelem(line, flag):
    return (line.startswith("ELEM") or flag) and (not line.startswith("SPEC"))


def inspec(line, flag):
    return (line.startswith("SPEC") or flag) and (not line.startswith("THER") and not line.startswith("REAC"))


def open_parse_themo_data(path_to_file, species_dict):
    """
    tries opening thermodynamic data file and parsing it
    :param path_to_file: path to themdodynamic '.dat' file
    :param species_dict: dictionary of species for thermodynamic data to be appended
    """
    try:
        with open(path_to_file, "r") as thermo_data:
            linenum = 0
            spec = ''
            for line in thermo_data:
                linenum += 1
                spec = parse_therm_line(line, spec, species_dict)
    except OSError:
        fileerror(path_to_file)


def main():
    try:
        with open("input/KINETIC_MECHANISM_EVANS.DAT", "r") as chemkin_input:
            try:
                chemkin_parser(chemkin_input, elements_from_json())
            except Error as err:
                print(err)
    except OSError:
        print('Cannot open chemkin input file')
        pass


main()
