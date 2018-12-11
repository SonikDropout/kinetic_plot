# This module contains handlers for lines from chemkin input file
import re
from errors import inblockerror


def parse_elements_line(line, elements_dict, acceptable_elements_dict):
    """
    Handles elements data string, adds declared elements and thier atomic
    weights to elements_dict, raises error in case of invalid declaration
    """
    line = clearline(line)
    if line == '':
        return
    line = line.replace("/", " ")
    elements_list = re.findall(r'(\b[A-Z]{1,2}\b)(\s+\d+(\.\d+(E(-)?\d+)?)?)?', line)
    for elem in elements_list:
        atom = elem[0]
        weight = elem[1]
        if weight != '':
            weight = weight.lstrip()
            elements_dict[atom] = float(weight)
        elif atom in acceptable_elements_dict:
            elements_dict[atom] = acceptable_elements_dict[atom]


def parse_spec_line(line, spec_dict, elem_dict):
    line = clearline(line)
    spec_list = line.split()
    spec_pattern = get_spec_pattern(elem_dict)
    for spec in spec_list:
        spec_dict[spec] = {}
        spec_dict[spec]['weight'] = 0
        elements_in_spec = spec_pattern.findall(spec)
        for elem in elements_in_spec:
            atom = elem[0]
            if elem[1]:
                coefficient = int(elem[1])
            else:
                coefficient = 1
            if atom in elem_dict:
                spec_dict[spec]['weight'] += elem_dict[atom] * coefficient
            else:
                del spec_dict[spec]


def get_spec_pattern(elements_dict):
    pattern = '('
    for atom in elements_dict:
        pattern += atom + '|'
    pattern = pattern.rstrip('|')
    pattern += ')(\\d{1,2})?'
    pattern = re.compile(pattern)
    return pattern


def parse_therm_line(line, spec, spec_dict):
    thermlinenum = int(line[-2])
    if thermlinenum == 1:
        spec = parse_first_therm_line(line, spec_dict)
    elif thermlinenum == 2:
        coefficients = parsecoefficients(line)
        spec_dict[spec]['coefficients_up'] = coefficients
    elif thermlinenum == 3:
        coefficients = parsecoefficients(line)
        spec_dict[spec]['coefficients_up'] += coefficients[0:2]
        spec_dict[spec]['coefficients_low'] = coefficients[2:]
    elif thermlinenum == 4:
        coefficients = parsecoefficients(line)
        spec_dict[spec]['coefficients_low'] += coefficients
    return spec


def parsecoefficients(line):
    str_coeffs = re.findall(r'-?[0-9]\.\d{8}E[+-]\d{2}', line)
    return list(map(float, str_coeffs))


def parse_first_therm_line(line, spec_dict):
        spec = re.match(r'^\w+[+-]?\s', line).group().rstrip()
        line = line[-37:]
        phase = re.search(r'[GLS]', line).group()
        spec_dict[spec]['phase'] = phase
        temperatures = re.findall(r'\d+\.\d*', line)
        temperatures = list(map(float, temperatures))
        temps_low = temperatures[0], temperatures[2]
        temps_up = temperatures[2], temperatures[1]
        spec_dict[spec]['temp_ranges'] = [temps_low, temps_up]
        return spec


def parse_units(line):
    return re.findall(r'[A-Z/]+', line.lstrip('REACTIONS'))


def parse_reac_line(line, linenum, reaction_index, reactions, spec_dict):
    """
    parses reaction string, joins current line with
    the previous if necessary, appends the dictionary,
    containing reaction data, to reactions list
    """
    line = line.strip()
    is_line_of_parameters = get_other_arrenius_coefficients(line, reaction_index, reactions)
    if is_line_of_parameters:
        return reaction_index
    reaction_dict = {}
    splited_line = line.split()
    coefficients = splited_line[-3:]
    reaction = ' '.join(splited_line[:-3])
    reaction = re.match(r'([\w+.\s]+)<?=>?([\w+.\s]+)', reaction)
    reactants = parsereactiongroup(reaction.group(1), linenum, spec_dict)
    products = parsereactiongroup(reaction.group(2), linenum, spec_dict)
    reaction_dict['FOR_ARR_COEFS'] = get_arrenius_coefficients(coefficients)
    reaction_dict['reactants'] = reactants
    reaction_dict['products'] = products
    reactions.append(reaction_dict)
    reaction_index += 1
    return reaction_index


def get_other_arrenius_coefficients(line, reaction_index, reactions):
    for key_word in ['REV', 'TROE', 'LOW']:
        if line.startswith(key_word):
            coefficients = re.findall(r'(-?\d+\.?(?:\d+(?:E[+-])?\d+)?)', line)
            coefficients = list(map(lambda tup: ''.join(tup), coefficients))
            reactions[reaction_index][key_word + '_ARR_COEFFS'] = get_arrenius_coefficients(coefficients)
            return True
    return False


def get_arrenius_coefficients(coefficients):
    index = 0
    parameters_dict = {}
    for parameter_name in ['A', 'beta', 'E']:
        parameter_value = float(coefficients[index])
        parameters_dict[parameter_name] = parameter_value
        index += 1
    return parameters_dict


def parsereactiongroup(reaction_part, linenum, spec_dict):
    """
    parses one of the reaction equation parts
    """
    reactants_list = []
    reactants = reaction_part.split('+')
    for reactant in reactants:
        reactant = reactant.strip()
        coefficient = re.match(r'[\d.]+', reactant)
        if coefficient:
            coefficient = coefficient.group()
        spec = reactant.lstrip(coefficient)
        if spec != 'M':
            if spec not in spec_dict:
                inblockerror('chemkin input', linenum)
            coefficient = float(coefficient or 1)
            reactants_list.append((spec, coefficient))
    return reactants_list


def clearline(line):
    """
    removes keywords and comments from chemkin line
    """
    pattern = r'((ELEM(ENTS)?\s+)|(SPEC(IES)?\s+)|(\s*END\s*)|(!.*))*'
    line = re.sub(pattern, '', line).strip()
    return line
