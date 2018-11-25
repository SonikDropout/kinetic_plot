# This module contains handlers for lines from chemkin input file
import re
from errors import inblockerror

def elemhandler(line, linenum, elements_dict, acceptable_elements_dict):
    """
    Handles elements data string, adds declared elements and thier atomic
    weights to elements_dict, raises error in case of invalid declaration
    """
    line = clearline(line)
    line = line.replace("/", " ")
    elements_list = re.findall(r'(\b[A-Z]{1,2}\b)(\s+\d+(\.\d+(E(-)?\d+)?)?)?', line)
    if elements_list == []:
       inblockerror('chemkin input', linenum)
    for elem in elements_list:
        atom = elem[0]
        weight = elem[1]
        if weight != '':
            weight = weight.lstrip()
            elements_dict[atom] = float(weight)
        elif atom in acceptable_elements_dict:
            elements_dict[atom] = acceptable_elements_dict[atom]
        else:
            inblockerror('chemkin input', linenum)


def spechandler(line, linenum, spec_dict, elem_dict):
    line = clearline(line)
    spec_list = line.split()
    spec_pattern = specpattern(elem_dict)
    for spec in spec_list:
        spec_dict[spec] = {}
        spec_dict[spec]['weight'] = 0 # TODO this may lead to dictionary clogging
        elements_in_spec = spec_pattern.findall(spec)
        for elem in elements_in_spec:
            atom = elem[0]
            coefficient = 1
            try:
                coefficient += int(elem[1]) - 1
            except ValueError:
                pass
            if atom in elem_dict:
                spec_dict[spec]['weight'] += elem_dict[atom] * coefficient
            else:
                del spec_dict[spec]
                inblockerror('chemkin input', linenum)


def specpattern(elements_dict):
    pattern = '('
    for atom in elements_dict:
        pattern += atom + '|'
    pattern = pattern.rstrip('|')
    pattern += ')(\\d{1,2})?'
    pattern = re.compile(pattern)
    return pattern


def therhandler(file, line, linenum, spec, spec_dict):
    try:
        thermlinenum = int(line[-2])
    except ValueError:
        inblockerror(file, linenum)
    if thermlinenum == 1:
        spec = parsefirstline(file, line, linenum, spec_dict)
    elif thermlinenum == 2:
        coefficients = parsecoefficients(file, line, linenum)
        spec_dict[spec]['coefficients_up'] = coefficients
    elif thermlinenum == 3:
        coefficients = parsecoefficients(file, line, linenum)
        spec_dict[spec]['coefficients_up'] += coefficients[0:2]
        spec_dict[spec]['coefficients_low'] = coefficients[2:]
    elif thermlinenum == 4:
        coefficients = parsecoefficients(file, line, linenum)
        spec_dict[spec]['coefficients_low'] += coefficients
    else:
        inblockerror(file, linenum)
    return spec


def parsecoefficients(file, line, linenum):
    str_coeffs = re.findall(r'-?0\.\d{8}E[+-]\d{2}', line)
    try:
        return list(map(float, str_coeffs))
    except ValueError:
        inblockerror(file, linenum)


def parsefirstline(file, line, linenum, spec_dict):
    try:
        spec = re.search(r'^\w+[+-]?', line).group()
        phase = re.search(r'G|L|S', line).group()
        spec_dict[spec]['phase'] = phase
        temperatures = re.findall(r'\d{4}\.\d{2}', line)
        map(float, temperatures)
        try:
            temps_low = float(temperatures[0]), float(temperatures[2])
            temps_up = float(temperatures[2]), float(temperatures[1])
            spec_dict[spec]['temp_ranges'] = [temps_low, temps_up]
        except IndexError:
            pass
        return spec
    except (ValueError, KeyError):
        inblockerror(file, linenum)


def reachandler(line, linenum, prevline, reactions, spec_dict):
    """
    parses reaction string, joins current line with
    the previous if necessary, appends the dictionary,
    containing reaction data, to reactions list
    :param line: line containing reaction data from chemkin input
    :param linenum: number of that line in file
    :param prevline: previous line from the file unless it didn't end with '&'
    :param reactions: reactions list to be complemented
    :return: prevline
    """
    line = clearline(line)
    if line[-2] == '&':
        prevline = line
        return prevline
    elif prevline != '':
        line = prevline.strip('&') + line
        prevline = ''
    reaction_dict = {}
    line.strip()
    reaction = re.match(r'([\w+.\s]+)<?=>?([\w+.\s]+[A-Z]\s)', line)
    reactants = parsereactiongroup(reaction.group(1), linenum, spec_dict)
    products = parsereactiongroup(reaction.group(2), linenum, spec_dict)
    coefficients = line.lstrip(reaction.group()).split()
    reaction_dict['A'] = coefficients[0]
    reaction_dict['beta'] = coefficients[1]
    reaction_dict['E'] = coefficients[2]
    reaction_dict['reactants'] = reactants
    reaction_dict['products'] = products
    reactions.append(reaction_dict)
    return prevline


def parsereactiongroup(reaction_part, linenum, spec_dict):
    """
    parses one of the reaction equation parts
    :param reaction_part: either left or right part of the reaction equation
    :return: list of reactants with tuples of reactants
             and their stoichiometric coefficients
    """
    reactants_list = []
    reactants = reaction_part.split('+')
    for reactant in reactants:
        reactant = reactant.strip()
        coefficient = re.match(r'[\d.]+', reactant)
        if coefficient: coefficient = coefficient.group()
        spec = reactant.lstrip(coefficient)
        if spec != 'M':
            if not spec in spec_dict:
                inblockerror('chemkin input', linenum)
            coefficient = float(coefficient or 1)
            reactants_list.append((spec, coefficient))
    return reactants_list


def clearline(line):
    """
    removes keywords and comments from chemkin line
    """
    pattern = r'((ELEM(ENT)?\s+)|(SPEC(IES)?\s+)|(\s*END\s*)|(!.*))*'
    line = re.sub(pattern, '', line)
    return line
