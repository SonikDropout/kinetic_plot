# This module contains handlers for lines from chemkin input file
import re

def elemhandler(line, linenum, elements_dict, acceptable_elements_dict):
    """
    Handles elements data string, adds declared elements and thier atomic
    weights to elements_dict, raises error in case of invalid declaration
    """
    line = clearline(line)
    line = line.replace("/", " ")
    elements_list = re.findall(r'(\b[A-Z]{1,2}\b)(\s+\d+(\.\d+(E(-)?\d+)?)?)?', line)
    if elements_list == []:
        raise Error(f'An error occured in chemkin input on line ({linenum}):' \
                    f' invalid element \'{line}\'')
    for elem in elements_list:
        atom = elem[0]
        weight = elem[1]
        if weight != '':
            weight = weight.lstrip()
            elements_dict[atom] = float(weight)
        elif atom in acceptable_elements_dict:
            elements_dict[atom] = acceptable_elements_dict[atom]
        else:
            raise Error(f'An error occured in chemkin input on line ({linenum}):' \
                        f' invalid element \'{elem}\'')


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
                raise Error(f'An error occured in chemkin input on line ({linenum}): \
                            spec \'{spec}\' contains undefined element \'{atom}\'')


def specpattern(elements_dict):
    pattern = '('
    for atom in elements_dict:
        pattern += atom + '|'
    pattern = pattern.rstrip('|')
    pattern += ')(\\d{1,2})?'
    pattern = re.compile(pattern)
    return pattern


def therhandler(line, spec, spec_dict):
    thermlinenum = int(line[-2])
    if thermlinenum == 1:
        spec = parsefirstline(line, spec_dict)
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
    str_coeffs = re.findall(r'-?0\.\d{8}E[+-]\d{2}', line)
    return list(map(float, str_coeffs))


def parsefirstline(line, spec_dict):
    spec = re.search(r'^\w+[+-]?', line).group()
    phase = re.search(r'G|L|S', line).group()
    spec_dict[spec]['phase'] = phase
    temperatures = re.findall(r'\d{4}\.\d{2}', line)
    map (float, temperatures)
    try :
        temps_low = float(temperatures[0]), float(temperatures[2])
        temps_up = float(temperatures[2]), float(temperatures[1])
        spec_dict[spec]['temp_ranges'] = [temps_low, temps_up]
    except IndexError:
        pass
    return spec



def clearline(line):
    """
    removes keywords and comments from chemkin line
    """
    patterns = [r'ELEM(ENT)?\s+', r'SPEC(IES)?\s+', r'\s*END\s*', r'!.*']
    for pattern in patterns:
        line = re.sub(pattern, '', line)
    return line
