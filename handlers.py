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
        spec_dict[spec] = 0 # TODO this may lead to dictionary clogging
        elements_in_spec = spec_pattern.findall(spec)
        for elem in elements_in_spec:
            atom = elem[0]
            coefficient = 1
            try:
                coefficient += int(elem[1]) - 1
            except ValueError:
                pass
            if atom in elem_dict:
                spec_dict[spec] += elem_dict[atom] * coefficient
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


def clearline(line):
    """
    removes keywords and comments from chemkin line
    """
    patterns = [r'ELEM(ENT)?\s+', r'SPEC(IES)?\s+', r'\s*END\s*', r'!.*']
    for pattern in patterns:
        line = re.sub(pattern, '', line)
    return line
