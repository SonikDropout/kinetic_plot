from json import loads
import re

class Error(Exception):
    pass

def elements_from_json():
    weights_json = open("elementsWeights.json", "r")
    weights_data = weights_json.read()
    return loads(weights_data)

def chemkin_parser(input_file, acceptable_elements_dict):
    """ 
    Parses chemkin input file line by line
    """
    elements_dict = {}
    species_dict = {}
    is_in_elem = False
    is_in_spec = False
    is_in_ther = False
    is_in_reac = False
    linenum = 0
    for line in input_file:
        line = line.upper()
        linenum += 1
        if line.startswith("!") or line.strip() == "":
            continue
        elif inelem(line, is_in_elem):
            is_in_elem = True
            elemhandler(line, linenum, elements_dict, acceptable_elements_dict)
        elif inspec(line, linenum, is_in_spec, is_in_elem):
            if is_in_elem:
                is_in_elem = False
#            else:
#                raise Error(f'An error occured in chemkin input on line ({linenum}):' \
#                    'SPECIES block declared before ELEMENTS')
                #TODO correct error raising
            is_in_spec = True
            spechandler(line, linenum, species_dict, elements_dict)
            
def elemhandler(line, linenum, elements_dict, acceptable_elements_dict):
    """
    Handles elements data string, adds declared elements and thier atomic
    weights to elements_dict, raises error in case of invalid declaration
    """
    line = clearline(line)
    line = line.replace("/", " ")
    elements_list = line.split()
    for i in range(len(elements_list)):
        element = elements_list[i]
        if re.fullmatch(r'[A-Z]{1,2}', element):
            if element in acceptable_elements_dict:
                elements_dict[element] = acceptable_elements_dict[element]
            else:
                try:
                    re.fullmatch(r"\d+(\.\d+(E(-)?\d+)?)?", elements_list[i+1])
                except IndexError:
                    raise Error(f'An error occured in chemkin input on ({linenum}): \
                                  isotop \'{element}\' has undecrlared weight')
        elif re.fullmatch(r"\d+(\.\d+(E(-)?\d+)?)?", element):
            try:
                if re.fullmatch(r'[A-Z]{1,2}', elements_list[i-1]):
                    elements_dict[elements_list[i-1]] = float(element)
            except IndexError:
                raise Error(f'An error occured in chemkin input on line ({linenum}): \
                              invalid element \'{element}\'')
        else:
            raise Error(f'An error occured in chemkin input on line ({linenum}): \
                          invalid element \'{element}\'')
    
def spechandler(line, linenum, spec_dict, elem_dict):
    line = clearline(line)
    spec_list = line.split()
    for spec in spec_list:
        spec_dict[spec] = 0 #TODO this may lead to dictionary clogging
        elements_in_spec = re.findall(r'([A-Z]{1,2})(\d{1,2})?', spec) #TODO compile another pattern form elements dictionary
        print(re.match(r'([A-Z]{1,2})(\d{1,2})?', spec))
        for elem in elements_in_spec:
            atom = elem[0]
            coefficient = 1
            try:
                coefficient += int(elem[1])
            except ValueError:
                pass
            if atom in elem_dict:
                spec_dict[spec] += elem_dict[atom] * coefficient
            else:
                del spec_dict[spec]
                raise Error(f'An error occured in chemkin input on line ({linenum}): \
                            spec \'{spec}\' contains undefined element \'{atom}\'')
    
def inelem(line, flag):
    return (line.startswith("ELEM") or flag) and (not line.startswith("SPEC"))

def inspec(line, linenum, flag, prevflag):
    return (line.startswith("SPEC") or flag) and (not line.startswith("TERM"))
    
def clearline(line):
    """
    removes keywords and comments from chemkin line
    """
    patterns = [r'ELEM(ENT)?\s+', r'SPEC(IES)?\s+', r'\s*END\s*', r'!.*']
    for pattern in patterns:
        line = re.sub(pattern, '', line)
    return line

def main():
    chemkin_input = open("input/chemkinInput.txt", "r")        
    try:    
        chemkin_parser(chemkin_input, elements_from_json())
    except Error as err:
        print(err)

main()
