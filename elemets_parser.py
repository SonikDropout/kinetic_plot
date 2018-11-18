from json import loads
import re

class Error(Exception):
    pass

class InputError(Error):
    pass

def elements_from_json():
    weights_json = open("elementsWeights.json", "r")
    weights_data = weights_json.read()
    return loads(weights_data)

def chemkin_parser(acceptable_elements_dict):
    chemkin_input = open("input/chemkinInput.txt", "r")
    elements_dict = {}
    species_dict = {}
    is_in_elem = False
    is_in_spec = False
    is_in_ther = False
    is_in_reac = False
    for line in chemkin_input:
        if line.startswith("!") or line.strip() == "":
            continue
        elif inelem(line, is_in_elem):
            is_in_elem = True
            elemhandler(line, elements_dict, acceptable_elements_dict)
        elif inspec(line, is_in_spec, is_in_elem):
            is_in_spec = True
            spechandler(line, species_dict, elements_dict)
            
def elemhandler(elements_string, elements_dict, acceptable_elements_dict):
    elements_string = elements_string.replace("/", "")
    elements_list = elements_string.split()
    for i in range(len(elements_list)):
        element = elements_list[i]
        if element in acceptable_elements_dict:
            elements_dict[element] = acceptable_elements_dict[element]
        elif re.match(r"\d.*", element):
            elements_dict[elements_list[i-1]] = float(element) # TODO not float type handler
        # TODO error handling
    print(elements_dict)
    
def spechandler(line, spec_dict, elem_dict):
    pass
    # TODO handle species data
    
def inelem(line, flag):
    return (line.startswith("ELEM") or flag) and (not line.startswith("SPEC"))

def inspec(line, flag, prevflag):
    if prevflag:
        prevflag = False
    else:
        raise InputError
    return (line.startswith("SPEC") or flag) and (not line.startswith("TERM"))
    
def clearline(patterns, string):
    for pattern in patterns:
        string = re.sub(pattern, "", string)
        
try:    
    chemkin_parser(elements_from_json())
except InputError:
    print("Invalid input")
