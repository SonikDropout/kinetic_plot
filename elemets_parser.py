import json
import re

def elements_from_json():
    weights_json = open("elementsWeights.json", "r")
    weights_data = weights_json.read()
    return json.loads(weights_data)

def chemkin_parser(acceptable_elemets_dict):
    chemkin_input = open("input/chemkinInput.txt", "r")
    chemkin_content = chemkin_input.read()
    chemkin_lines = chemkin_content.split("\n")
    remove_waste_lines(chemkin_lines)
    for chemkin_line in chemkin_lines:
        if chemkin_line.startswith("ELEM"):
            is_in_elements_block = True
            process_elements_string(chemkin_line)
        elif is_in_elements_block:
            process_elemets_string(chemkin_line)
        elif chemkin_line.startswith("SPEC"):
            is_in_elements_block = False
    print(chemkin_lines)

def remove_waste_lines(lines_array):
    number_of_empty_lines = lines_array.count("")
    for i in range(number_of_empty_lines):
        lines_array.remove("")
    for line in lines_array:
        if line.startswith("!"):
            lines_array.remove(line)
            
def process_elements_string(elements_string):
    elements_array = re.findall(r"\b[A-Za-z]{1,2}\b", elements_string)
    print(elements_array)
    #todo data processing
        
    
chemkin_parser(elements_from_json())
