from json import loads
from handlers import elemhandler, spechandler


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
        elif inspec(line, is_in_spec):
            if is_in_elem:
                is_in_elem = False
#            else:
#                raise Error(f'An error occured in chemkin input on line ({linenum}):' \
#                    'SPECIES block declared before ELEMENTS')
                #TODO correct error raising
            is_in_spec = True
            spechandler(line, linenum, species_dict, elements_dict)


def inelem(line, flag):
    return (line.startswith("ELEM") or flag) and (not line.startswith("SPEC"))


def inspec(line, flag):
    return (line.startswith("SPEC") or flag) and (not line.startswith("TERM"))


def main():
    try:
        with open("input/chemkinInput.txt", "r") as chemkin_input:
            try:
                chemkin_parser(chemkin_input, elements_from_json())
            except Error as err:
                print(err)
    except OSError:
        print('Cannot open chemkin input file')


main()
