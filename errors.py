class Error(Exception):
    pass


def fileerror(file_path):
    raise Error(f'Can\'t open file \'{file_path}\'')


def inblockerror(file_name, linenum):
    raise Error(f'An error occurred in file \'{file_name}\' on line ({linenum})')


def missingblockerror(block_name):
    raise Error(f'An error occurred in ckemkin input: block \'{block_name}\' was not found')

