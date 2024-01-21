# myscript.py

import sys

def process_file(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    updated_lines = [line.replace('(',"").replace(')',"").replace(','," ") for line in lines]

    with open(file_name, 'w') as file:
        file.write(''.join(updated_lines))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python myscript.py <file_name>")
        sys.exit(1)

    file_name = sys.argv[1]

    process_file(file_name)

