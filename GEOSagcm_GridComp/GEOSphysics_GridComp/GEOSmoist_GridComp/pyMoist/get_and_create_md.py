import os

def get_python_module_names(directory):
    """
    Returns a list of Python module names (without .py extension) in the given directory.
    """
    return [
        os.path.splitext(f)[0]
        for f in os.listdir(directory)
        if f.endswith('.py') and os.path.isfile(os.path.join(directory, f))
    ]

# Example usage:
# module_names = get_python_module_names('/path/to/your/directory')

def create_md_files_from_list(strings, directory):
    """
    For each string in the list, creates a .md file named after the string,
    and writes the string as the content of the file.
    """
    os.makedirs(directory, exist_ok=True)
    strings.sort()
    for s in strings:
        if s != '__init__':
            filename = f"{s}.md"
            filepath = os.path.join(directory, filename)
            dir_name = filepath.split('/')[-2]
            with open(filepath, 'w') as f:
                f.write(f"# {s}\n\n")
                f.write(f"::: pyMoist.{dir_name}.{s}")
            f.close()

            print(f'- "{s}": {dir_name}/{s}.md')

# Example usage:
# names = ['foo', 'bar',

file_names = get_python_module_names('/Users/ckung/Documents/Code/MkDocs_playground/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/pyMoist/interface')
print(file_names)

create_md_files_from_list(file_names, '/Users/ckung/Documents/Code/MkDocs_playground/GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/docs/interface')