import numpy as np

def change_taxa_names(path_in:str):
    """This function replace the name of the nodes by the name of the taxas. 

    Args:
        path_in (str): Path to find the file. 
    """
    # read file
    with open(path_in, 'r') as f:
        text = f.read()

    list_lines = text.splitlines()

    # From 'translate' to ';'
    index_start = int(np.argwhere(np.array(list_lines)=='\tTRANSLATE')+1)
    index_stop = int(np.argwhere(np.array(list_lines[index_start:])=='\t;')+index_start)

    # Associate the number to the taxa
    corresponding_tables_index={}
    regex_number="[0-9][0-9]?"

    for i in range(index_start,index_stop):
        line = list_lines[i]
        matched_numbers = re.findall(regex_number,line)
        corresponding_tables_index[matched_numbers[0]]=f't{matched_numbers[1]}'

    # Tree
    tree=list_lines[-2]
    list_elt = tree.split(" ")
    list_values=list_elt[5][:-1]

    # Replace
    for index in list(corresponding_tables_index.keys()):
        list_values = re.sub(f'\({index}:',f'({corresponding_tables_index[index]}:',list_values)
        list_values = re.sub(f'\,{index}:',f',{corresponding_tables_index[index]}:',list_values)

    return(list_values)

def modify_branch_length(path:str, k:float, path_out:str):
    """This function multiply the branch length by a scalar.

    Args:
        path (str): path to the file.
        k (float): scalar.
    """
    # read file
    with open(path, 'r') as f:
        text = f.read()

    list_lines = text.splitlines()
    index_start = int(np.argwhere(np.array(list_lines)== '    <!-- specify the tree -->',)+1)

    # Find all the branch length and compute their new branch length
    corresponding_tables_index={}
    regex_number= "[0-9].[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]?" or "[0-9].[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]?" or  "[0-9].[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]?" or "[0-9].[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]?"
    line = list_lines[index_start]

    matched_numbers = re.findall(regex_number,line)

    for number in matched_numbers:
        corresponding_tables_index[number] = str(float(number)*k)
    
    # Replace
    for index in list(corresponding_tables_index.keys()):
        line = re.sub(index,corresponding_tables_index[index],line)
    
    # Ouvrez le fichier en mode écriture pour réécrire le contenu modifié
    list_lines[index_start] = line 

    list_lines = [ele + '\n' for ele in list_lines]

    print(path_out)
    with open(path_out, 'w') as new_file:
        new_file.writelines(list_lines)

    return
