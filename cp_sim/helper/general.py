"""
 Title:         Helper
 Description:   General helper functions
 Author:        Janzen Choi

"""

# Libraries
import numpy as np, math, os

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

def round_sf(value:float, sf:int) -> float:
    """
    Rounds a float to a number of significant figures

    Parameters:
    * `value`: The value to be rounded
    * `sf`:    The number of significant figures

    Returns the rounded number
    """
    format_str = "{:." + str(sf) + "g}"
    rounded_value = float(format_str.format(value))
    return rounded_value

def dict_to_csv(data_dict:dict, csv_path:str, include_header:bool=True) -> None:
    """
    Converts a dictionary to a CSV file
    
    Parameters:
    * `data_dict`:      The dictionary to be converted
    * `csv_path`:       The path that the CSV file will be written to
    * `include_header`: Whether to include the header
    """
    
    # Extract and check headers
    headers = data_dict.keys()
    if len(headers) == 0:
        return

    # Turn all values into lists
    for header in headers:
        if not isinstance(data_dict[header], list):
            data_dict[header] = [data_dict[header]]
    
    # Open CSV file and write headers
    csv_fh = open(csv_path, "w+")
    if include_header:
        csv_fh.write(",".join(headers) + "\n")
    
    # Write data and close
    max_list_size = max([len(data_dict[header]) for header in headers])
    for i in range(max_list_size):
        row_list = [str(data_dict[header][i]) if i < len(data_dict[header]) else "" for header in headers]
        row_str = ",".join(row_list)
        csv_fh.write(row_str + "\n")
    csv_fh.close()

def csv_to_dict(csv_path:str, delimeter:str=",") -> dict:
    """
    Converts a CSV file into a dictionary
    
    Parameters:
    * `csv_path`:  The path to the CSV file
    * `delimeter`: The separating character
    
    Returns the dictionary
    """

    # Read all data from CSV (assume that file is not too big)
    csv_fh = open(csv_path, "r")
    csv_lines = csv_fh.readlines()
    csv_fh.close()

    # Initialisation for conversion
    csv_dict = {}
    headers = csv_lines[0].replace("\n", "").split(delimeter)
    csv_lines = csv_lines[1:]
    for header in headers:
        csv_dict[header] = []

    # Start conversion to dict
    for csv_line in csv_lines:
        csv_line_list = csv_line.replace("\n", "").split(delimeter)
        for i in range(len(headers)):
            value = csv_line_list[i]
            if value == "":
                continue
            try:
                value = float(value)
            except:
                pass
            csv_dict[headers[i]].append(value)
    
    # Convert single item lists to items and things multi-item lists
    for header in headers:
        if len(csv_dict[header]) == 1:
            csv_dict[header] = csv_dict[header][0]
    
    # Return
    return csv_dict

def get_sorted(value_list:list, reverse:bool=True) -> tuple:
    """
    Gets the top values and indexes of a list of values
    
    Parameters:
    * `value_list`: The list of values
    
    Returns the list of top values and indexes
    """
    sorted_value_list = sorted(value_list, reverse=reverse)
    sorted_index_list = []
    for value in sorted_value_list:
        for i in range(len(value_list)):
            if value == value_list[i] and not i in sorted_index_list:
                sorted_index_list.append(i)
                break
    return sorted_value_list, sorted_index_list

def transpose(list_of_lists:list) -> list:
    """
    Transposes a 2D list of lists
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the transposed list of lists
    """
    transposed = np.array(list_of_lists).T.tolist()
    return transposed

def flatten(list_of_lists:list) -> list:
    """
    Flattens a 2D list into a 1D list
    
    Parameters:
    * `list_of_lists`: A list of lists (i.e., a 2D grid)
    
    Returns the flattened list
    """
    return [item for sublist in list_of_lists for item in sublist]

def get_thinned_list(unthinned_list:list, density:int) -> list:
    """
    Gets a thinned list

    Parameters:
    * `unthinned_list`: The list before thinning
    * `density`:        The goal density of the thinned list

    Returns the thinned list
    """
    src_data_size = len(unthinned_list)
    step_size = src_data_size / density
    thin_indexes = [math.floor(step_size*i) for i in range(1, density - 1)]
    thin_indexes = [0] + thin_indexes + [src_data_size - 1]
    thinned_list = [unthinned_list[i] for i in thin_indexes]
    return thinned_list

def safe_mkdir(dir_path:str) -> None:
    """
    For safely making a directory

    Parameters:
    * `dir_path`: The path to the directory
    """
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        pass

