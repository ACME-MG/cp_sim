# Libraries
import numpy as np, os

# Constants
X_HEADER        = "strain"
PHI_HEADERS     = ["phi_1", "Phi", "phi_2"]
GRAIN_INDEX     = 1 # 0
NUM_POINTS      = 10
PARAM_NAME_LIST = ["tau_sat", "b", "tau_0", "gamma_0", "n", "cd", "beta"]

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
        else:
            csv_dict[header] = csv_dict[header]
    
    # Return
    return csv_dict

def dict_to_csv(data_dict:dict, csv_path:str) -> None:
    """
    Converts a dictionary to a CSV file
    
    Parameters:
    * `data_dict`: The dictionary to be converted
    * `csv_path`: The path that the CSV file will be written to
    """
    
    # Extract headers and turn all values into lists
    headers = data_dict.keys()
    for header in headers:
        if not isinstance(data_dict[header], list):
            data_dict[header] = [data_dict[header]]
    
    # Open CSV file and write headers
    csv_fh = open(csv_path, "w+")
    csv_fh.write(",".join(headers) + "\n")
    
    # Write data and close
    max_list_size = max([len(data_dict[header]) for header in headers])
    for i in range(max_list_size):
        row_list = [str(data_dict[header][i]) if i < len(data_dict[header]) else "" for header in headers]
        row_str = ",".join(row_list)
        csv_fh.write(row_str + "\n")
    csv_fh.close()

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

def append_list_dict(list_dict:dict, value_dict:dict) -> dict:
    """
    Appends values from a dictionary to lists in another dictionary
    
    Parameters:
    * `list_dict`:  The dictionary of lists
    * `value_dict`: The dictionary of values

    Returnns the new dictionary
    """
    for key in value_dict.keys():
        list_dict[key].append(value_dict[key])
    return list_dict

# Initialise success dictionary
success_keys = PARAM_NAME_LIST + ["phi"] + [f"y_{i+1}" for i in range(NUM_POINTS)]
success_dict = {}
for key in success_keys:
    success_dict[key] = []

# Read all CSV files and iterate through them
# results_dir = "../results"
results_dir = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/20240422 (tensile)"
csv_file_list = [file for file in os.listdir(results_dir) if file.endswith(".csv")]

# # Only retrieve a subset of the CSVs (for debugging)
# import random
# csv_file_list = list(random.sample(csv_file_list, 100))

# Iterate through CSV files
for csv_file in csv_file_list:
    
    # Convert csv file to dictionary
    data_dict = csv_to_dict(f"{results_dir}/{csv_file}")

    # Get parameter informationn
    param_dict = {}
    for param_name in PARAM_NAME_LIST:
        param_dict[param_name] = data_dict[param_name]

    # Check whether the simulation failed or timed out
    failed = False
    for keyword in ["failed", "timeout"]:
        if keyword in csv_file:
            failed = True
    if failed:
        continue

    # Add parameter information
    for _ in range(3):
        success_dict = append_list_dict(success_dict, param_dict)

    # Get the trajectories for one of the grains
    poly_str_list = [data_dict[header][GRAIN_INDEX] for header in PHI_HEADERS]
    poly_list = [[float(coef) for coef in poly_str.split(" ")] for poly_str in poly_str_list]
    
    # Get points on each trajectory and add to success dictionary
    for i in range(len(poly_list)):
        x_end = len(data_dict[X_HEADER])
        x_list = [x_end/(NUM_POINTS-1)*j for j in range(NUM_POINTS)]
        y_list = np.polyval(poly_list[i], x_list)
        success_dict["phi"].append(i+1)
        for i in range(NUM_POINTS):
            success_dict[f"y_{i+1}"].append(round_sf(y_list[i], 5))

# Write results
dict_to_csv(success_dict, "phi_success.csv")
