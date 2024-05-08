# Libraries
import numpy as np, os, math
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep

# Constants
X_HEADER        = "strain"
Y_HEADER        = "stress"
NUM_POINTS      = 30
PARAM_NAME_LIST = ["tau_sat", "b", "tau_0", "gamma_0", "n"]

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

# The Interpolator Class
class Interpolator:

    def __init__(self, x_list:list, y_list:list, resolution:int=50, smooth:bool=False):
        """
        Class for interpolating two lists of values

        Parameters:
        * `x_list`:     List of x values
        * `y_list`:     List of y values
        * `resolution`: The resolution used for the interpolation
        * `smooth`:     Whether to smooth the interpolation
        """
        x_list, indices = np.unique(np.array(x_list), return_index=True)
        y_list = np.array(y_list)[indices]
        if len(x_list) > resolution:
            x_list = get_thinned_list(list(x_list), resolution)
            y_list = get_thinned_list(list(y_list), resolution)
        smooth_amount = resolution if smooth else 0
        self.spl = splrep(x_list, y_list, s=smooth_amount)

    def evaluate(self, x_list:list) -> list:
        """
        Run the interpolator for specific values

        Parameters
        * `x_list`: The list of x values

        Returns the evaluated values
        """
        return list(splev(x_list, self.spl))

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

    Returns the new dictionary
    """
    for key in value_dict.keys():
        list_dict[key].append(value_dict[key])
    return list_dict

# Initialise success dictionary
success_keys = PARAM_NAME_LIST + ["x_end"] + [f"y_{i+1}" for i in range(NUM_POINTS)]
success_dict = {}
for key in success_keys:
    success_dict[key] = []

# Initialise failure dictionary
failure_keys = PARAM_NAME_LIST + ["reason"]
failure_dict = {}
for key in failure_keys:
    failure_dict[key] = []

# Read all CSV files and iterate through them
results_dir = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/20240508 (tensile bcc)"
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
            failure_dict = append_list_dict(failure_dict, param_dict)
            failure_dict["reason"].append(keyword)
            failed = True
    if failed:
        continue

    # Check whether the number of datapoints are sufficient
    if len(data_dict[X_HEADER]) < 5:
        failure_dict = append_list_dict(failure_dict, param_dict) 
        failure_dict["reason"].append("insufficient")
        continue

    # Format curve informationn
    interpolator = Interpolator(data_dict[X_HEADER], data_dict[Y_HEADER])
    x_end  = max(data_dict[X_HEADER])
    x_list = [x_end/(NUM_POINTS-1)*j for j in range(NUM_POINTS)]
    y_list = interpolator.evaluate(x_list)

    # Plot curves
    plt.scatter(data_dict[X_HEADER], data_dict[Y_HEADER], color="silver")
    plt.plot(x_list, y_list)

    # Add curve information to success dictionary
    success_dict = append_list_dict(success_dict, param_dict)
    success_dict["x_end"].append(round_sf(x_end, 5))
    for i in range(NUM_POINTS):
        success_dict[f"y_{i+1}"].append(round_sf(y_list[i], 5))

# Write results
dict_to_csv(failure_dict, "failures.csv")
dict_to_csv(success_dict, "tc.csv")

# Format and save plot for fits
plt.xlim(0.0, 0.5)
plt.ylim(0.0, None)
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.savefig("tc_plot.png")
plt.clf()
