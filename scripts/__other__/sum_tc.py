# Libraries
import numpy as np, os, math
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
import sys; sys.path += ["../.."]
from cp_sim.helper import csv_to_dict, dict_to_csv, round_sf

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
results_dir = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/20240515 (cp_bcc 1ss)"
csv_file_list = [file for file in os.listdir(results_dir) if file.endswith(".csv")]

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
dict_to_csv(failure_dict, "summary/failures.csv")
dict_to_csv(success_dict, "summary/tc_bcc.csv")

# Format and save plot for fits
plt.xlim(0.0, 0.5)
plt.ylim(0.0, None)
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.savefig("summary/tc_plot.png")
plt.clf()
