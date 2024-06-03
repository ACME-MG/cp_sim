"""
 Title:         Summarise
 Description:   Summarises tensile simulation runs into small CSV files
 Author:        Janzen Choi

"""

# Libraries
import os, math
import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sim.helper.general import csv_to_dict, dict_to_csv, round_sf
from cp_sim.helper.interpolator import Interpolator

# Constants
RESULTS_DIR     = "results"
SUMMARY_ID      = "617_s1"
SAMPLE_PATH     = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/2024-06-03 (617_s1_activated)"
PARAM_NAME_LIST = ["tau_sat", "b", "tau_0", "gamma_0", "n"]
STRAIN_HEADER   = "strain"
STRESS_HEADER   = "stress"
NUM_STRAINS     = 30

def append_to_dict(list_dict:dict, key:str, value:float) -> dict:
    """
    Appends a value to a dictionary;
    if key doesn't exist, initialises a list for it

    Parameters:
    * `list_dict`: The dictionary
    * `key`:      The key to add the value to
    * `value`:    The value to append
    """
    if not key in list_dict.keys():
        list_dict[key] = []
    list_dict[key].append(value)
    return list_dict

def append_list_dict(list_dict:dict, value_dict:dict) -> dict:
    """
    Appends values from a dictionary to lists in another dictionary;
    if key doesn't exist, initialises a list for it
    
    Parameters:
    * `list_dict`:  The dictionary of lists
    * `value_dict`: The dictionary of values

    Returns the new dictionary
    """
    for key in value_dict.keys():
        append_to_dict(list_dict, key, value_dict[key])
    return list_dict

# Initialise summary dictionaries
sum_tc_dict = {}
sum_phi_dict = {}
sum_fail_dict = {}

# Initialise plot
plt.figure(figsize=(6,6))

# Iterate through CSV files
csv_file_list = [file for file in os.listdir(SAMPLE_PATH) if file.endswith(".csv")]
for csv_file in csv_file_list:
    
    # Convert csv file to dictionary
    data_dict = csv_to_dict(f"{SAMPLE_PATH}/{csv_file}")

    # Get parameter informationn
    param_dict = {}
    for param_name in PARAM_NAME_LIST:
        param_dict[param_name] = data_dict[param_name]

    # Check whether the simulation failed or timed out
    failed = False
    for keyword in ["failed", "timeout", "unoriented"]:
        if keyword in csv_file:
            sum_fail_dict = append_list_dict(sum_fail_dict, param_dict)
            sum_fail_dict = append_to_dict(sum_fail_dict, "reason", keyword)
            failed = True
    if failed:
        continue

    # Check whether the number of datapoints are sufficient
    if len(data_dict["strain"]) < 5:
        sum_fail_dict = append_list_dict(sum_fail_dict, param_dict) 
        sum_fail_dict = append_to_dict(sum_fail_dict, "reason", "insufficient")
        continue

    # Add parameter information
    sum_tc_dict = append_list_dict(sum_tc_dict, param_dict)
    sum_phi_dict = append_list_dict(sum_phi_dict, param_dict)

    # Add to tensile curve dictionary
    interpolator = Interpolator(data_dict[STRAIN_HEADER], data_dict[STRESS_HEADER], resolution=100)
    x_end  = max(data_dict[STRAIN_HEADER])
    x_list = [x_end/(NUM_STRAINS-1)*j for j in range(NUM_STRAINS)]
    y_list = interpolator.evaluate(x_list)
    sum_tc_dict = append_to_dict(sum_tc_dict, "x_end", round_sf(x_end, 5))
    for i in range(NUM_STRAINS):
        sum_tc_dict = append_to_dict(sum_tc_dict, f"y_{i+1}", round_sf(y_list[i], 5))

    # Add to orientation dictionary
    for i, grain_id in enumerate(data_dict["grain_id"]):
        for strain_interval in ["0p2", "0p4", "0p6", "0p8", "1p0"]:
            for phi in ["phi_1", "Phi", "phi_2"]:
                label = f"{strain_interval}_{phi}"
                value = data_dict[label][i]
                value = round_sf(value if value > 0 else value + 2*math.pi, 5)
                sum_phi_dict = append_to_dict(sum_phi_dict, f"g{int(grain_id)}_{label}", value)

    # Plot tensile curves
    plt.scatter(data_dict[STRAIN_HEADER], data_dict[STRESS_HEADER], color="silver")
    plt.plot(x_list, y_list)

# Write results
dict_to_csv(sum_fail_dict, f"{RESULTS_DIR}/{SUMMARY_ID}_fail.csv")
dict_to_csv(sum_tc_dict,   f"{RESULTS_DIR}/{SUMMARY_ID}_tc.csv")
dict_to_csv(sum_phi_dict,  f"{RESULTS_DIR}/{SUMMARY_ID}_phi.csv")

# Format and save plot for tensile curves
plt.xlim(0.0, 0.3)
plt.ylim(0.0, None)
plt.xlabel("Strain (mm/mm)", fontsize=12)
plt.ylabel("Stress (MPa)", fontsize=12)
plt.savefig(f"{RESULTS_DIR}/tc_plot.png")
plt.clf()