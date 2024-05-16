# Libraries
import os, math
import sys; sys.path += ["../.."]
from cp_sim.helper import csv_to_dict, dict_to_csv, round_sf

# Constants
PARAM_NAME_LIST = ["tau_sat", "b", "tau_0", "gamma_0", "n"]
GRAIN_INDEXES   = list(range(54))

# Initialise success dictionary
success_keys = PARAM_NAME_LIST + [f"g{i}_{label}" for i in GRAIN_INDEXES for label in ["phi_1", "Phi", "phi_2"]]
success_dict = {}
for key in success_keys:
    success_dict[key] = []

# Read all CSV files and iterate through them
results_dir = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/20240515 (cp_bcc)"
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
            failed = True
    if failed:
        continue

    # Add parameter information
    for key in param_dict.keys():
        success_dict[key].append(param_dict[key])

    # Gets the start and end points of the trajectory and store
    for i in GRAIN_INDEXES:
        for label in ["phi_1", "Phi", "phi_2"]:
            value = data_dict[label][i]
            value = round_sf(value if value > 0 else value + 2*math.pi, 5)
            success_dict[f"g{i}_{label}"].append(value)

# Write results
dict_to_csv(success_dict, f"results/phi_bcc.csv")
