# Libraries
import os, math
import sys; sys.path += ["../.."]
from cp_sampler.helper import csv_to_dict, dict_to_csv, round_sf

# Constants
NUM_GRAINS      = 5
PARAM_NAME_LIST = ["tau_sat", "b", "tau_0", "gamma_0", "n"]

# Initialise success dictionary
success_keys = PARAM_NAME_LIST + [f"g{i+1}_{label}_{pos}" for i in range(NUM_GRAINS)
                                  for label in ["phi_1", "Phi", "phi_2"] for pos in ["start", "end"]]
success_dict = {}
for key in success_keys:
    success_dict[key] = []

# Read all CSV files and iterate through them
results_dir = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/results/cp_neml/20240512 (tensile bcc wide)"
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
    for key in param_dict.keys():
        success_dict[key].append(param_dict[key])

    # Gets the start and end points of the trajectory and store
    for i in range(NUM_GRAINS):
        for label in ["phi_1", "Phi", "phi_2"]:
            for pos in ["start", "end"]:
                value = data_dict[f"{label}_{pos}"][i]
                value = round_sf(value if value > 0 else value + 2*math.pi, 5)
                success_dict[f"g{i+1}_{label}_{pos}"].append(value)

# Write results
dict_to_csv(success_dict, f"phi_bcc.csv")
