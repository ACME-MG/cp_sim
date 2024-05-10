"""
 Title:         CP Sampler
 Description:   Sampler for the Crystal Plasticity model
 Author:        Janzen Choi

"""

# Libraries
import threading
import sys; sys.path += [".."]
from cp_sampler.models.cp import Model, STRAIN_RATE
from cp_sampler.helper import round_sf, dict_to_csv, get_combinations, csv_to_dict
import math

# Constants
MAX_TIME     = 300 # seconds
GRAINS_PATH  = "data/grain_p91.csv"
MAPPING_PATH = "data/mapping_p91.csv"
LATTICE      = 1.0
TOP_GRAINS   = 10

def get_grain_dict(history:list, indexes:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `history`: The orientation history
    * `indexes`: The grain indexes to include in the dictionary
    
    Returns the dictionary of euler-bunge angles (rads)
    """
    grain_dict = {"phi_1_start": [], "phi_1_end": [], "Phi_start": [],
                  "Phi_end": [], "phi_2_start": [], "phi_2_end": []}
    domain = lambda x_list : [x if x>0 else x+2*math.pi for x in x_list]
    for i in indexes:
        for j, key in enumerate(["phi_1_start", "Phi_start", "phi_2_start"]):
            grain_dict[key].append(domain(history[0][i][j]))
        for j, key in enumerate(["phi_1_end", "Phi_end", "phi_2_end"]):
            grain_dict[key].append(domain(history[-1][i][j]))
    return grain_dict

# Define parameter domains
index_1 = int(sys.argv[1])
index_2 = int(sys.argv[2])
all_params_dict = {
    "tau_sat": [[50, 100, 200, 400, 800, 1600][index_2]],
    "b":       [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6],
    "tau_0":   [50, 100, 200, 400, 800],
    "gamma_0": [round_sf(STRAIN_RATE/3, 4)],
    "n":       [[2, 4, 8, 16][index_1]],
}

# Initialise model and top grain weights
model = Model(GRAINS_PATH, "bcc", 1.0)
map_dict = csv_to_dict(MAPPING_PATH)
sorted_indexes = list(map_dict["start"])
sorted_indexes = [int(si)-1 for si in sorted_indexes]

# Iterate through the parameters
combinations = get_combinations(all_params_dict)
param_names = list(all_params_dict.keys())
for i in range(len(combinations)):

    # Initialise
    index_str = str(i+1).zfill(3)
    param_dict = dict(zip(param_names, combinations[i]))
    results_path = f"results/{index_1}_{index_2}_{index_str}"

    # Prepare the thread for the function
    model.define_params(**param_dict)
    thread = threading.Thread(target=model.run_cp)
    thread.start()
    thread.join(timeout=MAX_TIME)

    # Runs the model for the parameter set
    if thread.is_alive():
        dict_to_csv(param_dict, f"{results_path}_timeout.csv")
        continue
    model_output = model.get_results()
    if model_output == None:
        dict_to_csv(param_dict, f"{results_path}_failed.csv")
        continue
    sc_model, pc_model, results = model_output

    # Get tensile curve
    strain_list = [round_sf(s[0], 5) for s in results["strain"]]
    stress_list = [round_sf(s[0], 5) for s in results["stress"]]
    data_dict = {
        "strain": strain_list,
        "stress": stress_list,
    }

    # Get grain and stress information
    history = model.get_orientation_history()
    grain_dict = get_grain_dict(history, sorted_indexes[:TOP_GRAINS])

    # Compile results and write to CSV file
    combined_dict = {**param_dict, **data_dict, **grain_dict}
    dict_to_csv(combined_dict, f"{results_path}.csv")
