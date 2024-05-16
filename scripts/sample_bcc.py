"""
 Title:         CP Sampler
 Description:   Sampler for the Crystal Plasticity model
 Author:        Janzen Choi

"""

# Libraries
import threading
import sys; sys.path += [".."]
from cp_sim.models.cp import Model
from cp_sim.helper import round_sf, dict_to_csv, get_combinations, csv_to_dict
import math

# Constants
MAX_TIME      = 300 # seconds
SAMPLE_INDEX  = 0
EXP_PATH      = f"data/tensile_p91_{SAMPLE_INDEX}.csv"
GRAINS_PATH   = f"data/grain_p91_{SAMPLE_INDEX}.csv"
MAPPING_PATH  = f"data/mapping_p91_{SAMPLE_INDEX}.csv"
LATTICE       = 1.0

def get_grain_dict(history:list, indexes:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `history`: The orientation history
    * `indexes`: The grain indexes to include in the dictionary
    
    Returns the dictionary of euler-bunge angles (rads)
    """
    grain_dict = {"phi_1": [], "Phi": [], "phi_2": []}
    domain = lambda x : x if x>0 else x+2*math.pi
    for i in indexes:
        for j, key in enumerate(grain_dict.keys()):
            grain_dict[key].append(domain(history[-1][i][j]))
    return grain_dict

# Gets the experimental data
exp_dict = csv_to_dict(EXP_PATH)
num_grains = len([field for field in exp_dict.keys() if "phi_1" in field])
strain_rate = round_sf(max(exp_dict["strain"])/max(exp_dict["time"]), 5)

# Define parameter domains
index_1 = int(sys.argv[1])
index_2 = int(sys.argv[2])
all_params_dict = {
    "tau_sat": [[100, 200, 400, 800][index_1]],
    "b":       [0.5, 1, 2, 4, 8, 16],
    "tau_0":   [[100, 200, 400, 800][index_2]],
    "gamma_0": [round_sf(strain_rate/3, 4)],
    "n":       [1, 2, 4, 8, 16, 32],
}

# Initialise model
model = Model(
    grains_path = GRAINS_PATH,
    structure   = "bcc",
    lattice_a   = 1.0,
    num_threads = 12,
    strain_rate = strain_rate,
    max_strain  = 0.30,
    youngs      = 190000,
    poissons    = 0.28,
)

# Initialise specific grains
map_dict = csv_to_dict(MAPPING_PATH)
sorted_indexes = list(map_dict["start"])
sorted_indexes = [int(si)-1 for si in sorted_indexes]

# Iterate through the parameters
combinations = get_combinations(all_params_dict)
param_names = list(all_params_dict.keys())
for i, combination in enumerate(combinations):

    # Initialise
    index_str = str(i+1).zfill(3)
    param_dict = dict(zip(param_names, combination))
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
    data_dict = {"strain": strain_list, "stress": stress_list}

    # Get grain and stress information
    history = model.get_orientation_history()
    grain_dict = get_grain_dict(history, sorted_indexes)

    # Compile results and write to CSV file
    combined_dict = {**param_dict, **data_dict, **grain_dict}
    dict_to_csv(combined_dict, f"{results_path}.csv")
