"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import sys; sys.path += [".."]
from cp_sampler.models.cp import Model, STRAIN_RATE
from cp_sampler.helper import round_sf, dict_to_csv, get_top

# Constants
MAX_TIME    = 300 # seconds
GRAINS_PATH = "data/grain_data.csv"
LATTICE     = 1.0
TOP_GRAINS  = 10

def get_grain_dict(pc_model:dict, history:dict, indexes:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `pc_model`: The polycrystal model
    * `history`:  The history of the model simulation
    * `indexes`:  The grain indexes to include in the dictionary
    
    Returns the dictionary of euler-bunge angles (rads)
    """
    
    # Iterate through each grain
    grain_dict = {}
    for i in range(len(indexes)):
        
        # Initialise
        grain_dict[f"g{i+1}_phi_1"] = []
        grain_dict[f"g{i+1}_Phi"]   = []
        grain_dict[f"g{i+1}_phi_2"] = []

        # Get the trajectory of each grain throughout history
        euler_list = [[], [], []]
        for state in history:
            orientations = pc_model.orientations(state)
            euler = list(orientations[indexes[i]].to_euler(angle_type="radians", convention="bunge"))
            for j in range(len(euler_list)):
                euler_value = euler[j] if euler[j] > 0 else euler[j]+2*np.pi
                euler_list[j].append(euler_value)

        # Store the trajectories
        grain_dict[f"g{i+1}_phi_1"] = euler_list[0]
        grain_dict[f"g{i+1}_Phi"]   = euler_list[1]
        grain_dict[f"g{i+1}_phi_2"] = euler_list[2]
    
    # Return dictionary
    return grain_dict

# Define parameter domains
param_dict = {
    "tau_sat": 448.56,
    "b":       16.339,
    "tau_0":   527.98,
    "gamma_0": round_sf(STRAIN_RATE/3, 4),
    "n":       3.8835,
}
# all_params_dict = {
#     "tau_sat": 700,
#     "b":       15,
#     "tau_0":   300,
#     "gamma_0": round_sf(STRAIN_RATE/3, 4),
#     "n":       5.5,
# }

# Define and run the model
model = Model(GRAINS_PATH, "fcc", 1.0)
model.define_params(**param_dict)
model.run_cp()
sc_model, pc_model, results = model.get_results()

# Get tensile curve
strain_list = [round_sf(s[0], 5) for s in results["strain"]]
stress_list = [round_sf(s[0], 5) for s in results["stress"]]
data_dict = {
    "strain": strain_list,
    "stress": stress_list,
}

# Get grain and stress information
history = np.array(results["history"])
_, top_indexes = get_top(model.get_weights(), TOP_GRAINS)
grain_dict = get_grain_dict(pc_model, history, top_indexes)

# Compile results and write to CSV file
combined_dict = {**param_dict, **data_dict, **grain_dict}
dict_to_csv(combined_dict, "results/once.csv")
