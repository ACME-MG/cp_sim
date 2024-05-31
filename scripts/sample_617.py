"""
 Title:         617 Sampler
 Description:   Runs simulation for 617
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import cp_sim.simulate as sim 
from cp_sim.helper.general import csv_to_dict, dict_to_csv, round_sf, get_thinned_list

# Paths
EXP_PATH     = f"data/617_s1_exp.csv"
GRAINS_PATH  = f"data/617_s1_grains.csv"
RESULTS_PATH = "results"

# Constants
NUM_THREADS = 12
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.2
MAX_TIME    = 300 # seconds
THIN_AMOUNT = 100

# Get grain IDs
exp_dict = csv_to_dict(EXP_PATH)
grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in exp_dict.keys() if "phi_1" in key]

# Get model
model = sim.get_model(
    model_name   = "cp",
    lattice      = sim.get_lattice("fcc"),
    orientations = [sim.get_orientations(GRAINS_PATH)[26]],
    weights      = [sim.get_weights(GRAINS_PATH)[26]],
    num_threads  = NUM_THREADS,
    strain_rate  = STRAIN_RATE,
    max_strain   = MAX_STRAIN,
    youngs       = 211000,
    poissons     = 0.3,
)

# Define parameter domains
index_1 = int(sys.argv[1])
index_2 = int(sys.argv[2])
all_params_dict = {
    "tau_sat": [[100, 200, 400, 800][index_1]],
    "b":       [0.5, 1, 2, 4, 8, 16],
    "tau_0":   [[100, 200, 400, 800][index_2]],
    "gamma_0": [round_sf(STRAIN_RATE/3, 4)],
    "n":       [1, 2, 4, 8, 16, 32],
}

# Iterate through the parameters
param_combinations = sim.get_combinations(all_params_dict)
param_names = list(all_params_dict.keys())
for i, combination in enumerate(param_combinations):

    # Initialise
    index_str = str(i+1).zfill(3)
    param_dict = dict(zip(param_names, combination))
    results_path = f"{RESULTS_PATH}/{index_1}_{index_2}_{index_str}"

    # Runs the model
    status = model.run(param_dict, max_time=MAX_TIME)
    if status != "success":
        dict_to_csv(param_dict, f"{results_path}_{status}.csv")
        continue
    _, pc_model, results = model.get_output()

    # Process results
    strain_list = get_thinned_list([round_sf(s[0], 5) for s in results["strain"]], THIN_AMOUNT)
    stress_list = get_thinned_list([round_sf(s[0], 5) for s in results["stress"]], THIN_AMOUNT)
    data_dict   = {"strain": strain_list, "stress": stress_list}
    history     = sim.get_orientation_history(pc_model, results)
    grain_dict  = sim.get_grain_dict(strain_list, history, grain_ids, MAX_STRAIN)

    # Check and save results
    if grain_dict == None:
        dict_to_csv(param_dict, f"{results_path}_unoriented.csv")
        continue
    dict_to_csv({**param_dict, **data_dict, **grain_dict}, f"{results_path}.csv")
