"""
 Title:         CP Sampler
 Description:   Sampler for the Crystal Plasticity model
 Author:        Janzen Choi

"""

# Libraries
import threading
import sys; sys.path += [".."]
from cp_sim.models.cp import Model
from cp_sim.pole_figure import get_trajectories
from cp_sim.helper import round_sf, dict_to_csv, get_combinations, csv_to_dict, get_thinned_list
import math

# Paths
EXP_PATH     = f"data/617_s1_exp.csv"
GRAINS_PATH  = f"data/617_s1_grains.csv"
RESULTS_PATH = "results"

# Model constants
MAX_TIME    = 300 # seconds
MAX_STRAIN  = 0.2
STRAIN_RATE = 1e-4
LATTICE     = 1.0
THIN_AMOUNT = 100

def quick_spline(x_list:list, y_list:list, x_value:float) -> float:
    """
    Conducts a quick evaluation using spline interpolation without
    conducting the whole interpolation; assumes that the x_value is
    between min(x_list) and max(x_list) and that x_list is sorted

    Parameters:
    * `x_list`:  The list of x values
    * `y_list`:  The list of y values
    * `x_value`: The x value to evaluate
    
    Returns the evaluated y value
    """
    for i in range(len(x_list)-1):
        if x_list[i] <= x_value and x_value <= x_list[i+1]:
            gradient = (y_list[i+1]-y_list[i])/(x_list[i+1]-x_list[i])
            y_value = gradient*(x_value - x_list[i]) + y_list[i]
            return y_value
    return None

def get_grain_dict(strain_list:list, history:list, grain_ids:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `strain_list`: The list of strain values
    * `history`:     The orientation history
    * `grain_ids`:   The grain indexes to include in the dictionary (starts at 1)
    
    Returns the dictionary of euler-bunge angles (rads)
    """

    # Reformat orientations
    index_list = [grain_id-1 for grain_id in grain_ids] # starts at 0
    trajectories = get_trajectories(history, index_list)

    # Initialise grain dictionary (20%, 40$, .., 100% of max strain)
    strain_intervals = ["0p2", "0p4", "0p6", "0p8", "1p0"]
    grain_dict = {"grain_id": grain_ids}
    for strain_interval in strain_intervals:
        for euler_value in ["phi_1", "Phi", "phi_2"]:
            grain_dict[f"{strain_interval}_{euler_value}"] = []
    
    # Define domain of orientations
    domain = lambda x : round_sf(x, 5) if x>0 else round_sf(x+2*math.pi, 5)

    # Get orientations at different strain intervals
    num_intervals = len(strain_intervals)
    strain_values = [MAX_STRAIN*(i+1)/num_intervals for i in range(num_intervals)]
    for trajectory in trajectories:

        # Get orientations
        phi_1_list = [domain(t[0]) for t in trajectory]
        Phi_list   = [domain(t[1]) for t in trajectory]
        phi_2_list = [domain(t[2]) for t in trajectory]
        
        # Calculate reduced orientations
        reduced_phi_1_list = [quick_spline(strain_list, phi_1_list, strain) for strain in strain_values]
        reduced_Phi_list   = [quick_spline(strain_list, Phi_list, strain)   for strain in strain_values]
        reduced_phi_2_list = [quick_spline(strain_list, phi_2_list, strain) for strain in strain_values]
        
        # Check that the orientations are valid
        if None in reduced_phi_1_list + reduced_Phi_list + reduced_phi_2_list:
            return None

        # Store reduced orientations
        for i, strain_interval in enumerate(strain_intervals):
            grain_dict[f"{strain_interval}_phi_1"].append(reduced_phi_1_list[i])
            grain_dict[f"{strain_interval}_Phi"].append(reduced_Phi_list[i])
            grain_dict[f"{strain_interval}_phi_2"].append(reduced_phi_2_list[i])

    # Return the dictionary
    return grain_dict

# Define parameter domains
index_1 = int(sys.argv[1])
index_2 = int(sys.argv[2])
all_params_dict = {
    "tau_sat": [[200, 400, 800, 1000][index_1]],
    "b":       [0.5, 1, 2, 4, 8, 16],
    "tau_0":   [[200, 400, 800, 1000][index_2]],
    "gamma_0": [round_sf(STRAIN_RATE/3, 4)],
    "n":       [1, 2, 4, 8, 16, 32],
}

# Initialise model
model = Model(
    grains_path = GRAINS_PATH,
    structure   = "fcc",
    lattice_a   = 1.0,
    num_threads = 12,
    strain_rate = STRAIN_RATE,
    max_strain  = MAX_STRAIN,
    youngs      = 211000,
    poissons    = 0.30,
)

# Initialise specific grains
exp_dict = csv_to_dict(EXP_PATH)
grain_ids = [int(key.replace("g","").replace("_phi_1",""))-1 for key in exp_dict.keys() if "phi_1" in key]

# Iterate through the parameters
combinations = get_combinations(all_params_dict)
param_names = list(all_params_dict.keys())
for i, combination in enumerate(combinations):

    # Initialise
    index_str = str(i+1).zfill(3)
    param_dict = dict(zip(param_names, combination))
    results_path = f"{RESULTS_PATH}/{index_1}_{index_2}_{index_str}"

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
    strain_list = get_thinned_list([round_sf(s[0], 5) for s in results["strain"]], THIN_AMOUNT)
    stress_list = get_thinned_list([round_sf(s[0], 5) for s in results["stress"]], THIN_AMOUNT)
    data_dict = {"strain": strain_list, "stress": stress_list}

    # Get grain and stress information
    history = model.get_orientation_history()
    grain_dict = get_grain_dict(strain_list, history, grain_ids)

    # Checks that the orientations were obtained properly
    if grain_dict == None:
        dict_to_csv(param_dict, f"{results_path}_unoriented.csv")
        continue

    # Compile results and write to CSV file
    combined_dict = {**param_dict, **data_dict, **grain_dict}
    dict_to_csv(combined_dict, f"{results_path}.csv")
