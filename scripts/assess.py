"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import sys; sys.path += [".."]
from cp_sampler.models.cp import Model
from cp_sampler.helper import round_sf, csv_to_dict
from cp_sampler.pole_figure import IPF, get_trajectories
from cp_sampler.plotter import Plotter, save_plot, define_legend

# Constants
MAX_TIME      = 300 # seconds
EXP_PATH      = "data/tensile_p91.csv"
GRAINS_PATH   = "data/grain_p91.csv"
MAPPING_PATH  = "data/mapping_p91.csv"
CALIB_INDEXES = [0,1,2] # starts at 0
VALID_INDEXES = []  # starts at 0

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
        grain_dict[f"g{i}_phi_1"] = []
        grain_dict[f"g{i}_Phi"]   = []
        grain_dict[f"g{i}_phi_2"] = []

        # Get the trajectory of each grain throughout history
        euler_list = [[], [], []]
        for state in history:
            orientations = pc_model.orientations(state)
            euler = list(orientations[indexes[i]].to_euler(angle_type="radians", convention="bunge"))
            for j in range(len(euler_list)):
                euler_value = euler[j] if euler[j] > 0 else euler[j]+2*np.pi
                euler_list[j].append(euler_value)

        # Store the trajectories
        grain_dict[f"g{i}_phi_1"] = euler_list[0]
        grain_dict[f"g{i}_Phi"]   = euler_list[1]
        grain_dict[f"g{i}_phi_2"] = euler_list[2]
    
    # Return dictionary
    return grain_dict

# Gets the experimental data
exp_dict = csv_to_dict(EXP_PATH)
num_grains = len([field for field in exp_dict.keys() if "phi_1" in field])
strain_rate = round_sf(max(exp_dict["strain"])/max(exp_dict["time"]), 5)

# Initialises model
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

# Initialise plotter
direction = [[1,0,0], [0,1,0], [0,0,1]][1]
ipf = IPF(model.get_lattice())

# Initialise indexes for grains to capture
map_dict = csv_to_dict(MAPPING_PATH)
start_indexes = [int(si)-1 for si in list(map_dict["start"])]
end_indexes = [int(ei)-1 for ei in list(map_dict["end"])]

# Gets experimental history
exp_history = [[] for _ in range(2)] # start and end
for i in range(num_grains):
    phi_1 = exp_dict[f"g{i}_phi_1"]
    Phi   = exp_dict[f"g{i}_Phi"]
    phi_2 = exp_dict[f"g{i}_phi_2"]
    exp_history[0].append([phi_1[0], Phi[0], phi_2[0]])
    exp_history[1].append([phi_1[-1], Phi[-1], phi_2[-1]])

# Get simulated results
param_names = ["tau_sat", "b", "tau_0", "gamma_0", "n"]
param_str = """
202.56	5	200.94	3.33E-05	8
"""
param_list = [float(p) for p in param_str.split("\t")]
param_dict = dict(zip(param_names, param_list))
_, _, sim_results = model.get_results_direct(param_dict, try_run=False)
sim_history = model.get_orientation_history()

# Plot the tensile curves
plotter = Plotter(x_label="strain", y_label="stress")
plotter.scat_plot(exp_dict)
plotter.line_plot({
    "strain": [round_sf(s[0], 5) for s in sim_results["strain"]],
    "stress": [round_sf(s[0], 5) for s in sim_results["stress"]]
})
define_legend(["darkgray", "green"], ["Experimental", "Calibration"], [7, 1.5], ["scatter", "line"])
save_plot("plot_ss.png")

# Plot the experimental reorientation trajectories
exp_indexes = CALIB_INDEXES+VALID_INDEXES
exp_trajectories = get_trajectories(exp_history, exp_indexes)
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
for i, et in enumerate(exp_trajectories):
    ipf.plot_ipf_trajectory([[et[0]]], direction, "scatter", {"color": "darkgray", "s": 8**2})
    ipf.plot_ipf_trajectory([[et[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": start_indexes[exp_indexes[i]]+1})
    ipf.plot_ipf_trajectory([[et[-1]]], direction, "text", {"color": "black", "fontsize": 8, "s": end_indexes[exp_indexes[i]]+1})

# Plot the calibration reorientation trajectories
sim_trajectories = get_trajectories(sim_history, [start_indexes[i] for i in CALIB_INDEXES])
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "green", "linewidth": 1})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2})

# Plot the validation reorientation trajectories
sim_trajectories = get_trajectories(sim_history, [start_indexes[i] for i in VALID_INDEXES])
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "red", "linewidth": 1})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "red", "head_width": 0.0075, "head_length": 0.0075*1.5})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "red", "s": 6**2})

# Format and save
define_legend(["darkgray", "green", "red"], ["Experimental", "Calibration", "Validation"], [7, 1.5, 1.5], ["scatter", "line", "line"])
save_plot(f"plot_ipf_{''.join([str(d) for d in direction])}.png")
