"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import sys; sys.path += [".."]
from cp_sampler.models.cp import Model, STRAIN_RATE
from cp_sampler.helper import round_sf, get_sorted, csv_to_dict
from cp_sampler.pole_figure import IPF, get_trajectories
from cp_sampler.plotter import Plotter, save_plot, define_legend

# Constants
MAX_TIME     = 300 # seconds
GRAINS_PATH  = "data/grain_p91.csv"
CALIB_GRAINS = 5
VALID_GRAINS = 10

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

# Initialises the model and plotter
model = Model(GRAINS_PATH, "bcc", 1.0)
_, sorted_indexes = get_sorted(model.get_weights())
direction = [[1,0,0], [0,1,0], [0,0,1]][1]
ipf = IPF(model.get_lattice())

# Gets the experimental data
exp_path = "data/tensile_p91.csv"
exp_dict = csv_to_dict(exp_path)
num_grains = len([field for field in exp_dict.keys() if "phi_1" in field])

# Gets experimental history
exp_history = [[] for _ in range(2)] # start and end
for i in range(1,num_grains+1):
    phi_1 = exp_dict[f"g{i}_phi_1"]
    Phi   = exp_dict[f"g{i}_Phi"]
    phi_2 = exp_dict[f"g{i}_phi_2"]
    exp_history[0].append([phi_1[0], Phi[0], phi_2[0]])
    exp_history[1].append([phi_1[-1], Phi[-1], phi_2[-1]])

# Get simulated results
param_dict = {
    "tau_sat": 50.039,
    "b":       25.486,
    "tau_0":   182.49,
    "gamma_0": round_sf(STRAIN_RATE/3, 4),
    "n":       12.054,
}
_, _, sim_results = model.get_results_direct(param_dict)
sim_history = model.get_orientation_history()

# Plot the tensile curves
plotter = Plotter(x_label="strain", y_label="stress")
plotter.scat_plot(exp_dict)
data_dict = {"strain": [round_sf(s[0], 5) for s in sim_results["strain"]], "stress": [round_sf(s[0], 5) for s in sim_results["stress"]]}
plotter.line_plot(data_dict)
define_legend(["darkgray", "green"], ["Experimental", "Calibration"], [7, 1.5], ["scatter", "line"])
save_plot("plot_ss.png")

# Plot the experimental reorientation trajectories
exp_trajectories = get_trajectories(exp_history, list(range(VALID_GRAINS)))
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
ipf.plot_ipf_trajectory([[et[0]] for et in exp_trajectories], direction, "scatter", {"color": "darkgray", "s": 8**2})

# Plot the calibration reorientation trajectories
sim_trajectories = get_trajectories(sim_history, sorted_indexes[:CALIB_GRAINS])
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "green", "linewidth": 1})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2})

# Plot the validaiton reorientation trajectories
sim_trajectories = get_trajectories(sim_history, sorted_indexes[CALIB_GRAINS:VALID_GRAINS])
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "red", "linewidth": 1})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "red", "head_width": 0.0075, "head_length": 0.0075*1.5})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "red", "s": 6**2})

# Format and save
define_legend(["darkgray", "green", "red"], ["Experimental", "Calibration", "Validation"], [7, 1.5, 1.5], ["scatter", "line", "line"])
save_plot(f"plot_ipf_{''.join([str(d) for d in direction])}.png")
