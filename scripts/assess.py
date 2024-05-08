"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import sys; sys.path += [".."]
from cp_sampler.model import Model, STRAIN_RATE
from cp_sampler.helper import round_sf, get_top, csv_to_dict
from cp_sampler.pole_figure import IPF, get_trajectories
from cp_sampler.plotter import Plotter, save_plot, define_legend

# Constants
MAX_TIME    = 300 # seconds
GRAINS_PATH = "data/grain_p91.csv"
TOP_GRAINS  = 20

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

# Gets the experimental data
exp_path = "/mnt/c/Users/Janzen/OneDrive - UNSW/PhD/data/20240506 (ondrej_P91)/tensile_p91.csv"
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

# Define parameter domains
param_dict_1 = {
    "tau_sat": 700,
    "b":       15,
    "tau_0":   300,
    "gamma_0": round_sf(STRAIN_RATE/3, 4),
    "n":       5.5,
}

# Initialises the model and plotter
model = Model(GRAINS_PATH, "bcc", 1.3)
_, top_indexes = get_top(model.get_weights(), TOP_GRAINS)
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
ipf = IPF(model.get_lattice())

# # Run both sets of parameters
# _, _, results_1 = model.get_results_direct(param_dict_1)
# history_1 = model.get_orientation_history()

# # Plot the tensile curves
# plotter = Plotter(x_label="strain", y_label="stress")
# data_dict_1 = {"strain": [round_sf(s[0], 5) for s in results_1["strain"]], "stress": [round_sf(s[0], 5) for s in results_1["stress"]]}
# plotter.scat_plot(data_dict_1)
# define_legend(["darkgray", "green"], ["Experimental", "Calibration"], [7, 1.5], ["scatter", "line"])
# save_plot("plot_ss.png")


# Plot the experimental reorientation trajectories
exp_trajectories = get_trajectories(exp_history, list(range(10)))
ipf.plot_ipf_trajectory(exp_trajectories, direction, {"color": "darkgray", "linewidth": 3}, scatter=False)
exp_trajectories = [[et[0]] for et in exp_trajectories]
ipf.plot_ipf_trajectory(exp_trajectories, direction, {"color": "darkgray", "s": 8**2}, scatter=True)

# Plot the simulated reorientation trajectories
# trajectories = get_trajectories(history_2, top_indexes[:5])
# ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, {"color": "green"}, scatter=True)
# ipf.plot_ipf_trajectory(trajectories, direction, {"color": "green"}, scatter=False)
# trajectories = get_trajectories(history_2, top_indexes[5:])
# ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, {"color": "red"}, scatter=True)
# ipf.plot_ipf_trajectory(trajectories, direction, {"color": "red"}, scatter=False)
# define_legend(["darkgray", "green", "red"], ["Experimental", "Calibration", "Validation"], [7, 1.5, 1.5], ["scatter", "line", "line"])
save_plot(f"plot_ipf_{''.join([str(d) for d in direction])}.png")
