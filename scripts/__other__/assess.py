"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import sys; sys.path += [".."]
from cp_sampler.models.cp import Model, STRAIN_RATE
from cp_sampler.helper import round_sf, get_top
from cp_sampler.pole_figure import IPF, get_trajectories
from cp_sampler.plotter import Plotter, save_plot, define_legend

# Constants
MAX_TIME    = 300 # seconds
GRAINS_PATH = "data/grain_data.csv"
LATTICE     = 1.0
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

# Define parameter domains
param_dict_1 = {
    "tau_sat": 700,
    "b":       15,
    "tau_0":   300,
    "gamma_0": round_sf(STRAIN_RATE/3, 4),
    "n":       5.5,
}
param_dict_2 = {
    "tau_sat": 448.56,
    "b":       16.339,
    "tau_0":   527.98,
    "gamma_0": round_sf(STRAIN_RATE/3, 4),
    "n":       3.8835,
}

# Initialises the model and plotter
model = Model(GRAINS_PATH, "fcc", 1.0)
_, top_indexes = get_top(model.get_weights(), TOP_GRAINS)
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
ipf = IPF(model.get_lattice())

# Run both sets of parameters
_, _, results_1 = model.get_results_direct(param_dict_1)
history_1 = model.get_orientation_history()
_, _, results_2 = model.get_results_direct(param_dict_2)
history_2 = model.get_orientation_history()

# Plot the tensile curves
plotter = Plotter(x_label="strain", y_label="stress")
data_dict_1 = {"strain": [round_sf(s[0], 5) for s in results_1["strain"]], "stress": [round_sf(s[0], 5) for s in results_1["stress"]]}
data_dict_2 = {"strain": [round_sf(s[0], 5) for s in results_2["strain"]], "stress": [round_sf(s[0], 5) for s in results_2["stress"]]}
plotter.scat_plot(data_dict_1)
plotter.line_plot(data_dict_2)
define_legend(["darkgray", "green"], ["Experimental", "Calibration"], [7, 1.5], ["scatter", "line"])
save_plot("plot_ss.png")

# Plot the reorientation trajectories
trajectories = get_trajectories(history_1, top_indexes)
ipf.plot_ipf_trajectory(trajectories, direction, {"color": "darkgray"}, scatter=True)
trajectories = get_trajectories(history_2, top_indexes[:5])
ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, {"color": "green"}, scatter=True)
ipf.plot_ipf_trajectory(trajectories, direction, {"color": "green"}, scatter=False)
trajectories = get_trajectories(history_2, top_indexes[5:])
ipf.plot_ipf_trajectory([[t[0]] for t in trajectories], direction, {"color": "red"}, scatter=True)
ipf.plot_ipf_trajectory(trajectories, direction, {"color": "red"}, scatter=False)
define_legend(["darkgray", "green", "red"], ["Experimental", "Calibration", "Validation"], [7, 1.5, 1.5], ["scatter", "line", "line"])
save_plot(f"plot_ipf_{''.join([str(d) for d in direction])}.png")
