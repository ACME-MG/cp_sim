"""
 Title:         CP Runner
 Description:   Runs the Crystal Plasticity model once
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from neml.cp import crystallography
import sys; sys.path += [".."]
from cp_sim.helper import round_sf, csv_to_dict
from cp_sim.pole_figure import IPF, get_trajectories
from cp_sim.plotter import Plotter, save_plot, define_legend

# Constants
MAX_TIME = 300 # seconds
EXP_PATH = f"data/reorientation.csv"
INCLUDE = None

# Gets the experimental data
exp_dict = csv_to_dict(EXP_PATH)
mappable_grain_ids = [int(field.replace("g","").replace("_phi_1","")) for field in exp_dict.keys() if "phi_1" in field]
exit()

num_grains = len()
strain_rate = round_sf(max(exp_dict["strain"])/max(exp_dict["time"]), 5)

# Initialise plotter
lattice = crystallography.CubicLattice(1.0)
lattice.add_slip_system([1,1,0], [1,1,1])
# lattice.add_slip_system([1,1,1], [1,1,0])
# lattice.add_slip_system([1,1,1], [1,2,3])
# lattice.add_slip_system([1,1,1], [1,1,2])
ipf = IPF(lattice)

# Initialise indexes for grains to capture
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
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
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "green", "linewidth": 1, "zorder": 3})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})

# Plot the validation reorientation trajectories
sim_trajectories = get_trajectories(sim_history, [start_indexes[i] for i in VALID_INDEXES])
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "red", "linewidth": 1, "zorder": 3})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "red", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "red", "s": 6**2, "zorder": 3})

# Format and save
define_legend(["darkgray", "green", "red"], ["Experimental", "Calibration", "Validation"], [7, 1.5, 1.5], ["scatter", "line", "line"])
save_plot(f"plot_ipf_{''.join([str(d) for d in direction])}.png")
