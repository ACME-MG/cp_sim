"""
 Title:         617 Runner
 Description:   Runs simulation for 617
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import cp_sim.simulate as sim 
from cp_sim.helper.general import csv_to_dict, round_sf, get_thinned_list, transpose
from cp_sim.io.pole_figure import IPF
from cp_sim.io.plotter import Plotter, define_legend, save_plot

# Paths
EXP_PATH     = f"data/617_s1_exp.csv"
GRAINS_PATH  = f"data/617_s1_grains.csv"
RESULTS_PATH = "results"

# Constants
NUM_THREADS = 12
STRAIN_RATE = 1e-4
MAX_STRAIN  = 0.2
MAX_TIME    = 3000 # seconds
THIN_AMOUNT = 100

# Get grain IDs
exp_dict = csv_to_dict(EXP_PATH) # passive
grain_ids = [56, 346, 463, 568, 650] # [75, 189, 314, 338, 463]
# grain_ids = [int(key.replace("g","").replace("_phi_1","")) for key in exp_dict.keys() if "phi_1" in key]
# grain_ids = [i+1 for i in range(len(sim.get_orientations(GRAINS_PATH)))]

# Get model
model = sim.get_model(
    model_name   = "cp",
    lattice      = sim.get_lattice("fcc"),
    orientations = sim.get_orientations(GRAINS_PATH, angle_type="radians", is_passive=True), # passive to active
    weights      = sim.get_weights(GRAINS_PATH),
    num_threads  = NUM_THREADS,
    strain_rate  = STRAIN_RATE,
    max_strain   = MAX_STRAIN,
    youngs       = 211000,
    poissons     = 0.3,
)

# Run the model
param_dict = {"tau_sat": 108.35, "b": 0.5840, "tau_0": 120.21, "gamma_0": round_sf(STRAIN_RATE/3, 4), "n": 2.5832}
status = model.run(param_dict, max_time=MAX_TIME)
_, pc_model, results = model.get_output() # active to passive

# Get simulation results
sim_strain_list  = [round_sf(s[0], 5) for s in results["strain"]]
sim_stress_list  = [round_sf(s[0], 5) for s in results["stress"]]
sim_history      = sim.get_orientation_history(pc_model, results, inverse=False) # passive
sim_trajectories = sim.get_trajectories(sim_history, [grain_id-1 for grain_id in grain_ids])
sim_dict = {"strain": get_thinned_list(sim_strain_list, THIN_AMOUNT), "stress": get_thinned_list(sim_stress_list, THIN_AMOUNT)}

# Plot stress-strain curve
plotter = Plotter("strain", "stress")
plotter.scat_plot(exp_dict)
plotter.line_plot(sim_dict)
save_plot("plot_ss.png")

# Initialise IPF plotter
direction = [[1,0,0], [1,1,0], [1,1,1]][0]
ipf = IPF(sim.get_lattice("fcc"))

# Plot experimental IPF
exp_trajectories = [transpose([exp_dict[f"g{grain_id}_{phi}"] for phi in ["phi_1", "Phi", "phi_2"]]) for grain_id in grain_ids]
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
for grain_id, et in zip(grain_ids, exp_trajectories):
    ipf.plot_ipf_trajectory([[et[0]]], direction, "scatter", {"color": "darkgray", "s": 8**2})
    ipf.plot_ipf_trajectory([[et[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": grain_id})

# Plot simulated IPF
ipf.plot_ipf_trajectory(sim_trajectories, direction, "plot", {"color": "green", "linewidth": 1, "zorder": 3})
ipf.plot_ipf_trajectory(sim_trajectories, direction, "arrow", {"color": "green", "head_width": 0.0075, "head_length": 0.0075*1.5, "zorder": 3})
ipf.plot_ipf_trajectory([[st[0]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2, "zorder": 3})

# Format and save
define_legend(["darkgray", "green"], ["Experimental", "Calibration"], [7, 1.5], ["scatter", "line"])
save_plot("plot_ipf.png")

