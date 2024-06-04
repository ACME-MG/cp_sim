# Libraries
import sys; sys.path += ["../.."]
import cp_sim.simulate as sim
from cp_sim.helper.general import csv_to_dict
from cp_sim.io.pole_figure import IPF
from cp_sim.io.plotter import save_plot, define_legend

# Paths
SUM_ID    = "p91_s3"
SIM_PATH  = f"results/{SUM_ID}_phi.csv"
EXP_PATH  = f"../data/{SUM_ID}_exp.csv"

# Parameters
# GRAIN_IDS = [75, 189, 314, 346, 463] # 56, 346, 463, 568*, 650
GRAIN_IDS = [43, 44, 45, 47, 48] # 4, 29, 33
PHI_LIST  = ["phi_1", "Phi", "phi_2"]
STRAINS   = ["0p2", "0p4", "0p6", "0p8", "1p0"]

# Read experimental and simulated data
exp_dict = csv_to_dict(EXP_PATH)
exp_size = len(exp_dict[f"g{GRAIN_IDS[0]}_{PHI_LIST[0]}"])
sim_dict = csv_to_dict(SIM_PATH)
num_sims = len(sim_dict[f"g{GRAIN_IDS[0]}_{STRAINS[0]}_{PHI_LIST[0]}"])

# Get experimental trajectories
exp_trajectories = []
for grain_id in GRAIN_IDS:
    exp_trajectory = []
    for i in range(exp_size):
        exp_trajectory.append([exp_dict[f"g{grain_id}_{phi}"][i] for phi in PHI_LIST])
    exp_trajectories.append(exp_trajectory)

# Get simulation trajectories
sim_trajectories = []
for grain_id in GRAIN_IDS:
    for i in range(num_sims):
        sim_trajectory = []
        for strain in STRAINS:
            sim_trajectory.append([sim_dict[f"g{grain_id}_{strain}_{phi}"][i] for phi in PHI_LIST]) 
        sim_trajectories.append(sim_trajectory)

# Initialise plotter
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
ipf = IPF(sim.get_lattice("fcc"))

# Plot experimental trajectories
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
for i, et in enumerate(exp_trajectories):
    ipf.plot_ipf_trajectory([[et[0]]], direction, "scatter", {"color": "darkgray", "s": 8**2})
    ipf.plot_ipf_trajectory([[et[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": GRAIN_IDS[i]})

# Plot the simulated trajectories
ipf.plot_ipf_trajectory(sim_trajectories, direction, "scatter", {"color": "green", "s": 6**2})

# Save
define_legend(["darkgray", "green"], ["Experimental", "Simulation"], [2, 6], ["line", "scatter"])
save_plot("results/plot_ipf.png")
