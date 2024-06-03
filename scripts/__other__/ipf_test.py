# Libraries
import sys; sys.path += ["../.."]
import cp_sim.simulate as sim
from cp_sim.helper.general import csv_to_dict
from cp_sim.io.pole_figure import IPF
from cp_sim.io.plotter import save_plot

# Constants
EXP_PATH  = "../data/617_s1_exp.csv"
GRAIN_IDS = [75, 189, 314, 346, 463]
PHI_LIST  = ["phi_1", "Phi", "phi_2"]
STRAINS   = ["0p2", "0p4", "0p6", "0p8", "1p0"]

# Read experimental data
exp_dict = csv_to_dict(EXP_PATH)
exp_size = len(exp_dict[f"g{GRAIN_IDS[0]}_{PHI_LIST[0]}"])

# Get experimental trajectories
exp_trajectories = []
for grain_id in GRAIN_IDS:
    exp_trajectory = []
    for i in range(exp_size):
        exp_trajectory.append([exp_dict[f"g{grain_id}_{phi}"][i] for phi in PHI_LIST])
    exp_trajectories.append(exp_trajectory)
for et in exp_trajectories:
    print(et[0])

# # Reorient active to passive rotations
# for i in range(len(exp_trajectories)):
#     for j in range(len(exp_trajectories[i])):
#         exp_trajectories[i][j] = sim.reorient(exp_trajectories[i][j])
# print(exp_trajectories[0][:10])

# Initialise plotter
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
ipf = IPF(sim.get_lattice("fcc"))

# Plot experimental trajectories
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
for i, et in enumerate(exp_trajectories):
    ipf.plot_ipf_trajectory([[et[0]]], direction, "scatter", {"color": "darkgray", "s": 8**2})
    ipf.plot_ipf_trajectory([[et[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": GRAIN_IDS[i]})
save_plot("results/ipf_test.png")
