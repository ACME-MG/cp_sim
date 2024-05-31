# Libraries
import sys; sys.path += ["../.."]
from cp_sim.helper.general import csv_to_dict
from cp_sim.io.pole_figure import PF
from cp_sim.io.plotter import save_plot, define_legend
from cp_sim.simulate import get_lattice

# Paths
SUM_ID    = "617_s1"
SIM_PATH  = f"results/{SUM_ID}_phi.csv"
EXP_PATH  = f"../data/{SUM_ID}_exp.csv"

# Parameters
GRAIN_IDS = [75, 189, 314, 346, 463]
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
pf = PF(get_lattice("fcc"))

# Plot simulated and experimental pole figures
pf.plot_pf([et[-1] for et in exp_trajectories], direction)
pf.plot_pf([st[-1] for st in sim_trajectories], direction)
save_plot("plot_pf.png")
