# Libraries
import sys; sys.path += ["../.."]
import numpy as np
from cp_sim.helper.general import csv_to_dict, transpose
from cp_sim.io.pole_figure import PF
from cp_sim.io.plotter import save_plot
from cp_sim.simulate import get_lattice

# Paths
SIM_PATH = f"../data/617_s1_sim_all_1p0.csv"
# EXP_PATH = f"../data/617_s1_grains_final.csv"
EXP_PATH = f"../data/617_s1_grains.csv"

# Get experimental final orientations
exp_list = np.loadtxt(EXP_PATH, delimiter=",")
exp_orientations = [list(e[:3]) for e in exp_list]

# Get simulation final orientations
sim_dict = csv_to_dict(SIM_PATH)
sim_orientations = [sim_dict[key] for key in ["0p2_phi_1", "0p2_Phi", "0p2_phi_2"]]
# sim_orientations = [sim_dict[key] for key in ["1p0_phi_1", "1p0_Phi", "1p0_phi_2"]]
sim_orientations = transpose(sim_orientations)

# Initialise plotter
direction = [[1,0,0], [0,1,0], [0,0,1]][0]
pf = PF(get_lattice("fcc"))

# Plot simulated and experimental pole figures
pf.plot_pf(exp_orientations, direction)
save_plot("plot_exp_pf.png")
pf.plot_pf(sim_orientations, direction)
save_plot("plot_sim_pf.png")
