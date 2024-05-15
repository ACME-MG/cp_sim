# Libraries
import sys; sys.path += ["../.."]
from cp_sampler.models.cp import Model
from cp_sampler.helper import csv_to_dict
from cp_sampler.pole_figure import IPF, get_trajectories
from cp_sampler.helper import csv_to_dict
from cp_sampler.plotter import save_plot, define_legend

# Constants
GRAIN_INDEXES = [0, 1, 2]
SIM_PATH      = "summary/phi_bcc.csv"
EXP_PATH      = "../data/tensile_p91.csv"
GRAINS_PATH   = "../data/grain_p91.csv"
MAPPING_PATH  = "../data/mapping_p91.csv"

# Initialises model
model = Model(
    grains_path = GRAINS_PATH,
    structure   = "bcc",
    lattice_a   = 1.0,
    num_threads = 12,
    strain_rate = 1.0e-4,
    max_strain  = 0.30,
    youngs      = 190000,
    poissons    = 0.28,
)

# Initialise plotter
direction = [[1,0,0], [0,1,0], [0,0,1]][1]
ipf = IPF(model.get_lattice())

# Initialise grain mapping
map_dict = csv_to_dict(MAPPING_PATH)
start_indexes = [int(si)-1 for si in list(map_dict["start"])]
end_indexes = [int(ei)-1 for ei in list(map_dict["end"])]

# Read experimental data
exp_dict = csv_to_dict(EXP_PATH)
num_grains = len([field for field in exp_dict.keys() if "phi_1" in field])
exp_history = [[] for _ in range(2)] # start and end
for i in range(num_grains):
    phi_1 = exp_dict[f"g{i}_phi_1"]
    Phi   = exp_dict[f"g{i}_Phi"]
    phi_2 = exp_dict[f"g{i}_phi_2"]
    exp_history[0].append([phi_1[0], Phi[0], phi_2[0]])
    exp_history[1].append([phi_1[-1], Phi[-1], phi_2[-1]])

# Plot experimental trajectories
exp_trajectories = get_trajectories(exp_history, GRAIN_INDEXES)
ipf.plot_ipf_trajectory(exp_trajectories, direction, "plot", {"color": "darkgray", "linewidth": 2})
ipf.plot_ipf_trajectory(exp_trajectories, direction, "arrow", {"color": "darkgray", "head_width": 0.01, "head_length": 0.015})
for i, et in enumerate(exp_trajectories):
    ipf.plot_ipf_trajectory([[et[0]]], direction, "scatter", {"color": "darkgray", "s": 8**2})
    ipf.plot_ipf_trajectory([[et[0]]], direction, "text", {"color": "black", "fontsize": 8, "s": start_indexes[GRAIN_INDEXES[i]]+1})
    ipf.plot_ipf_trajectory([[et[-1]]], direction, "text", {"color": "black", "fontsize": 8, "s": end_indexes[GRAIN_INDEXES[i]]+1})

# Read simulated data
sim_dict = csv_to_dict(SIM_PATH)
num_sims = len(sim_dict[list(sim_dict.keys())[0]])
sim_history_list = []
for i in range(num_sims):
    sim_history = [[] for _ in range(2)]
    for j in range(num_grains):
        fields = [f"g{j}_phi_1", f"g{j}_Phi", f"g{j}_phi_2"]
        sim_history[0].append(exp_history[0][j])
        sim_history[1].append([sim_dict[field][i] for field in fields])
    sim_history_list.append(sim_history)

# Plot simulated trajectories
for i, sim_history in enumerate(sim_history_list):
    sim_trajectories = get_trajectories(sim_history, GRAIN_INDEXES)
    ipf.plot_ipf_trajectory([[st[-1]] for st in sim_trajectories], direction, "scatter", {"color": "green", "s": 6**2})

# Save
define_legend(["darkgray", "black"], ["Experimental", "Simulation"], [7, 1.5], ["scatter", "line"])
save_plot("plot_ipf.png")
