# Libraries
import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sim.helper.general import csv_to_dict
from cp_sim.io.plotter import save_plot

# Constants
GRAIN_IDS = [75, 189, 314, 346, 463]
SUM_ID   = "617_s1"
SIM_PATH = f"results/{SUM_ID}_phi.csv"
EXP_PATH = f"../data/{SUM_ID}_exp.csv"

def plot_boxplots(point_grid:list, ideal_grid:list, field_list:list) -> None:
    """
    Plots boxplots horizontally

    Parameters:
    * `point_grid`: The list of list of points to form the boxplots
    * `ideal_grid`: The list of list of points to plot atop the boxplots
    * `field_list`: The list of fields
    """

    # Initialise plots
    fig, axes = plt.subplots(ncols=len(point_grid), nrows=1, figsize=(len(point_grid)*2, 5))
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.8, wspace=0.6)

    # Plot all objective values
    for i, axis in enumerate(axes):
        box_plots = axis.boxplot(point_grid[i], patch_artist=True, showfliers=False, widths=0.8, zorder=1)
        box_plots["boxes"][0].set(linewidth=1, edgecolor="black", facecolor=[j/255 for j in (255, 193, 184)])
        box_plots["medians"][0].set_color("black")
        for value in ideal_grid[i]:
            axis.plot([0.5, 1.5], [value]*2, color="blue", zorder=3) 
        
        # Apply general format
        axis.set_title(field_list[i])
        axis.tick_params(axis="y", labelsize=13)
        axis.yaxis.major.formatter._useMathText = True
        axis.set_xticks([])
        axis.set_xticklabels([])
        axis.grid(axis="y")

# Define fields
phi_list = ["phi_1", "Phi", "phi_2"]
strain_intervals = ["0p2", "0p4", "0p6", "0p8", "1p0"]

# Get experimental data
exp_dict = csv_to_dict(EXP_PATH)
exp_grid = []
for grain_id in GRAIN_IDS:
    for phi in phi_list:
        exp_grid.append(exp_dict[f"g{grain_id}_{phi}"])

# Get simulation data
sim_dict = csv_to_dict(SIM_PATH)
sim_grid = []
flatten = lambda nested_list : [item for sublist in nested_list for item in sublist]
for grain_id in GRAIN_IDS:
    for phi in phi_list:
        sim_list = [sim_dict[f"g{grain_id}_{strain_interval}_{phi}"] for strain_interval in strain_intervals]
        sim_grid.append(flatten(sim_list))

# Plot boxplots
field_list = [f"g{grain_id}_{phi}" for grain_id in GRAIN_IDS for phi in phi_list]
plot_boxplots(sim_grid, exp_grid, field_list)
save_plot("results/boxes.png")
