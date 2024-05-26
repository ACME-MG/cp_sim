# Libraries
import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sim.helper import csv_to_dict
from cp_sim.plotter import save_plot

# Constants
SUM_ID   = "617_s1"
SIM_FILE = f"results/{SUM_ID}_phi.csv"
EXP_FILE = f"../data/{SUM_ID}_exp.csv"

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
grain_indexes = [75, 189, 314, 338, 463]
phi_list = ["phi_1", "Phi", "phi_2"]
strain_intervals = ["0p2", "0p4", "0p6", "0p8", "1p0"]

# Get experimental data
exp_dict = csv_to_dict(EXP_FILE)
exp_grid = []
for grain_index in grain_indexes:
    exp_list = [exp_dict[f"g{grain_index}_{phi}"] for grain_index in grain_indexes for phi in phi_list]
    exp_grid.append(exp_list)

# Get simulation data
sim_dict = csv_to_dict(SIM_FILE)
sim_grid = []
for grain_index in grain_indexes:
    sim_list = [sim_dict[f"g{grain_index}_{strain_interval}_{phi}"] for grain_index in grain_indexes
                for strain_interval in strain_intervals for phi in phi_list]
    sim_grid.append(sim_list)

# Plot boxplots
field_list = [f"g{grain_index}_{phi}" for grain_index in grain_indexes for phi in phi_list]
plot_boxplots(sim_grid, exp_grid, field_list)
save_plot("results/boxes.png")
