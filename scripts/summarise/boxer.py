# Libraries
import numpy as np
import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sampler.helper import csv_to_dict
from cp_sampler.plotter import save_plot

# Constants
SIM_FILE = "summary/phi_bcc.csv"
EXP_FILE = "../data/tensile_p91.csv"

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

# Read experimental and simulation data
sim_dict = csv_to_dict(SIM_FILE)
exp_dict = csv_to_dict(EXP_FILE)

# Define fields to investigate
# grain_indexes = [0,1,2,3,4]
grain_indexes = [5,6,7,8,9]
ori = lambda i : [f"g{i}_phi_1", f"g{i}_Phi", f"g{i}_phi_2"]
field_list = [item for sublist in [ori(i) for i in grain_indexes] for item in sublist]
sim_grid = [sim_dict[field] for field in field_list]
exp_grid = [[exp_dict[field][i] for i in [-1]] for field in field_list]

# Plot boxplots
plot_boxplots(sim_grid, exp_grid, field_list)
save_plot("plot.png")
