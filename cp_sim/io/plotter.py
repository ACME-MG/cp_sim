"""
 Title:         Plotter
 Description:   For plotting data
 Author:        Janzen Choi
 
"""

# Libraries
import matplotlib.pyplot as plt
import matplotlib.colors as mcolours

# Constants
EXP_COLOUR   = "darkgray"
CAL_COLOUR   = "green"
VAL_COLOUR   = "red"
ALL_COLOURS  = list(mcolours.TABLEAU_COLORS) + list(mcolours.BASE_COLORS) + list(mcolours.CSS4_COLORS)

# Plotter class
class Plotter:

    def __init__(self, x_label:str="x", y_label:str="y"):
        """
        Class for plotting data

        Parameters:
        * `x_label`: The label for the x axis
        * `y_label`: The label for the y axis
        """
        self.x_label = x_label
        self.y_label = y_label

    def prep_plot(self, title:str="", size:int=12) -> None:
        """
        Prepares the plot
        
        Parameters:
        * `title`: The title of the plot
        * `size`:  The size of the font
        """

        # Set figure size and title
        plt.figure(figsize=(5,5))
        plt.title(title, fontsize=size+3, fontweight="bold", y=1.05)
        plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
        plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")

        # Set x and y labels
        plt.xlabel(f"{self.x_label.replace('_', ' ').capitalize()}", fontsize=size)
        plt.ylabel(f"{self.y_label.replace('_', ' ').capitalize()}", fontsize=size)
    
    def set_limits(self, x_limits:tuple=None, y_limits:tuple=None) -> None:
        """
        Sets the limits of the x and y scales

        Parameters:
        * `x_limits`: The upper and lower bounds of the plot for the x scale
        * `y_limits`: The upper and lower bounds bound of the plot for the y scale
        """
        if x_limits != None:
            plt.xlim(*x_limits)
        if y_limits != None:
            plt.ylim(*y_limits)

    def set_log_scale(self, x_log:bool=False, y_log:bool=False) -> None:
        """
        Changes the scale of the plot
        
        Parameters:
        * `x_log`: Whether to log the x scale
        * `y_log`: Whether to log the y scale
        """
        if x_log:
            plt.xscale("log")
        if y_log:
            plt.yscale("log")
    
    def scat_plot(self, data_dict:dict, colour:str=EXP_COLOUR, size:int=5, priority:int=1) -> None:
        """
        Plots the experimental data using a scatter plot

        Parameters:
        * `data_dict`: The dictionary to store the data
        * `colour`:    The colour to plot the data
        * `size`:      The size of the curve
        * `priority`:  The priority of the curve
        """
        x_list = data_dict[self.x_label]
        if self.x_label == "time":
            x_list = [x/3600 for x in x_list]
        plt.scatter(x_list, data_dict[self.y_label], s=size**2,
                    marker="o", color=colour, linewidth=1, zorder=priority)
        
    def line_plot(self, data_dict:dict, colour=CAL_COLOUR, priority:int=1) -> None:
        """
        Plots the experimental data using a line plot

        Parameters:
        * `data_dict`: The dictionary to store the data
        * `colour`:    The colour to plot the data
        * `priority`:  The priority of the curve
        """
        x_list = data_dict[self.x_label]
        if self.x_label == "time":
            x_list = [x/3600 for x in x_list]
        plt.plot(x_list, data_dict[self.y_label], colour, zorder=priority)

def define_legend(colour_list:list, label_list:list, size_list:list, type_list:list) -> None:
    """
    Defines the plot legend
    
    Parameters:
    * `colour_list`: The colours in the legend
    * `label_list`:  The keys to add to the legend
    * `size_list`:   The size of the icons in the legend
    * `type_list`:   The type of the icons in the legend
    """
    for i in range(len(colour_list)):
        if type_list[i] == "scatter":
            plt.scatter([], [], color=colour_list[i], label=label_list[i], s=size_list[i]**2)
        elif type_list[i] == "line":
            plt.plot([], [], color=colour_list[i], label=label_list[i], linewidth=size_list[i])
    plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white")

def save_plot(plot_name) -> None:
    """
    Saves the plot and clears the figure

    Parameters:
    * `plot_name`: The name of the plot
    """
    plt.savefig(plot_name)
    plt.cla()
    plt.clf()
    plt.close()
