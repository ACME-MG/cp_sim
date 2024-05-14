import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sampler.helper import csv_to_dict
from cp_sampler.plotter import save_plot

# SUMMARY_FILE = "../summary/phi_bcc.csv"
SUMMARY_FILE = "bbbbb.csv"

data_dict = csv_to_dict(SUMMARY_FILE)

grain_id = 0
input_list = ["gamma_0"]
output_list = [f"g{grain_id}_phi_1", f"g{grain_id}_Phi", f"g{grain_id}_phi_2"]

for input in input_list:
    for output in output_list:
        plt.figure(figsize=(5,5))
        plt.scatter(data_dict[input], data_dict[output])
        plt.xlabel(input)
        plt.ylabel(output)
        plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
        plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
        save_plot(f"{input}___{output}")
