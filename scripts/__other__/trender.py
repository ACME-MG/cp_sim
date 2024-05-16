import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sim.helper import csv_to_dict
from cp_sim.plotter import save_plot

SUMMARY_FILE = "summary/phi_bcc.csv"
GRAIN_INDEXES = [0]

data_dict = csv_to_dict(SUMMARY_FILE)
input_list = ["tau_sat", "n"]
ori = lambda i : [f"g{i}_phi_1", f"g{i}_Phi", f"g{i}_phi_2"]
output_list = [item for sublist in [ori(i) for i in GRAIN_INDEXES] for item in sublist]

for input in input_list:
    for output in output_list:
        plt.figure(figsize=(5,5))
        plt.scatter(data_dict[input], data_dict[output])
        plt.xlabel(input)
        plt.ylabel(output)
        plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
        plt.gca().grid(which="major", axis="both", color="SlateGray", linewidth=1, linestyle=":")
        save_plot(f"trends/{input}_{output}")
