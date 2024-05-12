import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sampler.helper import csv_to_dict
from cp_sampler.plotter import save_plot

SUMMARY_FILE = "../summary/phi_bcc.csv"

data_dict = csv_to_dict(SUMMARY_FILE)

input_list = ["b"]
output_list = ["g1_phi_1_end", "g1_Phi_end", "g1_phi_2_end",
               "g2_phi_1_end", "g2_Phi_end", "g2_phi_2_end",
               "g3_phi_1_end", "g3_Phi_end", "g3_phi_2_end",
               "g4_phi_1_end", "g4_Phi_end", "g4_phi_2_end",
               "g5_phi_1_end", "g5_Phi_end", "g5_phi_2_end"]

for input in input_list:
    for output in output_list:
        plt.scatter(data_dict[input], data_dict[output])
        save_plot(f"{input}___{output}")
