import matplotlib.pyplot as plt
import sys; sys.path += ["../.."]
from cp_sampler.helper import csv_to_dict
from cp_sampler.plotter import save_plot

SUMMARY_FILE = "../summary/phi_bcc.csv"

data_dict = csv_to_dict(SUMMARY_FILE)

input_list = ["n"]
output_list = ["g1_phi_1", "g1_Phi", "g1_phi_2"]

for input in input_list:
    for output in output_list:
        plt.scatter(data_dict[input], data_dict[output])
        save_plot(f"{input}___{output}")
