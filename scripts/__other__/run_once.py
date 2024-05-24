"""
 Title:         CP Sampler
 Description:   Sampler for the Crystal Plasticity model
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
from cp_sim.models.cp import Model
from cp_sim.helper import round_sf

# # Initialise model
# model = Model(
#     grains_path = "data/p91_s1_grains.csv",
#     structure   = "bcc",
#     lattice_a   = 1.0,
#     num_threads = 12,
#     strain_rate = 1e-4,
#     max_strain  = 0.05,
#     youngs      = 190000,
#     poissons    = 0.28,
# )

# param_names = ["tau_sat", "b", "tau_0", "gamma_0", "n"]
# param_str = """
# 192.35	1.9351	125.91	3.33E-05	2.7478
# """
# param_list = [float(p) for p in param_str.split("\t")]
# param_dict = dict(zip(param_names, param_list))
# _, _, sim_results = model.get_results_direct(param_dict, try_run=False)

# strain_list = [round_sf(s[0], 5) for s in sim_results["strain"]]
# stress_list = [round_sf(s[0], 5) for s in sim_results["stress"]]
# data_dict = {"strain": strain_list, "stress": stress_list}

import matplotlib.pyplot as plt
plt.figure(figsize=(5,5))
plt.gca().set_position([0.17, 0.12, 0.75, 0.75])
plt.plot([0],[0])
plt.scatter([], [], color="grey",  label="Experimental"),
plt.plot([], [],    color="green", label="Surrogate"),
plt.plot([], [],    color="red",   label="CPFE Model"),
legend = plt.legend(framealpha=1, edgecolor="black", fancybox=True, facecolor="white", fontsize=12, loc="upper left")
plt.gca().add_artist(legend)
plt.xlim(0,0.05)
plt.ylim(0,800)
plt.savefig("temp.png")
