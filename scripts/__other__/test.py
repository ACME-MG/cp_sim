"""
 Title:         617 Runner
 Description:   Runs simulation for 617
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
import cp_sim.simulate as sim
from cp_sim.helper.general import csv_to_dict, dict_to_csv, round_sf, get_thinned_list, transpose
from cp_sim.io.pole_figure import IPF
from cp_sim.io.plotter import save_plot

# Paths
EXP_PATH     = f"data/617_s1_exp.csv"
GRAINS_PATH  = f"data/617_s1_grains.csv"
RESULTS_PATH = "results"

# Constants
NUM_THREADS = 12
STRAIN_RATE = 1e-4
MAX_STRAIN  = 1.0
MAX_TIME    = 3000 # seconds
THIN_AMOUNT = 100

# Get microstructure information
exp_dict = csv_to_dict(EXP_PATH) # passive
orientations = sim.get_orientations(GRAINS_PATH, angle_type="radians", is_passive=True) # passive to active
weights = sim.get_weights(GRAINS_PATH)
grain_ids = [i+1 for i in range(len(orientations))]

# Get model
model = sim.get_model(
    model_name   = "cp",
    lattice      = sim.get_lattice("fcc"),
    orientations = orientations,
    weights      = weights,
    num_threads  = NUM_THREADS,
    strain_rate  = STRAIN_RATE,
    max_strain   = MAX_STRAIN,
    youngs       = 211000,
    poissons     = 0.3,
)

# Initialise
param_dict = {"tau_sat": 95, "b": 0.25, "tau_0": 820, "gamma_0": round_sf(STRAIN_RATE/3, 4), "n": 4.5}
results_path = f"{RESULTS_PATH}/617"

# Runs the model
status = model.run(param_dict, max_time=MAX_TIME)
if status != "success":
    dict_to_csv(param_dict, f"{results_path}_{status}.csv")
    exit()
sc_model, pc_model, results = model.get_output() # active to passive

# Get stress in individual crystals
n_hist = sc_model.nstore
history = np.array(results['history'])
final_stresses = list(history[:,pc_model.n*n_hist:pc_model.n*n_hist+6*pc_model.n:6][-1])

# Process results
strain_list = [round_sf(s[0], 5) for s in results["strain"]]
stress_list = [round_sf(s[0], 5) for s in results["stress"]]
history     = sim.get_orientation_history(pc_model, results, inverse=False) # passive
grain_dict  = sim.get_grain_dict(strain_list, history, grain_ids)
data_dict   = {"strain": get_thinned_list(strain_list, THIN_AMOUNT), "stress": get_thinned_list(stress_list, THIN_AMOUNT)}

# Plot IPF
direction = [[1,0,0], [1,1,0], [1,1,1]][0]
ipf = IPF(sim.get_lattice("fcc"))
orientation_list = [grain_dict[phi] for phi in ["1p0_phi_1", "1p0_Phi", "1p0_phi_2"]]
orientation_list = transpose(orientation_list)
ipf.plot_ipf(orientation_list, direction, colour_list=final_stresses) # passive
save_plot("617_rotated.png")
