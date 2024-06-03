"""
 Title:         617 Runner
 Description:   Runs simulation for 617
 Author:        Janzen Choi

"""

# Libraries
import sys; sys.path += [".."]
import numpy as np
import cp_sim.simulate as sim 
from cp_sim.helper.general import dict_to_csv, round_sf, transpose
from cp_sim.io.pole_figure import PF
from cp_sim.io.plotter import save_plot
from cp_sim.simulate import get_lattice

# Paths
GRAINS_PATH  = f"data/old_grains.csv"
RESULTS_PATH = "results"

# Constants
NUM_THREADS = 12
STRAIN_RATE = 1e-4
MAX_STRAIN  = 1.0
MAX_TIME    = 3000 # seconds
THIN_AMOUNT = 100

# Get grain IDs
grain_ids = [i+1 for i in range(len(sim.get_orientations(GRAINS_PATH)))]

# Get model
model = sim.get_model(
    model_name   = "cp",
    lattice      = sim.get_lattice("fcc"),
    orientations = sim.get_orientations(GRAINS_PATH, "degrees"),
    weights      = sim.get_weights(GRAINS_PATH),
    num_threads  = NUM_THREADS,
    strain_rate  = STRAIN_RATE,
    max_strain   = MAX_STRAIN,
    youngs       = 211000,
    poissons     = 0.3,
)

# Initialise
param_dict = {"tau_sat": 95, "b": 0.25, "tau_0": 820, "gamma_0": round_sf(STRAIN_RATE/3, 4), "n": 4.5}
results_path = f"{RESULTS_PATH}/old"

# Runs the model
status = model.run(param_dict, max_time=MAX_TIME)
if status != "success":
    dict_to_csv(param_dict, f"{results_path}_{status}.csv")
    exit()
_, pc_model, results = model.get_output()

# Save final orientations (1)
strain_list = [round_sf(s[0], 5) for s in results["strain"]]
history     = sim.get_orientation_history(pc_model, results, False)
grain_dict  = sim.get_grain_dict(strain_list, history, grain_ids)
orientation_list_1 = [grain_dict[phi] for phi in ["1p0_phi_1", "1p0_Phi", "1p0_phi_2"]]
orientation_list_1 = transpose(orientation_list_1)

# Save final orientations (2)
final_state = np.array(results["history"])[-1]
orientation_list_2 = []
for orientation in pc_model.orientations(final_state):
    euler = list(orientation.to_euler(angle_type="radians", convention="bunge"))
    # euler = sim.reorient(euler)
    orientation_list_2.append(euler)

# Plot experimental final orientations
direction = [[1,0,0], [1,1,0], [1,1,1]][0]
pf = PF(get_lattice("fcc"))
pf.plot_pf(orientation_list_1, direction)
save_plot("results/pf_1.png")
pf.plot_pf(orientation_list_2, direction)
save_plot("results/pf_2.png")

