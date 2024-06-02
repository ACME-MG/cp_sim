
# Libraries
import sys; sys.path += [".."]
import cp_sim.simulate as sim 
from cp_sim.helper.general import round_sf, transpose
from cp_sim.io.pole_figure import IPF
from cp_sim.io.plotter import save_plot
from neml.math import rotations

# Get angles
passive_euler_list = [
    [3.7164707703901287, 2.7248590459733864, 3.4311189115751706],
    [5.5520905264366425, 2.5398823155141548, 5.0600236101222285],
    [6.231773610449561, 2.426458798330133, 5.796235205595585],
    [3.62723709852502, 2.7967654114089644, 2.882792469738934],
    [5.310653498642153, 2.1909620483200536, 4.601299411608237],
]
passive_orientations  = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in passive_euler_list]
active_euler_list   = [sim.reorient(euler) for euler in passive_euler_list]
active_orientations = [rotations.CrystalOrientation(*euler, angle_type="radians", convention="bunge") for euler in active_euler_list]

# Get model
model = sim.get_model(
    model_name   = "cp",
    lattice      = sim.get_lattice("fcc"),
    orientations = active_orientations,
    weights      = [1]*len(passive_euler_list),
    num_threads  = 16,
    strain_rate  = 1e-4,
    max_strain   = 0.1,
    youngs       = 211000,
    poissons     = 0.3,
)

# Run model
param_dict = {"tau_sat": 95, "b": 0.25, "tau_0": 820, "gamma_0": round_sf(1e-4/3, 4), "n": 4.5}
model.run(param_dict, max_time=3600)
_, pc_model, results = model.get_output()

# Process results
strain_list = [round_sf(s[0], 5) for s in results["strain"]]
history     = sim.get_orientation_history(pc_model, results, inverse=False) # passive
orientation_list = history[-1]
print(history[0])

# orientation_list = [grain_dict[phi] for phi in ["1p0_phi_1", "1p0_Phi", "1p0_phi_2"]]
# orientation_list = transpose(orientation_list)
# orientation_list = [sim.reorient(orientation) for orientation in orientation_list] # active -> passive

# # Plot
# ipf = IPF(sim.get_lattice("fcc"))
# ipf.plot_ipf(passive_euler_list, [1,0,0])
# save_plot("initial.png")
# ipf.plot_ipf(orientation_list, [1,0,0])
# save_plot("active.png")
