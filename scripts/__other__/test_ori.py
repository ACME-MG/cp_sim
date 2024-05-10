from neml.math import rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml import elasticity, drivers
import numpy as np

# euler_1 = [math.pi/2, math.pi/4, math.pi]
# orientation_1 = [rotations.CrystalOrientation(euler_1[0], euler_1[1], euler_1[2], angle_type="radians", convention="bunge")]
orientation_1 = rotations.random_orientations(1)
lattice = crystallography.CubicLattice(1.0)
lattice.add_slip_system([1,1,0], [1,1,1])

e_model     = elasticity.IsotropicLinearElasticModel(211000, "youngs", 0.3, "poissons")
str_model   = slipharden.VoceSlipHardening(50, 0.1, 50)
slip_model  = sliprules.PowerLawSlipRule(str_model, 3.333e-05, 2)
i_model     = inelasticity.AsaroInelasticity(slip_model)
k_model     = kinematics.StandardKinematicModel(e_model, i_model)
sc_model    = singlecrystal.SingleCrystalModel(k_model, lattice, miter=16, max_divide=2, verbose=False)
pc_model    = polycrystal.TaylorModel(sc_model, orientation_1, nthreads=8)
results     = drivers.uniaxial_test(pc_model, 1e-4, emax=0.10, nsteps=200, rtol=1e-6, atol=1e-10, miter=25, verbose=False, full_results=True)
first_state = np.array(results["history"])[0]
orientation_2 = pc_model.orientations(first_state)[0]

euler_1 = list(orientation_1[0].to_euler(angle_type="radians", convention="bunge"))
euler_2 = list(orientation_2.inverse().to_euler(angle_type="radians", convention="bunge"))
print(euler_1)
print(euler_2)
