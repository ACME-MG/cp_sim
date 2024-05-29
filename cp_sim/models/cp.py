"""
 Title:         Crystal Plasticity Model
 Description:   Contains a Crystal Plasticity model implemented in NEML
 Author:        Janzen Choi

"""

# Libraries
from cp_sim.models.__model__ import __Model__
from neml.cp.crystallography import Lattice
from neml.cp import slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml import elasticity, drivers

# Model class
class Model(__Model__):

    def initialise(self, lattice:Lattice, orientations:list, weights:list, num_threads:int=5,
                   strain_rate:float=1.0e-4, max_strain:float=0.3, youngs:float=190000, poissons:float=0.28) -> None:
        """
        Initialises the model
        
        Parameters:
        * `lattice`:      Lattice object
        * `orientations`: List of initial orientation objects
        * `weights:`      List of weights for the Taylor model
        * `num_threads`:  Number of threads to use to run the model
        * `strain_rate`:  The strain rate
        * `max_strain`:   The maximum strain to run the driver to
        * `youngs`:       The elastic modulus
        * `poissons`:     The poissons ratio
        """
        self.lattice      = lattice
        self.orientations = orientations
        self.weights      = weights
        self.num_threads  = num_threads
        self.strain_rate  = strain_rate
        self.max_strain   = max_strain
        self.e_model      = elasticity.IsotropicLinearElasticModel(youngs, "youngs", poissons, "poissons")
        
    def run_model(self, tau_sat:float, b:float, tau_0:float, gamma_0:float, n:float) -> tuple:
        """
        Runs the model

        Parameters:
        * `tau_sat`: VoceSlipHardening parameter
        * `b`:       VoceSlipHardening parameter
        * `tau_0`:   VoceSlipHardening parameter
        * `gamma_0`: AsaroInelasticity parameter
        * `n`:       AsaroInelasticity parameter

        Returns the single crystal model, polycrystal model, and driver results
        """
        str_model  = slipharden.VoceSlipHardening(tau_sat, b, tau_0)
        slip_model = sliprules.PowerLawSlipRule(str_model, gamma_0, n)
        i_model    = inelasticity.AsaroInelasticity(slip_model)
        k_model    = kinematics.StandardKinematicModel(self.e_model, i_model)
        sc_model   = singlecrystal.SingleCrystalModel(k_model, self.lattice, miter=16, max_divide=2, verbose=False)
        pc_model   = polycrystal.TaylorModel(sc_model, self.orientations, nthreads=self.num_threads, weights=self.weights) # problem
        results    = drivers.uniaxial_test(pc_model, self.strain_rate, emax=self.max_strain, nsteps=500, rtol=1e-6,
                                           atol=1e-10, miter=25, verbose=False, full_results=True)
        return sc_model, pc_model, results
