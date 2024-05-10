"""
 Title:         Crystal Platicity Damage Model
 Description:   Contains a Crystal Plasticity model coupled with a damage model implemented in NEML
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from neml.math import rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml import elasticity, drivers

# Constants
NUM_THREADS = 8
STRAIN_RATE = 1.0e-4
MAX_STRAIN  = 0.5
YOUNGS      = 211000
POISSONS    = 0.30

# Model class
class Model:

    def __init__(self, grains_path:str, structure:str="fcc", lattice_a:int=1.0):
        """
        Constructor for the Model class

        Parameters:
        * `grains_path`: Path to the grains file (in euler-bunge notation)
        * `structure`:   Crystal structure ("bcc" or "fcc")
        * `lattice_a`:   The lattice parameter (a)
        """

        # Create grain information
        grain_stats = np.loadtxt(grains_path, delimiter=",")
        self.orientations = [rotations.CrystalOrientation(gs[0], gs[1], gs[2], angle_type="radians", convention="bunge") for gs in grain_stats]
        self.weights = [gs[3] for gs in grain_stats]
        
        # Create lattice
        self.lattice = crystallography.CubicLattice(lattice_a)
        if structure == "fcc":
            self.lattice.add_slip_system([1,1,0], [1,1,1])
        elif structure == "bcc":
            self.lattice.add_slip_system([1,1,1], [1,1,0])
            self.lattice.add_slip_system([1,1,1], [1,2,3])
            self.lattice.add_slip_system([1,1,1], [1,1,2])

        # Initialise results
        self.model_output = None

    def get_lattice(self) -> crystallography.CubicLattice:
        """
        Returns the lattice
        """
        return self.lattice
    
    def get_weights(self) -> list:
        """
        Returns the weights
        """
        return self.weights

    def get_elastic_model(self):
        """
        Returns the elastic model
        """
        e_model = elasticity.IsotropicLinearElasticModel(YOUNGS, "youngs", POISSONS, "poissons")
        return e_model

    def define_params(self, tau_sat:float, b:float, tau_0:float, gamma_0:float, n:float, cd:float, beta:float) -> None:
        """
        Defines the parameters for the model

        Parameters:
        * `tau_sat`: VoceSlipHardening parameter
        * `b`:       VoceSlipHardening parameter
        * `tau_0`:   VoceSlipHardening parameter
        * `gamma_0`: AsaroInelasticity parameter
        * `n`:       AsaroInelasticity parameter
        * `cd`:      Critical damage value
        * `beta`:    Abruptness of damage onset
        """
        self.tau_sat = tau_sat
        self.b = b
        self.tau_0 = tau_0
        self.gamma_0 = gamma_0
        self.n = n
        self.cd = cd
        self.beta = beta

    def run_cp(self) -> None:
        """
        Calibrates and runs the crystal plasticity and damage models;
        returns the single crystal damage model, crystal plasticity damage model,
        and the dictionary output of the NEML driver as a tuple
        """
        
        # Get the results
        try:
            e_model    = self.get_elastic_model()
            str_model  = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau_0)
            slip_model = sliprules.PowerLawSlipRule(str_model, self.gamma_0, self.n)
            i_model    = inelasticity.AsaroInelasticity(slip_model)
            dm_model   = crystaldamage.WorkPlaneDamage()
            dm_func    = crystaldamage.SigmoidTransformation(self.cd, self.beta)
            dm_planar  = crystaldamage.PlanarDamageModel(dm_model, dm_func, dm_func, self.lattice)
            dm_k_model = kinematics.DamagedStandardKinematicModel(e_model, i_model, dm_planar)
            sc_model   = singlecrystal.SingleCrystalModel(dm_k_model, self.lattice, miter=16, max_divide=2, verbose=False)
            pc_model   = polycrystal.TaylorModel(sc_model, self.orientations, nthreads=NUM_THREADS, weights=self.weights)
            results    = drivers.uniaxial_test(pc_model, STRAIN_RATE, emax=MAX_STRAIN, nsteps=200, rtol=1e-6, atol=1e-10, miter=25, verbose=False, full_results=True)
            self.model_output = (sc_model, pc_model, results)
        except:
            self.model_output = None

    def get_results(self) -> tuple:
        """
        Returns the single crystal model, polycrystal model, and driver
        results from the last model run
        """
        return self.model_output
    
    def get_results_direct(self, param_dict:dict) -> None:
        """
        Runs the model with parameters and returns the results

        Parameters:
        * `param_dict`: Dictionary of parameters

        Returns the single crystal model, polycrystal model, and driver results
        """
        self.define_params(**param_dict)
        self.run_cp()
        return self.get_results()

    def get_orientation_history(self) -> list:
        """
        Returns the orientation history in euler-bunge form (rads)
        """
        pc_model, results = self.model_output[1], self.model_output[2]
        history = np.array(results["history"])
        orientation_history = []
        for state in history:
            orientation_list = []
            for orientation in pc_model.orientations(state):
                inverse = orientation.inverse()
                orientation_list.append(list(inverse.to_euler(angle_type="radians", convention="bunge")))
            orientation_history.append(orientation_list)
        return orientation_history
