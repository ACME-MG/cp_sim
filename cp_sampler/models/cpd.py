"""
 Title:         Crystal Plasticity Model
 Description:   Contains a Crystal Plasticity model implemented in NEML
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
from neml.math import rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml import elasticity, drivers

# Model class
class Model:

    def __init__(self, grains_path:str, structure:str="fcc", lattice_a:int=1.0, num_threads:int=5,
                 strain_rate:float=1.0e-4, max_strain:float=0.3, youngs:float=190000, poissons:float=0.28):
        """
        Constructor for the Model class

        Parameters:
        * `grains_path`: Path to the grains file (in euler-bunge notation)
        * `structure`:   Crystal structure ("bcc" or "fcc")
        * `lattice_a`:   The lattice parameter (a)
        * `num_threads`: Number of threads to use to run the model
        * `strain_rate`: The strain rate
        * `max_strain`:  The maximum strain to run the driver to
        * `youngs`:      The elastic modulus
        * `poissons`:    The poissons ratio
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

        # Initialise other parameters
        self.num_threads  = num_threads
        self.strain_rate  = strain_rate
        self.max_strain   = max_strain
        self.youngs       = youngs
        self.poissons     = poissons
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
        e_model = elasticity.IsotropicLinearElasticModel(self.youngs, "youngs", self.poissons, "poissons")
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

    def run_cp_raw(self) -> None:
        """
        Calibrates and runs the crystal plasticity and damage models
        """
        e_model    = self.get_elastic_model()
        str_model  = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau_0)
        slip_model = sliprules.PowerLawSlipRule(str_model, self.gamma_0, self.n)
        i_model    = inelasticity.AsaroInelasticity(slip_model)
        dm_model   = crystaldamage.WorkPlaneDamage()
        dm_func    = crystaldamage.SigmoidTransformation(self.cd, self.beta)
        dm_planar  = crystaldamage.PlanarDamageModel(dm_model, dm_func, dm_func, self.lattice)
        dm_k_model = kinematics.DamagedStandardKinematicModel(e_model, i_model, dm_planar)
        sc_model   = singlecrystal.SingleCrystalModel(dm_k_model, self.lattice, miter=16, max_divide=2, verbose=False)
        pc_model   = polycrystal.TaylorModel(sc_model, self.orientations, nthreads=self.num_threads, weights=self.weights)
        results    = drivers.uniaxial_test(pc_model, self.strain_rate, emax=self.max_strain, nsteps=200, rtol=1e-6, atol=1e-10, miter=25, verbose=False, full_results=True)
        self.model_output = (sc_model, pc_model, results)

    def run_cp(self, try_run:bool=True) -> None:
        """
        Calibrates and runs the crystal plasticity and damage models
        
        Parameters:
        * `try_run`: Wraps a try and except around the code
        """
        if try_run:
            try:
                self.run_cp_raw()
            except:
                self.model_output = None
        else:
            self.run_cp_raw()

    def get_results(self) -> tuple:
        """
        Returns the single crystal model, polycrystal model, and driver
        results from the last model run
        """
        return self.model_output
    
    def get_results_direct(self, param_dict:dict, try_run:bool=True) -> None:
        """
        Runs the model with parameters and returns the results

        Parameters:
        * `param_dict`: Dictionary of parameters
        * `try_run`:    Wraps a try and except around the code

        Returns the single crystal model, polycrystal model, and driver results
        """
        self.define_params(**param_dict)
        self.run_cp(try_run)
        return self.get_results()

    def get_orientation_history(self, inverse:bool=True) -> list:
        """
        Gets the orientation history in euler-bunge form (rads)

        Parameters:
        * `inverse`: Whether to invert the orientations or not

        Returns the orientation history
        """
        pc_model, results = self.model_output[1], self.model_output[2]
        history = np.array(results["history"])
        orientation_history = []
        for state in history:
            orientation_list = []
            for orientation in pc_model.orientations(state):
                euler = list(orientation.to_euler(angle_type="radians", convention="bunge"))
                if inverse:
                    euler = reorient(euler)
                orientation_list.append(euler)
            orientation_history.append(orientation_list)
        return orientation_history

def reorient(euler:list) -> list:
    """
    Inverts the euler angle

    Parameters:
    * `euler`: The euler angle

    Returns the inverted euler angle
    """
    orientation = rotations.CrystalOrientation(euler[0], euler[1], euler[2], angle_type="radians", convention="bunge")
    inverse = orientation.inverse()
    new_euler = inverse.to_euler(angle_type="radians", convention="bunge")
    return new_euler
