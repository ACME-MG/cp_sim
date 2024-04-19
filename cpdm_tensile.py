
# Libraries
import numpy as np, itertools, sys
from neml.math import rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml import elasticity, drivers

# Constants
NUM_THREADS = 8
GRAINS_PATH = "orientations.csv"
STRAIN_RATE = 1.0e-4
MAX_STRAIN  = 1.0
YOUNGS      = 211000
POISSONS    = 0.30
LATTICE     = 1.0

# Model class
class Model:

    def __init__(self, grains_path:str, lattice_a:int=1, slip_plane:list=[1,1,1],
                 slip_direction:list=[1,1,0]):
        """
        Constructor for the RunModel class

        Parameters:
        * `grains_path`:    Path to the grains file (in euler-bunge notation)
        * `lattice_a`:      The lattice type (slip=0, twin=1)
        * `slip_plane`:     The plane of the slip system
        * `slip_direction`: The direction of the slip system
        """

        # Create grain information
        grain_stats = np.loadtxt(grains_path, delimiter=",")
        self.orientations = [rotations.CrystalOrientation(gs[0], gs[1], gs[2], angle_type="degrees", convention="bunge") for gs in grain_stats]
        self.weights = [gs[3] for gs in grain_stats]
        
        # Create lattice
        self.lattice = crystallography.CubicLattice(lattice_a)
        self.lattice.add_slip_system(slip_direction, slip_plane)

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

    def run_cp_dm(self, tau_sat:float, b:float, tau_0:float, gamma_0:float, n:float, cd:float, beta:float) -> tuple:
        """
        Calibrates and runs the crystal plasticity and damage models

        Parameters:
        * `b`:       VoceSlipHardening parameter
        * `tau_sat`: VoceSlipHardening parameter
        * `tau_0`:   VoceSlipHardening parameter
        * `gamma_0`: AsaroInelasticity parameter
        * `n`:       AsaroInelasticity parameter
        * `cd`:      Critical damage value
        * `beta`:    Abruptness of damage onset
        
        Returns the single crystal damage model, crystal plasticity damage model,
        and the dictionary output of the NEML driver as a tuple
        """
        e_model     = self.get_elastic_model()
        strength_model = slipharden.VoceSlipHardening(tau_sat, b, tau_0)
        slip_model  = sliprules.PowerLawSlipRule(strength_model, gamma_0, n)
        i_model     = inelasticity.AsaroInelasticity(slip_model)
        dm_model    = crystaldamage.WorkPlaneDamage()
        dm_func     = crystaldamage.SigmoidTransformation(cd, beta)
        dm_planar   = crystaldamage.PlanarDamageModel(dm_model, dm_func, dm_func, self.lattice)
        dm_k_model  = kinematics.DamagedStandardKinematicModel(e_model, i_model, dm_planar)
        sc_dm_model = singlecrystal.SingleCrystalModel(dm_k_model, self.lattice)
        cp_dm_model = polycrystal.TaylorModel(sc_dm_model, self.orientations, nthreads=NUM_THREADS, weights=self.weights)
        results     = drivers.uniaxial_test(cp_dm_model, STRAIN_RATE, emax=MAX_STRAIN, nsteps=200, rtol=1e-6, atol=1e-10, miter=25, verbose=False, full_results=True)
        return sc_dm_model, cp_dm_model, results

def round_sf(value:float, sf:int) -> float:
    """
    Rounds a float to a number of significant figures

    Parameters:
    * `value`: The value to be rounded
    * `sf`:    The number of significant figures

    Returns the rounded number
    """
    format_str = "{:." + str(sf) + "g}"
    rounded_value = float(format_str.format(value))
    return rounded_value

def dict_to_csv(data_dict:dict, csv_path:str) -> None:
    """
    Converts a dictionary to a CSV file
    
    Parameters:
    * `data_dict`: The dictionary to be converted
    * `csv_path`: The path that the CSV file will be written to
    """
    
    # Extract headers and turn all values into lists
    headers = data_dict.keys()
    for header in headers:
        if not isinstance(data_dict[header], list):
            data_dict[header] = [data_dict[header]]
    
    # Open CSV file and write headers
    csv_fh = open(csv_path, "w+")
    csv_fh.write(",".join(headers) + "\n")
    
    # Write data and close
    max_list_size = max([len(data_dict[header]) for header in headers])
    for i in range(max_list_size):
        row_list = [str(data_dict[header][i]) if i < len(data_dict[header]) else "" for header in headers]
        row_str = ",".join(row_list)
        csv_fh.write(row_str + "\n")
    csv_fh.close()

# Define the model
model = Model(GRAINS_PATH, 1.0, [1,1,1], [1,1,0])

# Grab command line arguments for parallelism
index_1 = int(sys.argv[1])
index_2 = int(sys.argv[2])

# Define parameter domains
all_params_dict = {
    "tau_sat": [1, 500, 1000, 1500, 2000],
    "b":       [[0.1, 1, 10, 100][index_1]],
    "tau_0":   [50, 150, 250, 350, 450],
    "gamma_0": [round_sf(STRAIN_RATE/3, 4)],
    "n":       [1, 2, 3, 4, 5],
    "cd":      [200, 400, 600, 800, 1000],
    "beta":    [[5, 10, 15, 20, 25, 30][index_2]],
}

# Get combinations of domains
param_list = list(all_params_dict.values())
combinations = list(itertools.product(*param_list))
combinations = [list(c) for c in combinations]

# Iterate through the parameters
param_names = list(all_params_dict.keys())
for i in range(len(combinations)):
    param_dict = dict(zip(param_names, combinations[i]))

    # Runs the model for the parameter set
    try:
        sc_model, cp_model, results = model.run_cp_dm(**param_dict)
    except:
        data_dict = {"params": combinations[i]}
        dict_to_csv(data_dict, f"results/{index_1}_{index_2}_{i}")
        continue

    # Compile results and write to CSV file
    strain_list = [s[0] for s in results["strain"]]
    stress_list = [s[0] for s in results["stress"]]
    data_dict = {"params": combinations[i], "strain": strain_list, "stress": stress_list}
    dict_to_csv(data_dict, f"results/{index_1}_{index_2}_{i}.csv")
