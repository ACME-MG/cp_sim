
# Libraries
import numpy as np, itertools
import threading
from neml.math import rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml import elasticity, drivers

# Constants
MAX_TIME    = 600 # seconds
NUM_THREADS = 12
GRAINS_PATH = "grain_data.csv"
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
            sc_model   = singlecrystal.SingleCrystalModel(dm_k_model, self.lattice)
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

def get_grain_dict(pc_model:dict, history:dict, indexes:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `pc_model`: The polycrystal model
    * `history`:  The history of the model simulation
    * `indexes`:  The grain indexes to include in the dictionary
    
    Returns the dictionary of euler-bunge angles (rads)
    """
    
    # Iterate through each grain
    grain_dict = {"phi_1": [], "Phi": [], "phi_2": []}
    for i in indexes:
        euler_list = [[], [], []]
        
        # Get the trajectory of each grain throughout history
        for state in history:
            orientations = pc_model.orientations(state)
            euler = list(orientations[i].to_euler(angle_type="radians", convention="bunge"))
            for j in range(len(euler_list)):
                euler_list[j].append(euler[j])

        # Store the trajectories
        for j in range(len(euler_list)):
            grain_dict[["phi_1", "Phi", "phi_2"][j]] = euler_list[j]
    
    # Return dictionary
    return grain_dict

def get_stress_dict(pc_model:dict, history:dict, indexes:list) -> dict:
    """
    Creates a dictionary of stress information

    Parameters:
    * `pc_model`: The polycrystal model
    * `history`:  The history of the model simulation
    * `indexes`:  The grain indexes to include in the dictionary
    
    Returns the dictionary of stress values
    """

    # Get stress information
    stress_history = history[:,pc_model.n*sc_model.nstore:pc_model.n*sc_model.nstore+6*pc_model.n:6]
    stress_grid    = [list(stress_list) for stress_list in stress_history]
    
    # Get stress trajectories for each grain
    stress_dict = {"grain_stress": []}
    for i in indexes:
        stress_list = [stress_grid[j][i] for j in range(len(history))]
        index_list = list(range(len(stress_list)))
        polynomial = list(np.polyfit(index_list, stress_list, 6))
        polynomial_str = " ".join([str(round_sf(coef, 5)) for coef in polynomial])
        stress_dict["grain_stress"].append(polynomial_str)
    
    # Returns the dictionary
    return stress_dict

def get_top(value_list:list, num_values:int) -> tuple:
    """
    Gets the top values and indexes of a list of values

    Parameters:
    * `value_list`: The list of values
    * `num_values`: The number of values to return

    Returns the list of top values and indexes
    """
    top_value_list = []
    top_index_list = []
    for i in range(len(value_list)):
        value = value_list[i]
        if len(top_value_list) == 0:
            top_value_list.append(value)
            top_index_list.append(i)
            continue
        if value < top_value_list[-1] and len(top_value_list) == num_values:
            continue
        for j in range(len(top_value_list)):
            if value > top_value_list[j]:
                top_value_list.insert(j, value)
                top_index_list.insert(j, i)
                break
        if len(top_value_list) >= num_values:
            top_value_list.pop(-1)
            top_index_list.pop(-1)
    return top_value_list, top_index_list

# Define the model
model = Model(GRAINS_PATH, 1.0, [1,1,1], [1,1,0])

# Define parameter domains
all_params_dict = {
    "tau_sat": [327.23],
    "b":       [0.9201],
    "tau_0":   [78.096],
    "gamma_0": [round_sf(STRAIN_RATE/3, 4)],
    "n":       [4.9745],
    "cd":      [543.36],
    "beta":    [18.909],
}

# Get combinations of domains
param_list = list(all_params_dict.values())
combinations = list(itertools.product(*param_list))
combinations = [list(c) for c in combinations]

# Get information about 10 biggest grains
top_weights, top_indexes = get_top(model.get_weights(), 10)

# Iterate through the parameters
param_names = list(all_params_dict.keys())
for i in range(len(combinations)):

    # Initialise
    index_str = str(i+1).zfill(3)
    param_dict = dict(zip(param_names, combinations[i]))
    results_path = f"results/once"

    # Prepare the thread for the function
    model.define_params(**param_dict)
    thread = threading.Thread(target=model.run_cp)
    thread.start()
    thread.join(timeout=MAX_TIME)

    # Runs the model for the parameter set
    if thread.is_alive():
        dict_to_csv(param_dict, f"{results_path}_timeout.csv")
        continue
    model_output = model.get_results()
    if model_output == None:
        dict_to_csv(param_dict, f"{results_path}_failed.csv")
        continue
    sc_model, pc_model, results = model_output

    # Get tensile curve
    strain_list = [round_sf(s[0], 5) for s in results["strain"]]
    stress_list = [round_sf(s[0], 5) for s in results["stress"]]
    data_dict = {
        "strain": strain_list,
        "stress": stress_list,
    }

    # Get grain and stress information
    history = np.array(results["history"])
    grain_dict = get_grain_dict(pc_model, history, top_indexes)
    # stress_dict = get_stress_dict(pc_model, history, top_indexes)

    # Compile results and write to CSV file
    combined_dict = {**param_dict, **data_dict, **grain_dict}
    dict_to_csv(combined_dict, f"{results_path}.csv")
