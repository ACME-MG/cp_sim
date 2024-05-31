"""
 Title:         Simulate
 Description:   Interface for running NEML models
 Author:        Janzen Choi

"""

# Libraries
import itertools, math, numpy as np
from neml.math import rotations
from neml.cp import crystallography
from cp_sim.helper.general import round_sf
from cp_sim.models.__model__ import create_model, __Model__

def get_combinations(params_dict:dict) -> list:
    """
    Returns a list of possible combinations of a set of parameters
    
    Parameters:
    * `params_dict`: Dictionary of parameter lists

    Returns the list of parameter combinations
    """
    param_list = list(params_dict.values())
    combinations = list(itertools.product(*param_list))
    combinations = [list(c) for c in combinations]
    return combinations

def get_model(model_name:str, **kwargs) -> __Model__:
    """
    Defines the model

    Parameters:
    * `model_name`: The name of the model
    """
    model = create_model(model_name, **kwargs)
    return model

def get_lattice(structure:str="fcc") -> crystallography.Lattice:
    """
    Gets the lattice object

    Parameters:
    * `structure`: The crystal structure

    Returns the lattice object
    """
    lattice = crystallography.CubicLattice(1.0)
    if structure == "bcc":
        lattice.add_slip_system([1,1,0], [1,1,1])
    elif structure == "fcc":
        lattice.add_slip_system([1,1,1], [1,1,0])
        lattice.add_slip_system([1,1,1], [1,2,3])
        lattice.add_slip_system([1,1,1], [1,1,2])
    else:
        raise ValueError(f"Crystal structure '{structure}' unsupported!")
    return lattice

def get_orientations(csv_path:str) -> list:
    """
    Given a path to a CSV file, loads the euler-bunge orientations (rads)

    Parameters:
    * `csv_path`: Path to CSV file of grain information

    Returns the list of orientation objects
    """
    grain_stats = np.loadtxt(csv_path, delimiter=",")
    orientations = [rotations.CrystalOrientation(gs[0], gs[1], gs[2], angle_type="radians", convention="bunge") for gs in grain_stats]
    return orientations

def get_weights(csv_path:str) -> list:
    """
    Given a path to a CSV file, loads the weights

    Parameters:
    * `csv_path`: Path to CSV file of grain information

    Returns the list of weight values
    """
    grain_stats = np.loadtxt(csv_path, delimiter=",")
    weights = [gs[3] for gs in grain_stats]
    return weights

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

def get_orientation_history(pc_model:dict, results:dict, inverse:bool=True) -> list:
    """
    Gets the orientation history in euler-bunge form (rads)

    Parameters:
    * `pc_model`: The polycrystal model
    * `results`:  The driver results
    * `inverse`:  Whether to invert the orientations or not
    
    Returns the orientation history
    """
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

def quick_spline(x_list:list, y_list:list, x_value:float) -> float:
    """
    Conducts a quick evaluation using spline interpolation without
    conducting the whole interpolation; assumes that the x_value is
    between min(x_list) and max(x_list) and that x_list is sorted

    Parameters:
    * `x_list`:  The list of x values
    * `y_list`:  The list of y values
    * `x_value`: The x value to evaluate
    
    Returns the evaluated y value
    """
    for i in range(len(x_list)-1):
        if x_list[i] <= x_value and x_value <= x_list[i+1]:
            gradient = (y_list[i+1]-y_list[i])/(x_list[i+1]-x_list[i])
            y_value = gradient*(x_value - x_list[i]) + y_list[i]
            return y_value
    return None

def get_trajectories(euler_history:list, index_list:list=None) -> list:
    """
    Converts a history of euler angles into a list of trajectories

    Parameters:
    * `euler_history`: The history of orientations in euler-bunge form (rads)
    * `index_list`:    The list of indexes to include; if undefined, includes
                       all the trajectories

    Returns the list of trajectories
    """
    trajectories = []
    for i in range(len(euler_history[0])):
        if index_list != None and not i in index_list:
            continue
        trajectory = [euler_state[i] for euler_state in euler_history]
        trajectories.append(trajectory)
    return trajectories

def get_grain_dict(strain_list:list, history:list, grain_ids:list) -> dict:
    """
    Creates a dictionary of grain information

    Parameters:
    * `strain_list`: The list of strain values
    * `history`:     The orientation history
    * `grain_ids`:   The grain indexes to include in the dictionary (starts at 1)
    
    Returns the dictionary of euler-bunge angles (rads)
    """

    # Reformat orientations
    index_list = [grain_id-1 for grain_id in grain_ids] # starts at 0
    trajectories = get_trajectories(history, index_list)

    # Initialise grain dictionary (20%, 40$, .., 100% of max strain)
    strain_intervals = ["0p2", "0p4", "0p6", "0p8", "1p0"]
    grain_dict = {"grain_id": grain_ids}
    for strain_interval in strain_intervals:
        for euler_value in ["phi_1", "Phi", "phi_2"]:
            grain_dict[f"{strain_interval}_{euler_value}"] = []
    
    # Initialise orientation interpolation
    domain = lambda x : round_sf(x, 5) if x>0 else round_sf(x+2*math.pi, 5)
    num_intervals = len(strain_intervals)
    max_strain    = max(strain_list)
    strain_values = [max_strain*(i+1)/num_intervals for i in range(num_intervals)]

    # Get orientations at different strain intervals
    for trajectory in trajectories:

        # Get orientations
        phi_1_list = [domain(t[0]) for t in trajectory]
        Phi_list   = [domain(t[1]) for t in trajectory]
        phi_2_list = [domain(t[2]) for t in trajectory]
        
        # Calculate reduced orientations
        reduced_phi_1_list = [quick_spline(strain_list, phi_1_list, strain) for strain in strain_values]
        reduced_Phi_list   = [quick_spline(strain_list, Phi_list, strain)   for strain in strain_values]
        reduced_phi_2_list = [quick_spline(strain_list, phi_2_list, strain) for strain in strain_values]
        
        # Check that the orientations are valid
        if None in reduced_phi_1_list + reduced_Phi_list + reduced_phi_2_list:
            return None

        # Store reduced orientations
        for i, strain_interval in enumerate(strain_intervals):
            grain_dict[f"{strain_interval}_phi_1"].append(reduced_phi_1_list[i])
            grain_dict[f"{strain_interval}_Phi"].append(reduced_Phi_list[i])
            grain_dict[f"{strain_interval}_phi_2"].append(reduced_phi_2_list[i])

    # Return the dictionary
    return grain_dict
