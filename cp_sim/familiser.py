"""
 Title:         Familiser
 Description:   Contains code for grouping grain orientations into grain families
 References:    http://www.ebsd.info/pdf/TablesOfTextureAnalysis.pdf
 Author:        Janzen Choi

"""

# Libraries
import math
from orientation import deg_to_rad
from csl import get_disorientation

def get_magnitude(vector:list) -> float:
    """
    Calculates the magnitude of a vector

    Parameters:
    * `vector`: The vector

    Returns the magnitude
    """
    square_sum = sum([math.pow(v, 2) for v in vector])
    magnitude = math.sqrt(square_sum)
    return magnitude

def miller_to_euler(hkl:list, uvw:list) -> list:
    """
    Converts a list of miller indices to their euler-bunge form;
    (plane)[directionn] = (hkl)[uvw] = {hkl}<uvw>

    Parameters:
    * `hkl`: The crystallographic plane
    * `uvw`: The crystallographic direction

    Returns the euler-bunge orientation (rads)
    """
    h, k, l = tuple(hkl)
    u, v, w = tuple(uvw)
    domain = lambda x : x + 2*math.pi if x < 0 else x
    phi_2 = math.atan2(h, k)
    Phi   = math.atan2(k, l*math.cos(phi_2))
    phi_1 = math.atan2(l*w, (k*u-h*v)*math.cos(Phi))
    phi_1 = domain(phi_1)
    Phi   = domain(Phi)
    phi_2 = domain(phi_2)
    return [phi_1, Phi, phi_2]

def get_grain_family(orientations:list, plane:list, direction:list, type:str="cubic",
                     threshold:float=10.0) -> list:
    """
    Groups a list of orientations to a family

    Parameters:
    * `orientations`: The list of orientations (neml.math.rotations.Orientation)
    * `plane`:        The plane
    * `direction`:    The direction
    * `type`:         The type of crystal structure
    * `threshold`:    The misorientation threshold for being part of a family (deg)

    Returns the indices of the grain family
    """
    
    # Calculate the miller orientation and convert threshold
    pd_orientation = miller_to_euler(plane, direction)
    rad_threshold = deg_to_rad(threshold)
    
    # Iterate through grains and add to family
    family_indices = []
    for i in range(len(orientations)):
        orientation = orientations[i].to_euler(angle_type="radians", convention="bunge")
        misorientation = get_disorientation(orientation, pd_orientation, type)
        if misorientation < rad_threshold:
            family_indices.append(i)
    return family_indices

# Testing
# hkl_uvw = [
#     [[0,0,1],  [1,0,0]],   # [0,0,90]
#     [[0,2,1],  [1,0,0]],   # [0,26,90]
#     [[0,1,1],  [1,0,0]],   # [0,45,90]
#     [[0,1,1],  [2,1,1]],   # [35,45,90]
#     [[1,2,3],  [6,3,4]],   # [59,37,63]
#     [[1,1,2],  [1,1,1]],   # [90,35,45]
#     [[4,4,11], [11,11,8]], # [90,27,45]
# ]
# from orientation import rad_to_deg
# for pd in hkl_uvw:
#     print(rad_to_deg(miller_to_euler(pd[0], pd[1])))
#     print("================")
