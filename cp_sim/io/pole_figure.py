"""
 Title:         Pole figure
 Description:   Contains code to plot orientations on PF and IPF plots;
 References:    https://github.com/Argonne-National-Laboratory/neml/blob/dev/neml/cp/polefigures.py 
 Author:        Janzen Choi

"""

# Libraries
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from neml.math import rotations, tensors
from neml.cp import crystallography
from cp_sim.helper.general import flatten

# Pole figure class
class PF:
    
    def __init__(self, lattice, sample_symmetry:list=[2,2,2], x_direction=[1.0,0,0], y_direction=[0,1.0,0]):
        """
        Constructor for PF class
        
        Parameters:
        * `lattice`:         The lattice object
        * `sample_symmetry`: Sample direction of the projection
        * `x_direction`:     X crystallographic direction of the projection
        * `y_direction`:     Y crystallographic direction of the projection
        """
        self.lattice = lattice
        sample_symmetry_str = "".join([str(ss) for ss in sample_symmetry])
        self.sample_symmetry = crystallography.symmetry_rotations(sample_symmetry_str)
        self.x_direction = [float(x) for x in x_direction]
        self.y_direction = [float(y) for y in y_direction]

    def get_equivalent_poles(self, plane:list) -> list:
        """
        Gets the equivalent poles

        Parameters:
        * `plane`: Plane of the projection

        Returns the list of equivalent poles
        """
        poles = self.lattice.miller2cart_direction(plane)
        eq_poles = self.lattice.equivalent_vectors(poles)
        eq_poles = [pole.normalize() for pole in eq_poles]
        eq_poles = [op.apply(p) for p in eq_poles for op in self.sample_symmetry]
        return eq_poles

    def initialise_polar_grid(self) -> plt.Axes:
        """
        Initialises the polar grid;
        returns the axis
        """
        _, axis = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "polar"})
        axis.grid(False)
        axis.get_yaxis().set_visible(False)
        plt.ylim([0, 1.0])
        plt.xticks([0, np.pi/2], ["  RD", "TD"], fontsize=15)
        return axis

    def get_polar_points(self, euler:list, eq_poles:list) -> list:
        """
        Converts the euler-bunge angles into polar points
        
        Parameters:
        * `euler`:    Orientation in euler-bunge angles (rads)
        * `eq_poles`: List of equivalent poles

        Returns the list of polar points
        """
        standard_rotation = rotations.Orientation(tensors.Vector(self.x_direction), tensors.Vector(self.y_direction))
        orientation = rotations.CrystalOrientation(euler[0], euler[1], euler[2], angle_type="radians", convention="bunge")
        points = [standard_rotation.apply(orientation.inverse().apply(pp)) for pp in eq_poles]
        points = [point for point in points if point[2] >= 0.0] # rid of points in the lower hemisphere
        cart_points = np.array([project_stereographic(v) for v in points])
        polar_points = np.array([cart2pol(cp) for cp in cart_points])
        return polar_points

    def plot_pf(self, euler_list:list, plane:list, colour_list:list=None, size_list:list=None) -> None:
        """
        Plots a standard pole figure using a stereographic projection;
        only works for cubic structures

        Parameters:
        * `euler_list`:  The list of orientations in euler-bunge form (rads)
        * `plane`:       Plane of the projection
        * `colour_list`: List of values to define the colouring scheme
        * `size_list`:   List of values to define the sizing scheme
        
        Plots the PF of the orientations
        """
        
        # Define the colouring and sizing scheme
        rgb_colours = get_colours(euler_list, colour_list)
        norm_size_list = get_sizes(euler_list, size_list)

        # Get the standard rotationn and equivalent poles, and creates the grid
        eq_poles = self.get_equivalent_poles(plane)
        axis = self.initialise_polar_grid()

        # Iterate through the orientations
        for i, euler in enumerate(euler_list):
            polar_points = self.get_polar_points(euler, eq_poles)
            size = norm_size_list[i] if size_list != None else 3
            colour = rgb_colours[i] if colour_list != None else None
            plot_points(axis, polar_points, size, colour)

    def plot_pf_density(self, euler_list:list, plane:list) -> None:
        """
        Plots a standard pole figure with contoured colours based on point density
        
        Parameters:
        * `euler_list`:  The list of orientations in euler-bunge form (rads)
        * `plane`:       Plane of the projection
        """
        eq_poles = self.get_equivalent_poles(plane)
        self.initialise_polar_grid()
        radius_list, theta_list = [], [] 
        for euler in euler_list:
            polar_points = self.get_polar_points(euler, eq_poles)
            radius_list += list(polar_points[:,0])
            theta_list += list(polar_points[:,1])
        sns.kdeplot(x=radius_list, y=theta_list, cmap="viridis", levels=5, thresh=0)

# Inverse pole figure class
class IPF:

    def __init__(self, lattice, sample_symmetry:list=[2,2,2], x_direction=[1,0,0], y_direction=[0,1,0]):
        """
        Constructor for IPF class

        Parameters:
        * `lattice`:         The lattice object
        * `sample_symmetry`: Sample direction of the projection
        * `x_direction`:     X crystallographic direction of the projection
        * `y_direction`:     Y crystallographic direction of the projection
        """
        self.lattice = lattice
        sample_symmetry_str = "".join([str(ss) for ss in sample_symmetry])
        self.sample_symmetry = crystallography.symmetry_rotations(sample_symmetry_str)
        self.vectors = (np.array([0,0,1.0]), np.array([1.0,0,1]), np.array([1.0,1,1])) # force float
        self.norm_vectors = [vector / np.linalg.norm(np.array(vector)) for vector in self.vectors]
        self.x_direction = [int(x) for x in x_direction]
        self.y_direction = [int(y) for y in y_direction]

    def project_ipf(self, quaternion:np.array, direction:list) -> None:
        """
        Projects a single sample direction onto a crystal

        Parameters:
        * `quaternion`:       Orientation in quaternion form
        * `direction`:        Direction of the projection

        Returns the projected points
        """

        # Normalise lattice directions
        norm_x  = self.lattice.miller2cart_direction(self.x_direction).normalize()
        norm_y  = self.lattice.miller2cart_direction(self.y_direction).normalize()
        if not np.isclose(norm_x.dot(norm_y), 0.0):
            raise ValueError("Lattice directions are not orthogonal!")
        norm_z = norm_x.cross(norm_y)
        trans  = rotations.Orientation(np.vstack((norm_x.data, norm_y.data, norm_z.data)))
        norm_d = tensors.Vector(np.array(direction)).normalize()

        # Populate the points
        points_list = []
        for rotation in self.sample_symmetry:
            sample_points = rotation.apply(norm_d)
            crystal_points = quaternion.apply(sample_points)
            points = [trans.apply(op.apply(crystal_points)).data for op in self.lattice.symmetry.ops]
            points_list += points

        # Format the points in the upper hemisphere and return
        points_list = np.array(points_list)
        points_list = points_list[points_list[:,2] > 0]
        return points_list

    def reduce_points_triangle(self, points:tuple) -> list:
        """
        Reduce points to a standard stereographic triangle

        Parameters:
        * `points`: The projected points
        
        Returns the reduced points
        """
        norm_0 = np.cross(self.norm_vectors[0], self.norm_vectors[1])
        norm_1 = np.cross(self.norm_vectors[1], self.norm_vectors[2])
        norm_2 = np.cross(self.norm_vectors[2], self.norm_vectors[0])
        reduced_points = [p for p in points if np.dot(p, norm_0) >= 0 and np.dot(p, norm_1) >= 0 and np.dot(p, norm_2) >= 0]
        return reduced_points

    def initialise_ipf(self) -> plt.Axes:
        """
        Initialises the format and border of the IPF plot

        Returns the axis
        """
        # Create the plot
        axis = plt.subplot(111)
        axis.axis("off")
        plt.text(0.10, 0.11, "[1 0 0]", transform=plt.gcf().transFigure)
        plt.text(0.86, 0.11, "[1 1 0]", transform=plt.gcf().transFigure)
        plt.text(0.74, 0.88, "[1 1 1]", transform=plt.gcf().transFigure)

        # Create the grid
        for i,j in ((0,1), (1,2), (2,0)):
            fs = np.linspace(0, 1, 100)
            points = np.array([project_stereographic((f*self.vectors[i]+(1-f)*self.vectors[j]) / 
                                                          np.linalg.norm(f*self.vectors[i]+(1-f)*self.vectors[j])) for f in fs])
            plt.plot(points[:,0], points[:,1], color="k", linewidth=1)

        # Returns the axis
        return axis

    def get_points(self, euler:list, direction:list) -> list:
        """
        Converts an euler orientation into stereo points

        Parameters:
        * `euler`:          The orientation in euler-bunge form (rads)
        * `direction`:      Direction of the projection

        Returns a list of stereo points
        """
        orientation = rotations.CrystalOrientation(euler[0], euler[1], euler[2], angle_type="radians", convention="bunge")
        projected_points = self.project_ipf(orientation, direction)
        projected_points = np.vstack(tuple(projected_points))
        reduced_points = self.reduce_points_triangle(projected_points)
        stereo_points  = np.array([project_stereographic(point) for point in reduced_points])
        return stereo_points

    def plot_ipf(self, euler_list:list, direction:list, colour_list:list=None, size_list:list=None) -> None:
        """
        Plot an inverse pole figure given a collection of discrete points;
        only works for cubic structures

        Parameters:
        * `euler_list`:  The list of orientations in euler-bunge form (rads)
        * `direction`:   Direction of the projection
        * `colour_list`: List of values to define the colouring scheme
        * `size_list`:   List of values to define the sizing scheme

        Plots the IPF of the orientations
        """

        # Initialise
        rgb_colours = get_colours(euler_list, colour_list)
        norm_size_list = get_sizes(euler_list, size_list)
        axis = self.initialise_ipf()

        # Iterate and plot the orientations
        for i, euler in enumerate(euler_list):
            points = self.get_points(euler, direction)
            size = norm_size_list[i] if size_list != None else 3
            colour = rgb_colours[i] if colour_list != None else None
            plot_points(axis, points, size, colour)

    def plot_ipf_trajectory(self, trajectories:list, direction:list, function:str="scatter", settings:dict={}) -> None:
        """
        Plot an inverse pole figure to display the reorientation trajectories

        Parameters:
        * `trajectories`: The list of reorientation trajectories
        * `direction`:    Direction of the projection
        * `function`      Function to plot points
        * `settings`:     The plotting settings
        
        Plots the IPF of the reorientation trajectories
        """
        axis = self.initialise_ipf()
        for trajectory in trajectories:
            points = np.array(flatten([self.get_points(euler, direction) for euler in trajectory]))
            if function == "arrow": # experimental
                axis.arrow(points[-3,0], points[-3,1], points[-1,0]-points[-3,0], points[-1,1]-points[-3,1], **settings)
            elif function == "text":
                axis.text(points[0,0], points[0,1], **settings)
            else:
                plotter = getattr(axis, function)
                plotter(points[:,0], points[:,1], **settings)

def get_colours(orientations:list, values:list) -> list:
    """
    Checks the colour list and returns a list of RGB colours

    Parameters:
    * `orientations`: The list of orientations
    * `values`:       The list of value

    Returns the list of colours
    """

    # Checks the values
    if not isinstance(values, list):
        return None
    if len(values) != len(orientations):
        raise ValueError("The 'colour_list' does not have the same number of values as the quaternions!")
    
    # Normalise values
    norm_values = np.array(values)
    norm_values = (norm_values - min(norm_values)) / (max(norm_values) - min(norm_values))
    
    # Define colours and return
    colour_map = plt.get_cmap("coolwarm")
    colours = [colour_map(norm_value) for norm_value in norm_values]
    return np.array(colours)

def get_sizes(orientations:list, values:list) -> list:
    """
    Checks the colour list and returns a list of RGB colours

    Parameters:
    * `orientations`: The list of orientations
    * `values`:       The list of value

    Returns the list of colours
    """

    # Checks values
    if values == None:
        return None
    if len(values) != len(orientations):
        raise ValueError("The 'size_list' does not have the same number of values as the quaternions!")
    
    # Normalise sizes and return
    norm_size_list = normalise(values)
    return norm_size_list

def project_stereographic(vector:np.array) -> np.array:
    """
    Stereographic projection of the given vector into a numpy array

    Parameters:
    * `vector`: Unprojected vector
    
    Returns the projected vector
    """
    return np.array([vector[0]/(1.0+vector[2]), vector[1]/(1.0+vector[2])])

def cart2pol(cart_point:np.array):
    """
    Convert a cartesian point into polar coordinates

    Parameters:
    * `cart_point`: Cartesian point

    Returns the polar coordinates
    """
    return np.array([np.arctan2(cart_point[1], cart_point[0]), np.linalg.norm(cart_point)])

def normalise(value_list:list, min_norm:float=1.0, max_norm:list=32.0) -> list:
    """
    Normalises a list of values

    Parameters:
    * `value_list`: The list of values
    * `min_norm`:   The minimum value for the normalised list of values
    * `max_norm`:   The maximum value for the normalised list of values

    Returns the normalised list
    """
    min_value = min(value_list)
    max_value = max(value_list)
    normalised = [min_norm+((value-min_value)/(max_value-min_value))*(max_norm-min_norm) for value in value_list]
    return normalised

def plot_contour(axis:plt.Axes, points:list) -> None:
    """
    Plots the points on a point with contours

    Parameters:
    * `axis`:    The axis to plot the points on
    * `points`:  The points to be plotted
    """
    axis.contourf(points[:,0], points[:,1], points[:,0])

def plot_points(axis:plt.Axes, points:list, size:float, colour:np.ndarray, contour:bool=False) -> None:
    """
    Plots the points on a plot

    Parameters:
    * `axis`:    The axis to plot the points on
    * `points`:  The points to be plotted
    * `size`:    The size of the points
    * `colour`:  The colour of the points; None if not defined
    """
    if not isinstance(colour, np.ndarray):
        axis.scatter(points[:,0], points[:,1], c="black", s=size**2)
    else:
        axis.scatter(points[:,0], points[:,1], color=colour, edgecolor="black", linewidth=0.25, s=size**2)
