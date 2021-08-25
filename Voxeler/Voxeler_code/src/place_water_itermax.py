# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated : January 2021
# ---------------------------------------------------------------------------- #


# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
# Allows Numpy array manipulation
import copy as cp
# Allows true copy of elements
import random as rd

# Program external resources

# Parameters
from config import global_parameters as gp

# Contains the global variabless

# Classes
# General library
# Specific modules

# ---------------------------------------------------------------------------- #


# Main function -------------------------------------------------------------- #

def place_water_itermax(d_parameters, o_system, o_structure):
    """
	Iteratively places water molecules by decreasing order of score and updates the grid
	:param d_parameters: Dictionary of the parameters used for the placement of water molecules
	:param o_system: Object containing various variables used for the placement of water molecules
	:param o_structure: The structure to solubilize
	"""

    # STEP 0 : Preparing variables ---------------------- #
    i_water_code = gp.D_ELEMENT_NUMBER["OOW"]  # Retrieves the element code for water
    f_water_radius = o_system.d_scaled_vdw["OOW"]  # Retrieves the VdW radius of a water molecule, in Angstrom
    l_i_radius_range = list(range(-f_water_radius,
                                  f_water_radius + 1))  # Generates the list of radius inside the WdV sphere of water molecules
    # END STEP 0 ---------------------------------------- #

    # STEP 1 : Determining the functions to be used ----- #
    d_geometry_sphere = {  # Dictionary matching the system geometry to the corresponding sphere formula
        "manhattan": determine_sphere_manhattan,
        "minkowski": determine_sphere_minkowski,
        "euclidean": determine_sphere_euclidean
    }
    d_geometry_distance = {  # Dictionary matching the system geometry to the corresponding distance formula
        "manhattan": compute_distance_manhattan,
        "minkowski": compute_distance_minkowski,
        "euclidean": compute_distance_euclidean
    }
    function_sphere = d_geometry_sphere[
        d_parameters["s_scoring_metric"].lower()]  # Loads the sphere function corresponding to the system geometry
    function_distance = d_geometry_distance[
        d_parameters["s_scoring_metric"].lower()]  # Loads the distance function corresponding to the system geometry
    # END STEP 1 ---------------------------------------- #

    # STEP 2 : Separating structure and water ----------- #
    a_water_grid = cp.deepcopy(o_structure.a_grid)  # Copies the grid containing the structure
    a_water_positions = np.array(np.where(  # Retrieves the position of potential water molecules
        o_structure.a_grid["element_symbol"] == i_water_code
    ))
    o_structure.a_grid[  # Deletes the water positions from the main grid
        a_water_positions[0],  # X coordinates for the potential water molecules
        a_water_positions[1],  # Y coordinates for the potential water molecules
        a_water_positions[2]  # Z coordinates for the potential water molecules
    ] = (0, 0, 0)

    a_structure_positions = np.array(np.where(  # Retrieves the coordinates of each structure atoms
        np.bitwise_and(  # Points must fulfil multiple conditions
            o_structure.a_grid["element_symbol"] > 0,  # Non-empty points
            o_structure.a_grid["element_symbol"] != i_water_code  # Points not containing water molecules
        )
    ))
    a_water_grid[  # Deletes the structure atoms from the water grid
        a_structure_positions[0],  # X coordinates of the structure atoms
        a_structure_positions[1],  # Y coordinates of the structure atoms
        a_structure_positions[2],  # Z coordinates of the structure atoms
    ] = (0, 0, 0)
    # END STEP 2 ---------------------------------------- #

    # STEP 3 : Filtering insufficient scores ------------ #
    a_low_scores = np.array(
        np.where(  # Retrieves the coordinates of each potential position with a score inferior to the minimal one
            a_water_grid["score"] <= max(0, d_parameters["f_min_water_score"])
        ))

    # If there is any position with a invalid score
    if np.any(a_low_scores):
        a_water_grid[  # Deletes the invalid positions
            a_low_scores[0],  # X coordinates of the positions with insufficient score
            a_low_scores[1],  # Y coordinates of the positions with insufficient score
            a_low_scores[2],  # Z coordinates of the positions with insufficient score
        ] = (0, 0, 0)
        a_water_positions = np.array(np.where(  # Updates the position of potential water molecules
            a_water_grid["element_symbol"] == i_water_code
        ))

    del a_low_scores  # Deletes the array of invalid positions
    # End if
    # END STEP 3 ---------------------------------------- #

    # STEP 4 : Sorting the positions by their scores ---- #
    l_all_scores, a_element_positions = sort_element_score(  # Sorts the possibles scores in the element grid
        a_element_grid=a_water_grid,  # The grid containing the element possible positions
        i_element_code=i_water_code,  # The code of the element to place
        f_min_score=d_parameters["f_min_water_score"]  # The minimal score for the placement
    )
    # END STEP 4 ---------------------------------------- #

    # STEP 5 : Determining the shape of the VdW sphere -- #
    a_sphere_shape = function_sphere(  # Generates the array of points in the VdW sphere
        l_i_radius_range=l_i_radius_range,  # All the distances from the atom center
        f_element_radius=f_water_radius  # The VdW radius from the atom center
    )
    # END STEP 5 ---------------------------------------- #

    # STEP 6 : Placing water molecules -- #
    i_score_index = 0  # The current index in the list of sorted scores
    i_repeat_count = 0  # If the list needs to be refined and the processes redone
    i_max_repeat = int(10 / o_system.f_grid_spacing)  # Max number of refining steps
    b_compute_distance = True  # If the distance needs to be calculated

    # While the refinement process is not finished
    while True:

        # Tries to retrieve the scores
        try:
            f_max_score = l_all_scores[i_score_index]  # Retrieves the actual maximum score

        # If all positions have been treated
        except IndexError:
            break  # Ends the loop

        a_max_coord = retrieve_element_coordinates(  # Retrieves the coordinates of points with the actual maximum score
            a_element_grid=a_water_grid,  # The grid containing the elements scores
            a_element_positions=a_water_positions,  # The list of coordinates of points able to contain the element
            f_max_score=f_max_score  # The current maximum score
        )

        # If there is at least one position with the current maximum score
        if np.any(a_max_coord):

            # If there is more than one position for this score
            if len(a_max_coord) > 1:

                # If the distance needs to be computed
                if b_compute_distance:

                    l_distances = [  # Computes the distance between each max point List of voxels with max score
                        tuple(function_distance(
                            a_coords_a=l_point,  # One of the max coord point
                            a_coords_b=a_max_coord  # The array of coordinates for points with the current maximum score
                        )) for l_point in a_max_coord  # For each point in the list of coordinates
                    ]
                    l_no_contact = l_distances > 2 * f_water_radius  # Retrieves the list of distances superior to the element VdW radius
                    # TODO
                    l_isolated_indexes = [
                        len(l_distances) - 1 == np.count_nonzero(a_position) for a_position in l_no_contact
                        # Retrieves the indexes of isolated positions
                    ]
                    l_isolated_coord = [
                        a_max_coord[a_position] for a_position in np.where(l_isolated_indexes)[0]
                        # Retrieves the coordinates of isolated positions
                    ]

                    # If there is any isolated coordinates
                    if np.any(l_isolated_coord):

                        # For each isolated coord
                        for a_atom_center in l_isolated_coord:
                            b_placed = place_element_volume(  # Clears the volume around the placed element
                                a_structure_grid=o_structure.a_grid,  # Grid of the main structure
                                a_element_grid=a_water_grid,  # Grid of the available positions for the element
                                a_atom_center=a_atom_center,  # Center of the element to place
                                a_sphere_shape=a_sphere_shape  # Spherical array of points
                            )
                            if b_placed:
                                update_dic(a_atom_center, f_max_score, gp.D_WATER_SCORING)
                                b_placed = False


                    # If there is no isolated coordinates
                    else:
                        a_atom_center = a_max_coord[
                            np.random.randint(len(a_max_coord))]  # Selects a random position for this element
                        b_placed = place_element_volume(  # Clears the volume around the placed element
                            a_structure_grid=o_structure.a_grid,  # Grid of the main structure
                            a_element_grid=a_water_grid,  # Grid of the available positions for the element
                            a_atom_center=a_atom_center,  # Center of the element to place
                            a_sphere_shape=a_sphere_shape  # Spherical array of points
                        )
                        if b_placed:
                            update_dic(a_atom_center, f_max_score, gp.D_WATER_SCORING)
                            b_placed = False
                    # End if

                    b_compute_distance = False  # Stops computing distances for the next iteration

                # If the distance does not need to be calculated
                else:
                    a_atom_center = a_max_coord[
                        np.random.randint(len(a_max_coord))]  # Selects a random position for this element
                    b_placed = place_element_volume(  # Clears the volume around the placed element
                        a_structure_grid=o_structure.a_grid,  # Grid of the main structure
                        a_element_grid=a_water_grid,  # Grid of the available positions for the element
                        a_atom_center=a_atom_center,  # Center of the element to place
                        a_sphere_shape=a_sphere_shape  # Spherical array of points
                    )
                    if b_placed:
                        update_dic(a_atom_center, f_max_score, gp.D_WATER_SCORING)
                        b_placed = False
            # End if

            # If there is only one position for this score
            else:
                a_atom_center = a_max_coord[0]  # Selects the only position for this score
                b_placed = place_element_volume(  # Clears the volume around the placed element
                    a_structure_grid=o_structure.a_grid,  # Grid of the main structure
                    a_element_grid=a_water_grid,  # Grid of the available positions for the element
                    a_atom_center=a_atom_center,  # Center of the element to place
                    a_sphere_shape=a_sphere_shape  # Spherical array of points
                )
                if b_placed:
                    update_dic(a_atom_center, f_max_score, gp.D_WATER_SCORING)
                    b_placed = False
                b_compute_distance = True  # Asks for computing distances in the next iteration
                i_score_index += 1  # Continues with the next score
        # End if

        # If there is no element for this distance
        else:
            b_compute_distance = True  # Asks for computing distances in the next iteration
            i_repeat_count += 1  # Requests a refining step by actualizing the list of scores

            # If at least a repetition is needed
            if i_repeat_count > i_max_repeat:
                l_all_scores, a_element_positions = sort_element_score(
                    # Sorts the possibles scores in the element grid
                    a_element_grid=a_water_grid,  # The grid containing the element possible positions
                    i_element_code=i_water_code,  # The code of the element to place
                    f_min_score=d_parameters["f_min_water_score"]  # The minimal score for the placement
                )
                i_score_index = 0  # Resets the determination of available positions
                i_repeat_count = 0  # Resets the number of needed repetition

            # If no repetition are needed
            else:
                i_score_index += 1  # Continues with the next score
        # End if
    # End if
    # End while
    # Assign score



    del a_water_grid  # Frees some memory


# END STEP 6 ---------------------------------------- #

# ---------------------------------------------------------------------------- #


# Auxiliary functions -------------------------------------------------------- #

# Placing water molecules

def sort_element_score(a_element_grid, i_element_code, f_min_score):
    """
	Sorts each score value in the element grid
	"""

    a_element_positions = np.array(np.where(
        a_element_grid["element_symbol"] == i_element_code  # Retrieves the positions where the element is present
    )) #TODO use this (and the rest)
    l_all_scores = np.sort(
        list(set(  # Retrieves the score of each point containing the element, save it into a sorted set list
            a_element_grid["score"][
                a_element_positions[0],
                a_element_positions[1],
                a_element_positions[2]
            ].flatten().tolist()
        )))[::-1]
    l_all_scores = [x for x in l_all_scores if x >= f_min_score]  # Deletes the scores below the minimal score

    # If there is a negative value for the minimal score
    if f_min_score <= 0:

        # If there is null scores
        if 0 in l_all_scores:
            l_all_scores.remove(0)  # Removes the null scores

    return l_all_scores, a_element_positions  # Returns the list of sorted scores and the array of elements coordinates


def retrieve_element_coordinates(a_element_grid, a_element_positions, f_max_score):
    """
	Retrieves the list of coordinates of the points containing the maximum score
	"""
    a_element_coord = np.array(  # Retrieves in the element grid the elements positions with the maximum score
        a_element_grid["score"][
            a_element_positions[0],
            a_element_positions[1],
            a_element_positions[2]
        ] == f_max_score
    )
    a_element_coord = a_element_positions.T[a_element_coord]  # Formats the coordinates
    return a_element_coord  # Returns the array of coordinates of elements with the maximum score


# End function ------------------------------------------ #

def place_element_volume(a_structure_grid, a_element_grid, a_atom_center, a_sphere_shape):
    """
	Saves the position of the placed element and clear its neighbouring space
	:return: boolean telling if the molecule was placed or not
	"""
    # TODO : tester juste sur la partie 1 le random
    b_place = random_creation()
    if b_place is True:
        a_structure_grid[  # Saves the placed element with the structure
            a_atom_center[0],
            a_atom_center[1],
            a_atom_center[2]
        ] = cp.deepcopy(a_element_grid[
                            a_atom_center[0],
                            a_atom_center[1],
                            a_atom_center[2]
                        ])
    a_element_grid[  # Removes anything in the radius of the placed element in the element specific grid
        a_atom_center[0] + a_sphere_shape[0],
        a_atom_center[1] + a_sphere_shape[1],
        a_atom_center[2] + a_sphere_shape[2]
    ] = (0, 0, 0)
    return b_place


# End function ------------------------------------------ #

# Sphere formulas

def determine_sphere_manhattan(l_i_radius_range, f_element_radius):
    """
	Determines the array of points contained within the VdW sphere of the element
	"""
    a_sphere = np.array([(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range
                         if abs(x) + abs(y) + abs(z) <= f_element_radius]).astype(np.int32).T
    return a_sphere  # Returns the relative positions of points included in the sphere


# End function ------------------------------------------ #

def determine_sphere_minkowski(l_i_radius_range, f_element_radius):
    """
	Determines the array of points contained within the VdW sphere of the element
	"""
    a_sphere = np.array(
        [(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range]).astype(np.int32).T
    return a_sphere  # Returns the relative positions of points included in the sphere


# End function ------------------------------------------ #

def determine_sphere_euclidean(l_i_radius_range, f_element_radius):
    """
	Determines the array of points contained within the VdW sphere of the element
	"""
    a_sphere = np.array([(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range
                         if x ** 2 + y ** 2 + z ** 2 <= f_element_radius ** 2]).astype(np.int32).T
    return a_sphere  # Returns the relative positions of points included in the sphere


# End function ------------------------------------------ #

# Distance formulas

def compute_distance_manhattan(a_coords_a, a_coords_b):
    """
	Computes a taxicab distance between two arrays of points
	"""
    a_distance = np.sum(np.abs(a_coords_a - a_coords_b), axis=1)  # Computes the distance between each points
    return a_distance  # Returns the distance


# End function ------------------------------------------ #

def compute_distance_minkowski(a_coords_a, a_coords_b):
    """
	Computes a uniform distance between two arrays of points
	"""
    a_distance = np.max(np.abs(a_coords_a - a_coords_b), axis=1)  # Computes the distance between each points
    return a_distance  # Returns the distance


# End function ------------------------------------------ #

def compute_distance_euclidean(a_coords_a, a_coords_b):
    """
	Computes a spherical distance between two arrays of points
	"""
    a_distance = np.rint(np.sqrt(np.sum(
        np.power(a_coords_a - a_coords_b, 2),  # Computes the distance between each points
        axis=1
    ))).astype(np.uint16)
    return a_distance  # Returns the distance


# End function ------------------------------------------ #


def random_creation():
    """
	Decides randomly if a molecule will be placed or not
    :param f_probability: the probability of the molecule being placed
    between O.0 and 1.0.
    :return: Boolean
    """
    if gp.D_PARAMETERS_SOLUBILIZATION['b_use_randomax'] == False:
        return True
    else :
        f_probability = gp.D_PARAMETERS_SOLUBILIZATION['f_random_threshold']
        b_choice = None
        f_choice = rd.uniform(0, 1)
        if f_choice < f_probability:
            b_choice = True
        else:
            b_choice = False
        return b_choice

def update_dic(a_atom_center, f_score, d_to_update):
    """
    	Update a dictionary, assigning a score to a position.
        :param a_atom_center: array, containing the position of the molecule
        :param f_score: the score associated to the position
        :param d_to_update: the dictionary to update
    """
    i_posi_x = a_atom_center[0]
    i_posi_y = a_atom_center[1]
    i_posi_z = a_atom_center[2]
    key_dic = str(i_posi_x) + "_" + str(i_posi_y) + "_" + str(i_posi_z)
    d_to_update[key_dic] = f_score


# End function ------------------------------------------ #

# ---------------------------------------------------------------------------- #


# Reference ------------------------------------------------------------------ #

# Importation
# from place_water_itermax import place_water_itermax
# Iteratively places water molecules by decreasing order of score and updates the grid
# In : (d) parameters used for the solubilization, (o) system containing variables used for the solubilization
# In : (o) structure to solubilize
# Out : None

# Usage
# place_water_itermax(				# Place water molecules around the structure with the Itermax method
# 	d_parameters=d_parameters,		# Dictionary of parameters used for the placement of water molecules
# 	o_system=o_system,				# The system containing the structures parameters and functions
# 	o_structure=o_structure			# The structure to analyse
# )

# ---------------------------------------------------------------------------- #
