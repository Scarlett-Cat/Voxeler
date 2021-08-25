# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
	# Allows Numpy array manipulation
import copy
	# Allows true copy of elements

# Parameters
import config.global_parameters as gp
	# Contains the global variables

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def determine_structure_surface(d_parameters, o_system, o_structure):
	"""
	Identifies neighbouring positions available for the placement of water molecules
	:param d_parameters: Parameters used for determination of available positions for water molecules
	:param o_system: The system containing the processed information needed to solubilize the structure
	:param o_structure: The structure to solubilize
	"""

	# STEP 0 : Preparing variables ---------------------- #
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Verifying the radius validity ------------ #
	f_min_radius = d_parameters["l_f_solubilization_radius"][0]		# Loads the maximal solubilization radius
	f_max_radius = d_parameters["l_f_solubilization_radius"][1]		# Loads the maximal solubilization radius
	i_water_code = gp.D_ELEMENT_NUMBER["OOW"]						# Retrieves the element code for water

	# If the minimal radius is not valid
	if f_min_radius < 0.1:
		f_min_radius = 0.1		# Sets the minimal radius to a minimal value

	# If the two radius are equals
	if f_min_radius == f_max_radius:
		f_max_radius *= 2		# Doubles the maximal radius value

	i_min_radius = int(f_min_radius / d_parameters["f_grid_spacing"])		# Converts the radius from angstroms to a grid distance
	i_max_radius = int(f_max_radius / d_parameters["f_grid_spacing"])		# Converts the radius from angstroms to a grid distance
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : ------------------------------------------ #
	# TODO Pockets
	# END STEP 2 ---------------------------------------- #


	# STEP 3 : Determining surface positions ------------ #
	o_water_structure = copy.deepcopy(o_structure)		# Copies the structure to solubilize in order to work with water
	o_system.generate_vdw_spheres(					# Extends the VdW volume
		o_structure=o_water_structure,				# The structure to solubilize
		d_parameters=d_parameters,					# The parameters used for the extension of VdW radius around the structure
		i_solubilization_radius=i_max_radius,		# The radius used for the VdW extension
	)
	o_system.generate_vdw_spheres(					# Removes the inner VdW volume
		o_structure=o_water_structure,				# The structure to solubilize
		d_parameters=d_parameters,					# The parameters used for the deletion of VdW radius in the structure
		i_solubilization_radius=i_min_radius,		# The radius used for the VdW deletion
		b_remove_volume=True						# Indicates the deletion of the selected volume
	)
	a_available_positions = np.array(np.where(		# Retrieves the number of non-empty points
		o_water_structure.a_grid["element_symbol"] > 0
	))
	del o_water_structure		# Frees memory
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : Saving available positions --------------- #
	# If there is any available position for water molecules
	if a_available_positions.any():

		o_structure.a_grid[					# Saves the available positions into the structure's grid
			a_available_positions[0],		# X coordinates of the available positions
			a_available_positions[1],		# Y coordinates of the available positions
			a_available_positions[2]		# Z coordinates of the available positions
		] = (
			i_water_code,		# Element code
			0,					# Position serial number
			0					# Position score
		)
	# END STEP 4 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from determine_structure_surface import determine_structure_surface
	# Identifies neighbouring positions available for the placement of water molecules
	# In : (d) parameters, (o) system containing all the information,
	# In : (o) the structure to solubilize
	# Out : None

# Usage
# 	determine_structure_surface(		# Identifies available positions for water molecules around the structure
# 		d_parameters=d_parameters,		# The solubilization parameters
# 		o_system=o_system,				# The system managing the solubilization
# 		o_structure=o_structure			# The structure to solubilize
# 	)

# ---------------------------------------------------------------------------- #
