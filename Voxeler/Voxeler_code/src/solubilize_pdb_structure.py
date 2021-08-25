# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
# Program external resources
# Parameters
# Classes
# General library
# Specific modules
from src.determine_structure_surface import determine_structure_surface
	# Identifies neighbouring positions available for the placement of water molecules
	# In : (d) parameters, (o) system containing all the information,
	# In : (o) the structure to solubilize
	# Out : None
from src.compute_position_scores import compute_position_scores
	# Computes the score of each possible position for water molecules
	# In : (d) parameters used for scoring, (o) the system holding the structures
	# In : (o) the structure containing the positions to score
	# Out : None
from src.place_water_itermax import place_water_itermax
	# Iteratively places water molecules by decreasing order of score and updates the grid
	# In : (d) parameters used for the solubilization, (o) system containing variables used for the solubilization
	# In : (o) structure to solubilize
	# Out : None
from src.convert_grid_pdb import convert_grid_pdb

from config import global_parameters as gp


# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def solubilize_pdb_structure(d_parameters, o_system, o_structure):
	"""
	TODO
	:param d_parameters:
	:param o_system:
	:param o_structure:
	:return:
	"""

	# STEP 0 : Preparing variables ---------------------- #
	d_solubilization_functions = {
		"ITERMAX": solubilization_itermax,
	}
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Solubilizing the structure --------------- #
	#TODO : try except? 
	d_solubilization_functions[
		d_parameters["s_solubilization_method"].upper()
	](
		d_parameters=d_parameters,
		o_system=o_system,
		o_structure=o_structure
	)
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Determining best places for water -------- #

	# END STEP 2 ---------------------------------------- #


	# STEP 3 : ------------------------------------------ #
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : ------------------------------------------ #
	# END STEP 4 ---------------------------------------- #


	# STEP 5 : ------------------------------------------ #
	# END STEP 5 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

def solubilization_itermax(d_parameters, o_system, o_structure):
	"""
	"""

	# Preparing variables

	#
	determine_structure_surface(		# Identifies available positions for water molecules around the structure
		d_parameters=d_parameters,		# The solubilization parameters
		o_system=o_system,				# The system managing the solubilization
		o_structure=o_structure			# The structure to solubilize
	)
	compute_position_scores(			# Computes the scores of each possible position for water molecules
		d_parameters=d_parameters,		# Dictionary of parameters used for the scoring of positions
		o_system=o_system,				# The system containing the structures parameters and functions
		o_structure=o_structure			# The structure to analyse
	)
	place_water_itermax(				# Place water molecules around the structure with the Itermax method
		d_parameters=d_parameters,		# Dictionary of parameters used for the placement of water molecules
		o_system=o_system,				# The system containing the structures parameters and functions
		o_structure=o_structure			# The structure to analyse
	)
	convert_grid_pdb(
		d_parameters=d_parameters,		# Dictionary of parameters used for saving the solubilized structure
		o_system=o_system,
		o_structure=o_structure
	)


# End function ------------------------------------------ #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from solubilize_pdb_structure import solubilize_pdb_structure
	#
	# In :
	# Out :

# Usage
# solubilize_pdb_structure(		#
# = ,		#
# )

# ---------------------------------------------------------------------------- #
