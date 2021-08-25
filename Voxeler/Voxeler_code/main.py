#!/usr/bin/python3

# Information ---------------------------------------------------------------- #
# Comparing Structures and Solubilizing Compounds CSaSC
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : June 2020
# Updated : February 2021
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import os
	# Allows file system operations
import time
	# Enables time manipulation

# Parameters
from config import global_parameters as gp
	# Contains the global variables

# General library
from lib.extract_valid_parameters import extract_valid_parameters
	# Extract parameters from a file and checks the validity of the values
	# In : (p) path to the parameters, (d) dictionary to complete,
	# In : (d) expected parameters, type and default values
	# Out : None

# Specific modules
from src.launch_structure_comparison import launch_structure_comparison
	# Manages the comparison of structures
	# In : None
	# Out : None
from src.launch_structure_solubilization import launch_structure_solubilization
#from src.save_density import create_dictionary


# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def main():
	"""
	TODO
	:param :
	:return:
	"""

	# STEP 0 : Loading parameters ----------------------- #
	gp.init()							# Initializes the global parameters variables
	gp.loads_default_parameters()		# Loads the base parameters

	extract_valid_parameters(										# Retrieves the program global parameters
		p_file=gp.D_PARAMETERS_GLOBAL["p_global_parameters"],		# The path to the file containing the parameters
		d_parameters=gp.D_PARAMETERS_GLOBAL,						# The dictionary container for the parameters
		d_expected_parameters=gp.D_EXPECTED_PARAMETERS_GLOBAL		# The dictionary guiding the extraction
	)
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Structure comparison --------------------- #
	# If a structure comparison must be done
	if gp.D_PARAMETERS_GLOBAL["b_run_comparison"]:

		t_start = time.time()												# Gets the actual time
		extract_valid_parameters(											# Retrieves the program structure comparison parameters
			p_file=gp.D_PARAMETERS_GLOBAL["p_comparison_parameters"],		# The path to the file containing the parameters
			d_parameters=gp.D_PARAMETERS_COMPARISON,						# The dictionary container for the parameters
			d_expected_parameters=gp.D_EXPECTED_PARAMETERS_COMPARISON		# The dictionary guiding the extraction
		)
		launch_structure_comparison()		# Manages the comparison of structures
		print("Comparison done in {:.1f} seconds".format(time.time() - t_start))
	# END STEP 1 ---------------------------------------- #

	# STEP 2 : Structure solubilization ----------------- #
	# If a structure solubilization must be done
	if gp.D_PARAMETERS_GLOBAL["b_run_solubilization"]:

		t_start = time.time()													# Gets the actual time
		extract_valid_parameters(												# Retrieves the program structure solubilization parameters
			p_file=gp.D_PARAMETERS_GLOBAL["p_solubilization_parameters"],		# The path to the file containing the parameters
			d_parameters=gp.D_PARAMETERS_SOLUBILIZATION,						# The dictionary container for the parameters
			d_expected_parameters=gp.D_EXPECTED_PARAMETERS_SOLUBILIZATION		# The dictionary guiding the extraction
		)

		launch_structure_solubilization()		# Manages the comparison of structures
		print("Solubilization done in {:.1f} seconds".format(time.time() - t_start))
	# END STEP 2 ---------------------------------------- #

	# STEP 3 : ------------------------------------------ #
	# END STEP 3 ---------------------------------------- #

	# STEP 4 : ------------------------------------------ #
	# END STEP 4 ---------------------------------------- #

	# STEP 5 : ------------------------------------------ #
	# END STEP 5 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

if __name__ == "__main__":
	main()		# Launches the main function

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from main import main
	#
	# In :
	# Out :

# Usage
# main(		#
# = ,		#
# )

# ---------------------------------------------------------------------------- #
