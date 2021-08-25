# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import multiprocessing as mp
	# Allows to work on multiple CPU threads
import numpy as np
	# Allows Numpy array manipulation
from sklearn.neighbors import KDTree
	# Allows the representation of the structure in a graph

# Parameters
from config import global_parameters as gp
	# Contains the global variabless

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def compute_position_scores(d_parameters, o_system, o_structure):
	"""
	Computes the score of each valid positions
	:param d_parameters: Dictionary of parameters to use for the point scoring
	:param o_system: The system containing the tasks and their parameters
	:param o_structure: The structure to score
	"""

	# STEP 0 : Preparing variables ---------------------- #
	i_water_code = gp.D_ELEMENT_NUMBER["OOW"]		# Retrieves the element code for water
	o_system.l_l_tasks = np.array(np.where(			# Retrieves the list of positions marked as available for water molecules
		o_structure.a_grid["element_symbol"] == i_water_code
	)).T
	o_system.l_l_tasks = (							# Retrieves the coordinates of the point to analyse and convert it in real space coordinates in Angstroms
		o_system.l_l_tasks.astype(np.float64)		# The coordinates in the grid
		- o_system.a_offset							# Updates with the system offset
	) * o_system.f_grid_spacing						# Converts the grid coordinates into Angstroms
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Creating the KDTree ---------------------- #
	l_l_atom_coord = []		# A list containing the coordinates of each atom of the structure

	# For each atom of the structure
	for a_atom in o_structure.a_atoms:
		l_l_atom_coord.append([a_atom[8], a_atom[9], a_atom[10]])		# Appends the coordinates of the atom

	a_atom_coord = np.array(l_l_atom_coord)		# Converts the atom coordinates into a numpy array

	o_structure.o_tree = KDTree(					# Generates a tree from the structure
		a_atom_coord,								# The atom coordinates to represent
		leaf_size=50000,							# Number of nodes in the tree
		metric=d_parameters["s_scoring_metric"]		# The metric to use to determine the neighbouring atoms
	)
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Defining settings ------------------------ #
	# If the user has defined a number of CPU
	if gp.D_PARAMETERS_GLOBAL["i_cpu_allocated"] is not None:
		i_cpu_count = gp.D_PARAMETERS_GLOBAL["i_cpu_allocated"]

	# If the user has not defined a number of CPU to use
	else:
		i_cpu_count = mp.cpu_count()		# Retrieves the amount of available CPU
	# END STEP 2 ---------------------------------------- #


	# STEP 3 : Running the tasks ------------------------ #
	s_text = "Scoring water for structure {} / {}".format(		# Defining the text to prompt in the progression bar
		o_system.i_current_structure + 1,
		len(o_system.l_o_structures)
	)
	o_system.setup_progress(		# Initializes the progression bar
		s_text=s_text
	)
	o_task_manager = mp.Pool(		# Creates a task manager
		processes=i_cpu_count		# Allocates a number of CPU to the task manager
	)
	l_l_results = o_task_manager.map(		# Launches the tasks
		func=compute_score_thread,				# The function to run on multiple CPU
		iterable=range(len(o_system.l_l_tasks)),		# Allocates a task index to each task
		chunksize=5000
	)
	o_task_manager.close()			# Closes the pool of tasks
	o_task_manager.join()			# Waits for the results
	o_task_manager.terminate()		# Kills the pool of tasks
	o_system.close_progress()		# Closes the progression bar
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : Saving the results ----------------------- #
	a_results = np.array(list(l_l_results)).T		# Extracts the scores from the results list
	a_results = a_results.astype(np.float32)		# Converts the results to float

	# Tries to retrieve the max score from the results
	try:
		f_max_score = np.max(a_results)					# Retrieves the maximal score

	# If there is no score
	except ValueError:
		f_max_score = 0		# Sets a default score
	a_results /= f_max_score						# Normalizes the scores

	a_positions = np.array(np.where(		# Retrieves the list of positions marked as available for water molecules
		o_structure.a_grid["element_symbol"] == i_water_code
	)).T

	l_x_coords = []		# List used to store the x coordinates
	l_y_coords = []		# List used to store the y coordinates
	l_z_coords = []		# List used to store the z coordinates

	# For each scored position
	for a_position in a_positions:
		l_x_coords.append(int(a_position[0]))		# Saves the X coordinate
		l_y_coords.append(int(a_position[1]))		# Saves the Y coordinate
		l_z_coords.append(int(a_position[2]))		# Saves the Z coordinate

	o_structure.a_grid["score"][		# Saves the scores in the grid
		l_x_coords,
		l_y_coords,
		l_z_coords
	] = a_results
	# END STEP 4 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

def compute_score_thread(i_position_index):
	"""
	Computes the score of a point around the structure
	:param i_position_index: The index of the point to analyse
	:return: The score of this position
	"""

	# Preparing variables
	d_parameters = gp.D_PARAMETERS_SOLUBILIZATION							# Loads the parameters used for the scoring
	o_system = gp.O_SYSTEM_SOLUBILIZATION									# Loads the system containing the structure
	o_structure = o_system.l_o_structures[o_system.i_current_structure]		# Loads the structure to work with

	a_atom_coords = [		# Retrieves the coordinates of the point to analyse
		(o_system.l_l_tasks[i_position_index])
	]

	# Searching neighbouring atoms
	a_i_index, a_f_distance = o_structure.o_tree.query_radius(		# Retrieves the index and distance of atom neighbouring this space point
		a_atom_coords,												# The coordinates of the point to analyse
		d_parameters["f_max_neighbor_distance"],					# The maximal search radius for neighbors
		return_distance=True,										# Returns the distance with each neighbour
		sort_results=True											# Sorts the neighbours by distances
	)

	# Scoring positions
	l_f_score = []						# Initializes the score
	i_skipped_atoms = 0					# The number of atoms ignored
	l_s_encountered_residues = []		# List of residues already encountered

	# For each neighbouring atom
	for i_neighbor in range(len(a_i_index[0])):

		# If only one atom per residue needs to be considered
		if d_parameters["b_scoring_per_residue"]:

			s_residue_key = "{}_{}".format(												# Creates a key
				o_structure.a_atoms["residue_serial"][a_i_index[0][i_neighbor]],		# Retrieves the serial number of the current residue
				o_structure.a_atoms["chain_id"][a_i_index[0][i_neighbor]]				# Retrieves the chain identifier of the current residue
			)

			# If the atom is part of an already considered residue
			if s_residue_key in l_s_encountered_residues:
				i_skipped_atoms += 1		# Keeps tracks of the skipped atoms
				continue					# Skip this neighbouring atom and process the next one

			# If the residue has not been encountered yet
			else:
				l_s_encountered_residues.append(s_residue_key)		# Keep tracks of the encountered atom
		# End if

		s_interaction = "OOW_{}_{}".format(
			o_structure.a_atoms["custom_type"][a_i_index[0][i_neighbor]],
			i_neighbor + 1
		)
		f_score_buffer = o_system.retrieve_nearest_score(		# Gets the score of the interaction
			s_interaction=s_interaction,						# The type of the interaction
			f_distance=a_f_distance[0][i_neighbor]				# The distance between the elements
		)

		# If an error has occurred
		if f_score_buffer == 0:
			i_skipped_atoms += 1		# Keeps tracks of the skipped atoms

		# If the scoring process worked fine
		else:
			l_f_score.append(f_score_buffer)		# Saves the score

		# If only the nearest neighbour needs to be considered
		if d_parameters["b_only_first_neighbor"]:

			# If a single interaction has been considered
			if i_neighbor - i_skipped_atoms + 1 == 1:
				break		# Stops the loop
		# End if

		# If enough neighbouring atoms have been considered
		if i_neighbor - i_skipped_atoms + 1 == d_parameters["i_max_neighbor_number"]:
			break		# Stops the loop
	# End for

	o_system.update_progress()		# Updates the progression bar

	# If the function needs to use mean scores and there is at least 1 score
	if d_parameters["b_mean_score"] and len(l_f_score) != 0:

		return sum(l_f_score) / len(l_f_score)		# Returns the position's mean score

	# If the function needs to use the summed scores
	else:
		return sum(l_f_score)		# Returns the sum of scores
# End function ------------------------------------------ #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from compute_position_scores import compute_position_scores
	# Computes the score of each possible position for water molecules
	# In : (d) parameters used for scoring, (o) the system holding the structures
	# In : (o) the structure containing the positions to score
	# Out : None

# Usage
# compute_position_scores(			# Computes the scores of each possible position for water molecules
# 	d_parameters=d_parameters,		# Dictionary of parameters used for the scoring of positions
# 	o_system=o_system,				# The system containing the structures parameters and functions
# 	o_structure=o_structure			# The structure to analyse
# )

# ---------------------------------------------------------------------------- #
