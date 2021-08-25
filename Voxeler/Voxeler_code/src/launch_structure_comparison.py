# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : June 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
	# Allows Numpy array manipulation
import time
	# Enables time manipulation
import psutil
	# Allows the gathering of system information
import multiprocessing as mp
	# Allows the execution of code on multiple threads

# Parameters
from config import global_parameters as gp
	# Contains the global variables

# Classes
from cla.system_comparison import SystemComparison
	# Singleton object containing every structure to be compared
from cla.tree_plot import TreePlot
	# Generates, renders and saves trees

# General library
from lib.retrieve_specific_files import retrieve_specific_files
	# Retrieves files path, recursively or not, matching a specific pattern, or not
	# In : (p) directory to retrieve files from, (s) pattern to match,
	# In : (b) if the search needs to be recursive, (i) minimum number of match,
	# In : (i) maximum number of match
	# Out : (l(p)) a list of the file paths founds
from lib.parse_pdb_file import parse_pdb_file
	# Extracts a PDB structure from a file and applies filters
	# In : (p) PDB file to extract, (d) parsing filters to apply
	# Out : (o) the object containing the structure
from lib.terminate_program_process import terminate_program_process
	# Stops the program and prints content
	# In : (l(s)) content to prompt
	# Out : None

# Specific modules
from src.compute_grid_similarity import compute_grid_similarity
	# Computes the similarity between two grids
	# In : (d) comparison parameters, (o) the first structure to compare,
	# In : (o) the second structure to compare
	# Out : (f) the percentage of similarity

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def launch_structure_comparison():
	"""
	Manages the comparison of structures
		Retrieves the PDB input files
		Extracts the PDB structures
		Compares the structures
		Generates a tree
	"""

	# STEP 0 : Preparing variables ---------------------- #
	l_s_logs = []									# A list for log
	gp.O_SYSTEM_COMPARISON = SystemComparison()		# Sets a system object as a global variable
	gp.O_SYSTEM_COMPARISON.initialize_system(		# Initializes system fields
		d_parameters=gp.D_PARAMETERS_COMPARISON		# The parameters used for the system initial configuration
	)
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Find PDB files --------------------------- #
	l_p_input_pdb = retrieve_specific_files(						# Retrieves PDB paths
		p_directory=gp.D_PARAMETERS_COMPARISON["p_input_pdb"],		# Path to the input directory
		s_pattern="*.pdb",											# Pattern to match within the directories
		b_recursive=True,											# Also searches in the subdirectories
		i_min_match=1,												# Minimum number of files to retrieve
		i_max_match=9999											# Maximum number of files to retrieve
	)
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Extracts PDB structures ------------------ #
	d_parsing_parameters = {																					# Dictionary of parameters required for the PDB parsing
		"b_discard_atom": gp.D_PARAMETERS_COMPARISON["b_discard_atom"],											# Discards atoms lines
		"b_discard_hetatm": gp.D_PARAMETERS_COMPARISON["b_discard_hetatm"],										# Discards hetero atoms lines
		"b_discard_hydrogen": gp.D_PARAMETERS_COMPARISON["b_discard_hydrogen"],									# Discards Hydrogen atoms
		"b_discard_water": gp.D_PARAMETERS_COMPARISON["b_discard_water"],										# Discards water molecules
		"b_discard_alternative": gp.D_PARAMETERS_COMPARISON["b_discard_alternative"],							# Discards alternative positions
		"l_c_chain_white": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_c_chain_white"]],			# List of chains to keep, discards others
		"l_c_chain_black": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_c_chain_black"]],			# List of chains to discard
		"l_s_residue_white": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_s_residue_white"]],		# List of residues to keep, discards others
		"l_s_residue_black": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_s_residue_black"]],		# List of residues to discard
		"l_i_residue_white": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_i_residue_white"]],		# List of residues ID to keep, discards others
		"l_i_residue_black": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_i_residue_black"]],		# List of residues ID to discard
		"l_s_atom_white": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_s_atom_white"]],				# List of atom type to keep, discards others
		"l_s_atom_black": [item.upper() for item in gp.D_PARAMETERS_COMPARISON["l_s_atom_black"]],				# List of atom type to discard
	}

	# For each PDB file to parse
	for p_pdb in l_p_input_pdb:
		o_structure = parse_pdb_file(			# Extracts a PDB structure into an object
			p_file=p_pdb,						# The PDB file to extract
			d_filters=d_parsing_parameters		# The parsing filters to apply
		)
		gp.O_SYSTEM_COMPARISON.l_o_structures.append(o_structure)		# Registers the structure in the system

	gp.O_SYSTEM_COMPARISON.actualize_properties()		# Actualizes the system properties depending of the structures
	# END STEP 2 ---------------------------------------- #


	# STEP 3 : Security check --------------------------- #
	# If there is not enough file
	if len(gp.O_SYSTEM_COMPARISON.l_o_structures) < 2:
		l_s_logs.append("ERROR : There is not enough valid PDB structure to run a comparison, at least 2 are required")		# Defines the error message
		terminate_program_process(		# Stops the program
			l_s_content=l_s_logs		# Content to save to logs
		)
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : Determining the normalisation method ----- #
	s_normalisation = gp.D_PARAMETERS_COMPARISON["s_comparison_normalisation"].upper()		# Retrieves the normalisation method to use

	# If the program needs to retrieve the number of non empty points for each structure
	if s_normalisation == "GLOBAL_MIN" or s_normalisation == "GLOBAL_MAX":

		l_structure_size = []		# List containing the size of each structure

		# For each loaded structure
		for o_structure in gp.O_SYSTEM_COMPARISON.l_o_structures:
			gp.O_SYSTEM_COMPARISON.generate_grid(			# Creates a grid, loads the structure into it, place VdW volumes
				o_structure=o_structure,					# The structure to place into a grid
				d_parameters=gp.D_PARAMETERS_COMPARISON		# Parameters used for the VdW volumes generation
			)
			l_structure_size.append(np.count_nonzero(o_structure.a_grid))		# Counts the number of points containing atoms
			o_structure.a_grid = None		# Frees some memory
		# End for

	# If the similarity needs to be based on the minimal number of empty points
	if s_normalisation == "GLOBAL_MIN":
		gp.D_PARAMETERS_COMPARISON["i_atom_total"] = min(l_structure_size)		# Retrieves the minimal number of non empty points

	if s_normalisation == "GLOBAL_MAX":
		gp.D_PARAMETERS_COMPARISON["i_atom_total"] = max(l_structure_size)		# Retrieves the maximal number of non empty points

	if s_normalisation == "MAX" or s_normalisation == "MIN":
		pass

	else:
		l_s_logs.append("ERROR : Unknown normalisation method. Known methods are 'Min' and 'Max'")
		terminate_program_process(		# Stops the program
			l_s_content=l_s_logs		# Content to save in the logs
		)
	# End if
	# END STEP 4 ---------------------------------------- #


	# STEP 5 : Running the comparison ------------------- #
	i_total_tasks = int(										# Computes the number of combination to process
		len(gp.O_SYSTEM_COMPARISON.l_o_structures) *
		(len(gp.O_SYSTEM_COMPARISON.l_o_structures) - 1)
	) / 2
	l_i_indexes = range(len(									# Creates an iterable over the structures to process
		gp.O_SYSTEM_COMPARISON.l_o_structures
	))

	# For each structure to compare
	for i_first_index in l_i_indexes:

		# For each other structure to compare
		for i_second_index in l_i_indexes[i_first_index+1:]:
			gp.O_SYSTEM_COMPARISON.l_l_tasks.append([i_first_index, i_second_index])		# Saves the parameters of the each task
	# End for

	# If the user has defined a number of CPU
	if gp.D_PARAMETERS_GLOBAL["i_cpu_allocated"] is not None:
		i_cpu_count = gp.D_PARAMETERS_GLOBAL["i_cpu_allocated"]

	# If the user has not defined a number of CPU to use
	else:
		i_cpu_count = mp.cpu_count()		# Retrieves the amount of available CPU

	# Running the tasks
	gp.O_SYSTEM_COMPARISON.setup_progress()		# Initializes the progression bar
	o_task_manager = mp.Pool(					# Creates a task manager
		processes=i_cpu_count					# Allocates a number of CPU to the task manager
	)
	l_l_results = o_task_manager.map(								# Launches the tasks
		func=compare_grid_multithreaded,							# The function to run on multiple CPU
		iterable=range(len(gp.O_SYSTEM_COMPARISON.l_l_tasks)),		# Allocates a task index to each task
		chunksize=5000
	)
	o_task_manager.close()									# Closes the pool of tasks
	o_task_manager.join()									# Joins the results
	o_task_manager.terminate()								# Kills the pool of tasks, reallocates resources
	gp.O_SYSTEM_COMPARISON.close_progress()					# Closes the progression bar
	del gp.O_SYSTEM_COMPARISON								# Deletes the object and frees memory
	# END STEP 5 ---------------------------------------- #


	# STEP 6 : Generating the tree ---------------------- #
	o_tree = TreePlot()		# Creates a tree object

	# For each result to save
	for l_result in l_l_results:
		o_tree.add_score(
			s_first_node=l_result[0],		# The name of the first node to save
			s_second_node=l_result[1],		# The name of the second node to save
			f_score=l_result[2]				# The score to save
		)
	o_tree.generate_tree(							# Computes all the similarity scores and generates a tree
		d_parameters=gp.D_PARAMETERS_COMPARISON		# The parameters used for the tree generation
	)

	# If the user requested to display the tree
	if gp.D_PARAMETERS_COMPARISON["b_show_tree"]:
		o_tree.show_tree()		# Displays the tree

	o_tree.save_tree(								# Saves the tree alongside the comparison parameters
		d_parameters=gp.D_PARAMETERS_COMPARISON		# The parameters used for saving the tree
	)
	del o_tree										# Deletes the tree and frees memory
	# END STEP 6 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

def compare_grid_multithreaded(i_task):
	"""
	Generates and compares grids of atoms
	:param i_task: The index of the task to process
	:return: The score of similarity between the grids
	"""

	l_task = gp.O_SYSTEM_COMPARISON.l_l_tasks[i_task]							# Loads the list of structure indexes to use
	o_first_structure = gp.O_SYSTEM_COMPARISON.l_o_structures[l_task[0]]		# Loads the first structure
	o_second_structure = gp.O_SYSTEM_COMPARISON.l_o_structures[l_task[1]]		# Loads the second structure

	# If the program is not allowed to use all the available memory
	if gp.D_PARAMETERS_GLOBAL["f_memory_allocated"] > 0:

		# While there is not enough memory to process the comparison
		while (psutil.virtual_memory().used / 1073741824) - 2.0 > gp.D_PARAMETERS_GLOBAL["f_memory_allocated"]:
			time.sleep(1)		# Waits until there is available memory

	# Generating grids
	gp.O_SYSTEM_COMPARISON.generate_grid(			# Creates a grid, loads the structure into it, place VdW volumes
		o_structure=o_second_structure,				# The structure to place into a grid
		d_parameters=gp.D_PARAMETERS_COMPARISON		# Parameters used for the VdW volumes generation
	)
	gp.O_SYSTEM_COMPARISON.generate_grid(			# Creates a grid, loads the structure into it, place VdW volumes
		o_structure=o_second_structure,				# The structure to place into a grid
		d_parameters=gp.D_PARAMETERS_COMPARISON		# Parameters used for the VdW volumes generation
	)

	# Computing the grid similarities
	f_similarity = compute_grid_similarity(				# Computes the similarity between two structures
		d_parameters=gp.D_PARAMETERS_COMPARISON,		# Comparison parameters
		o_first_structure=o_first_structure,			# The first structure to compare
		o_second_structure=o_second_structure,			# The second structure to compare
	)

	# Freeing memory
	o_first_structure.delete_grid()			# Frees some memory some memory
	o_second_structure.delete_grid()		# Frees some memory some memory

	# Display
	gp.O_SYSTEM_COMPARISON.update_progress()		# Updates the progression bar

	return [o_first_structure.s_name, o_second_structure.s_name, f_similarity]		# Returning the comparison results
# End function ------------------------------------------ #
# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from launch_structure_comparison import launch_structure_comparison
	# Manages the comparison of structures
	# In : None
	# Out : None

# Usage
# launch_structure_comparison()		# Manages the comparison of structures

# ---------------------------------------------------------------------------- #
