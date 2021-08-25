# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
# Program external resources
# Parameters
from config import global_parameters as gp
	# Contains the global variables

# Classes
from cla.system_solubilization import SystemSolubilization

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

# Specific modules
from src.solubilize_pdb_structure import solubilize_pdb_structure
from src.convert_dic_pdb import thresh_positions
from src.convert_dic_pdb import convert_dic_pdb


# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def launch_structure_solubilization():
	"""
	TODO
	:param :
	:return:
	"""

	# STEP 0 : Preparing variables ---------------------- #
	l_s_logs = []											# A list for log
	gp.O_SYSTEM_SOLUBILIZATION = SystemSolubilization()			# Sets a system object as a global variable
	gp.O_SYSTEM_SOLUBILIZATION.initialize_system(				# Initializes system fields
		d_parameters=gp.D_PARAMETERS_SOLUBILIZATION			# The parameters used for the system initial configuration
	)
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Find PDB files --------------------------- #
	l_p_input_pdb = retrieve_specific_files(							# Retrieves PDB paths
		p_directory=gp.D_PARAMETERS_SOLUBILIZATION["p_input_pdb"],		# Path to the input directory
		s_pattern="*.pdb",												# Pattern to match within the directories
		b_recursive=True,												# Also searches in the subdirectories
		i_min_match=1,													# Minimum number of files to retrieve
		i_max_match=9999												# Maximum number of files to retrieve
	)
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Extracts PDB structures ------------------ #
	d_parsing_parameters = {																					# Dictionary of parameters required for the PDB parsing
		"b_discard_atom": gp.D_PARAMETERS_SOLUBILIZATION["b_discard_atom"],											# Discards atoms lines
		"b_discard_hetatm": gp.D_PARAMETERS_SOLUBILIZATION["b_discard_hetatm"],										# Discards hetero atoms lines
		"b_discard_hydrogen": gp.D_PARAMETERS_SOLUBILIZATION["b_discard_hydrogen"],									# Discards Hydrogen atoms
		"b_discard_water": gp.D_PARAMETERS_SOLUBILIZATION["b_discard_water"],										# Discards water molecules
		"b_discard_alternative": gp.D_PARAMETERS_SOLUBILIZATION["b_discard_alternative"],							# Discards alternative positions
		"l_c_chain_white": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_c_chain_white"]],			# List of chains to keep, discards others
		"l_c_chain_black": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_c_chain_black"]],			# List of chains to discard
		"l_s_residue_white": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_s_residue_white"]],		# List of residues to keep, discards others
		"l_s_residue_black": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_s_residue_black"]],		# List of residues to discard
		"l_i_residue_white": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_i_residue_white"]],		# List of residues ID to keep, discards others
		"l_i_residue_black": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_i_residue_black"]],		# List of residues ID to discard
		"l_s_atom_white": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_s_atom_white"]],				# List of atom type to keep, discards others
		"l_s_atom_black": [item.upper() for item in gp.D_PARAMETERS_SOLUBILIZATION["l_s_atom_black"]],				# List of atom type to discard
	}

	# For each PDB file to parse
	for p_pdb in l_p_input_pdb:
		o_structure = parse_pdb_file(			# Extracts a PDB structure into an object
			p_file=p_pdb,						# The PDB file to extract
			d_filters=d_parsing_parameters		# The parsing filters to apply
		)
		gp.O_SYSTEM_SOLUBILIZATION.l_o_structures.append(o_structure)		# Registers the structure in the system

	gp.O_SYSTEM_SOLUBILIZATION.actualize_properties()		# Actualizes the system properties depending of the structures
	# END STEP 2 ---------------------------------------- #


	# STEP 3 : ------------------------------------------ #
	# For each structure to solubilize
	for o_structure in gp.O_SYSTEM_SOLUBILIZATION.l_o_structures:
		# Running several times solub :

		if gp.D_PARAMETERS_SOLUBILIZATION["b_use_randomax"] == False:
			i_n_launch = 1
		else:
			i_n_launch = gp.D_PARAMETERS_SOLUBILIZATION['i_launch_number']
		curent_launch = 0
		while curent_launch < i_n_launch:
			print("Launch {}/{}".format((curent_launch +1), i_n_launch))
			gp.O_SYSTEM_SOLUBILIZATION.generate_grid(			# Generates a grid for the structure
				o_structure=o_structure,						# The structure to include in the grid
				d_parameters=gp.D_PARAMETERS_SOLUBILIZATION		# The parameters used for the structure incorporation
			)
			solubilize_pdb_structure(								# Places water molecules around the structure
				d_parameters=gp.D_PARAMETERS_SOLUBILIZATION,		# The parameters used for the solubilization
				o_system=gp.O_SYSTEM_SOLUBILIZATION,				# The system containing all the elements used for the solubilization
				o_structure=o_structure								# The structure to solubilize
			)
			curent_launch += 1
		#TODO formater calls
		if gp.D_PARAMETERS_SOLUBILIZATION["b_use_randomax"] is True:
			a_density = thresh_positions(						# Select the position with a sufficient number of occurences
				gp.D_WATER_POSITION
			)
			convert_dic_pdb(									# Write the positions in a pdb file
				a_density,
				o_structure,
				gp.D_PARAMETERS_SOLUBILIZATION
			)
		gp.O_SYSTEM_SOLUBILIZATION.i_current_structure += 1		# Allows to work on the next structure
		del o_structure		# Deletes the structure and frees memory
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : ------------------------------------------ #
	# END STEP 4 ---------------------------------------- #


	# STEP 5 : ------------------------------------------ #
	# END STEP 5 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #
# End function ------------------------------------------ #
# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from launch_structure_solubilization import launch_structure_solubilization
	#
	# In :
	# Out :

# Usage
# launch_structure_solubilization(		#
# = ,		#
# )

# ---------------------------------------------------------------------------- #
