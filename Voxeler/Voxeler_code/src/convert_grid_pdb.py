# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated : February 2021
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
	# Allows Numpy array manipulation

# Program external resources

# Parameters
from config import global_parameters as gp
	# Contains the global variabless

# Classes
# General library
# Specific modules

#from Voxeler.Voxeler_code.src.save_density import create_dictionary
from Voxeler.Voxeler_code.src.save_density import update_dictionary

from Voxeler.Voxeler_code.src.place_water_itermax import update_dic

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def convert_grid_pdb(d_parameters, o_system, o_structure):
	"""
	TODO
	:param :
	:return:
	"""

	# STEP 0 : Preparing variables ---------------------- #
	p_output = "{}/{}_{}.pdb".format(							# Generates the name of the output file
		d_parameters["p_output_solubilization"],				# Path to the output file
		d_parameters["s_solubilization_method"].lower(),		# The method used for the placement of water molecules
		o_structure.s_name										# The name of the structure
	)
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : ------------------------------------------ #
	# Pockets
	# TODO : pockets
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : ------------------------------------------ #

	i_water_code = gp.D_ELEMENT_NUMBER["OOW"]		# Retrieves the element code for water
	a_water_molecules = np.array(np.where(			# Retrieves the list of positions marked as available for water molecules
		o_structure.a_grid["element_symbol"] == i_water_code
	)).T

	l_placed_scores = []
	for i_obj in a_water_molecules:					# Retrieves the score assigned with each water position
		s_key = str(i_obj[0]) + "_" + str(i_obj[1]) + "_" + str(i_obj[2])
		f_score = gp.D_WATER_SCORING[s_key]
		l_placed_scores.append(f_score)				# Add the score in a list

	a_water_molecules = (							# Retrieves the coordinates of the point to analyse and convert it in real space coordinates in Angstroms
		a_water_molecules.astype(np.float64)		# The coordinates in the grid
		- o_system.a_offset							# Updates with the system offset
	) * o_system.f_grid_spacing						# Converts the grid coordinates into Angstroms


	a_water_score = np.array(l_placed_scores)		# Add the list of scores in the array
	a_water_score_shp = a_water_score.reshape(len(l_placed_scores), 1)

	a_water_molecules = np.append(a_water_molecules, a_water_score_shp, axis=1)
	for a_coord in a_water_molecules:				# Updates the dictionary associating pdb coord and score
		update_dic(
			a_coord,
			a_coord[3],
			gp.D_PDB_SCORING
		)

	if gp.D_PARAMETERS_SOLUBILIZATION["b_use_randomax"] is True:
		update_dictionary(							# Update the dictionary of positions and occurrences
			a_water_molecules,
			gp.D_WATER_POSITION
		)

	l_water_molecules = []

	for a_water in a_water_molecules:
		l_line_buffer = []
	
		l_line_buffer.append("HETATM")						# The element type
		l_line_buffer.append(65535 - i_water_code)		# The code of the water molecules
		l_line_buffer.append("OOW")							# The code of the water molecule
		l_line_buffer.append("")							# Possible alternative location
		l_line_buffer.append("HOH")							# Name of the residue
		l_line_buffer.append("")							# Water molecules are not part of a chain
		l_line_buffer.append(9999 - i_water_code)		# Serial number of the residue
		l_line_buffer.append("")							# Code for the insertion of a residue
		l_line_buffer.append(a_water[0])					# The occupancy of this atom at this position
		l_line_buffer.append(a_water[1])					# The occupancy of this atom at this position
		l_line_buffer.append(a_water[2])					# The occupancy of this atom at this position
		l_line_buffer.append(1.0)							# The occupancy of this atom at this position
		l_line_buffer.append(a_water[3])					# The score of this atome at this position
		l_line_buffer.append("O")							# The symbol of the element
		l_line_buffer.append(0)								# No charge for water

		l_water_molecules.append(l_line_buffer)
	# END STEP 2 ---------------------------------------- #

	# a_atoms = np.array(l_atoms).T
	# print(len(a_atoms))

	# print(len(o_structure.a_atoms["coord_x"]))
	# print(len(o_structure.a_atoms["grid_x"]))
	# o_structure.a_atoms["coord_x"] = (o_structure.a_atoms["grid_x"].astype(np.float64) - o_system.a_offset[0]) * o_system.f_grid_spacing
	# o_structure.a_atoms["coord_y"] = (o_structure.a_atoms["grid_y"].astype(np.float64) - o_system.a_offset[1]) * o_system.f_grid_spacing
	# o_structure.a_atoms["coord_z"] = (o_structure.a_atoms["grid_z"].astype(np.float64) - o_system.a_offset[2]) * o_system.f_grid_spacing

	# STEP 3 : ------------------------------------------ #
	write_pdb_output(
		o_structure=o_structure,
		l_water_molecules=l_water_molecules,
		p_output=p_output
	)
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : ------------------------------------------ #
	# END STEP 4 ---------------------------------------- #

	f_total_score = sum(l_placed_scores)/len(l_placed_scores)
	print("Solubilization mean score : {:.2f}".format(f_total_score))

	# STEP 5 : ------------------------------------------ #
	# END STEP 5 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

def write_pdb_output(o_structure, l_water_molecules, p_output):
	"""
	"""

	# Opening the output file
	with open(p_output, "w") as output_file:

		# For each atom line to write
		for a_atom in o_structure.a_atoms:
			s_output_line = convert_numpy_pdb(
				x_atom=a_atom
			)
			output_file.write(s_output_line + '\n')		# Writes the line containing the atom

		# For each water molecule to save
		for l_atom in l_water_molecules:

			s_output_line = convert_numpy_pdb(
				x_atom=l_atom
			)
			output_file.write(s_output_line + '\n')		# Writes the line containing the atom
# End function ------------------------------------------ #

def convert_numpy_pdb(x_atom):
	"""
	"""
	# If the a_atom variable is a list
	if isinstance(x_atom, list):

		s_output_line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
			str(x_atom[0]),		# ATOM or HETATM
			int(x_atom[1]),		# Atom serial number
			str(x_atom[2]),		# Atom name
			str(x_atom[3]),		# Alternate location indicator
			str(x_atom[4]),		# Residue name
			str(x_atom[5]),		# Chain identifier
			int(x_atom[6]),		# Residue sequence number
			str(x_atom[7]),		# Code for insertion of residues
			float(x_atom[8]),		# Orthogonal coordinates for X in Angstroms
			float(x_atom[9]),		# Orthogonal coordinates for Y in Angstroms
			float(x_atom[10]),		# Orthogonal coordinates for Z in Angstroms
			float(x_atom[11]),		# Occupancy
			float(x_atom[12]),		# Temperature factor
			str(x_atom[13]),		# Element symbol
			str(x_atom[14]),		# Charge on the atom
		)

	# If the a_atom variable is an array
	elif isinstance(x_atom, np.void):

		s_output_line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
			str(x_atom["element_type"]),				# ATOM or HETATM
			int(x_atom["atom_serial"]),				# Atom serial number
			str(x_atom["atom_name"]),				# Atom name
			str(x_atom["alternative_location"]),		# Alternate location indicator
			str(x_atom["residue_name"]),				# Residue name
			str(x_atom["chain_id"]),					# Chain identifier
			int(x_atom["residue_serial"]),			# Residue sequence number
			str(x_atom["residue_insertion"]),		# Code for insertion of residues
			float(x_atom["coord_x"]),					# Orthogonal coordinates for X in Angstroms
			float(x_atom["coord_y"]),					# Orthogonal coordinates for Y in Angstroms
			float(x_atom["coord_z"]),					# Orthogonal coordinates for Z in Angstroms
			float(x_atom["occupancy"]),				# Occupancy
			float(x_atom["temperature_factor"]),		# Temperature factor
			str(x_atom["element_symbol"]),			# Element symbol
			str(x_atom["element_charge"])			# Charge on the atom
		)
	# End if

	return s_output_line		# Returns the formatted atom line
# End function ------------------------------------------ #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from convert_grid_pdb import convert_grid_pdb
	#
	# In :
	# Out :

# Usage
# convert_grid_pdb(		#
# = ,		#
# )

# ---------------------------------------------------------------------------- #
