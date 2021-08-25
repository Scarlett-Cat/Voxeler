# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : June 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Classes
from cla.pdb_structure import PdbStructure
	# PDB structure and associated grid

# General library
from lib.read_file_content import read_file_content
	# Extracts the content of a file
	# In : (p) file's path
	# Out : (l(s)) the file's content
from lib.terminate_program_process import terminate_program_process
	# Stops the program and prints content
	# In : (l(s)) content to prompt
	# Out : None

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def parse_pdb_file(p_file, d_filters):
	"""
	Extracts a PDB structure from a file and applies filters
	:param p_file: Path to the PDB file to extract
	:param d_filters: Dictionary of parsing filters to apply
	:return: The PDB structure extracted and saved in an object
	"""

	# STEP 0 : Preparing variables ---------------------- #
	l_s_logs = []					# A list for logs

	l_s_leading_pdb = []			# List of lines before the atoms records
	l_s_trailing_pdb = []			# List of lines after the atom records
	b_atom_encountered = False		# If the first atom line has been encountered
	i_line_count = -1				# Index of the current line of the file

	d_atoms = {}					# A dictionary of valid elements from the parsed PDB
	l_s_keys = [					# A list of PDB fields to save
		"element_type", "atom_serial", "atom_name", "alternative_location", "residue_name",
		"chain_id", "residue_serial", "residue_insertion", "coord_x", "coord_y",
		"coord_z", "occupancy", "temperature_factor", "element_symbol", "element_charge"
	]

	# For each one of the 15 useful fields in the PDB file
	for s_key in l_s_keys:
		d_atoms[s_key] = []		# Creates an empty list for each field
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Extracting the file content -------------- #
	l_s_content = read_file_content(		# Retrieves the content of the file
		p_file=p_file						# Path to the file to extract
	)
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Security check --------------------------- #
	# If the PDB file is empty
	if len(l_s_content) < 2:
		l_s_logs.append("ERROR : The '{}' PDB file is empty".format(p_file))		# Defines the error message
		terminate_program_process(		# Stops the program
			l_s_content=l_s_logs		# Content to save to logs
		)
	# END STEP 2 ---------------------------------------- #


	# STEP 3 : Parsing the PDB file --------------------- #
	# For each line in the PDB file
	for s_line in l_s_content:

		i_line_count += 1		# Counts the lines

		# If the line concerns ATOMs or HETATMs
		if s_line.startswith("ATOM") or s_line.startswith("HETATM"):

			b_atom_encountered = True		# If the first atom line has been encountered

			l_s_atom_properties = [			# Contains the element properties
				s_line[0:6].strip(),		# ATOM or HETATM
				s_line[6:11].strip(),		# Atom serial number
				s_line[12:16].strip(),		# Atom name
				s_line[16].strip(),			# Alternate location indicator
				s_line[17:20].strip(),		# Residue name
				s_line[21].strip(),			# Chain identifier
				s_line[22:26].strip(),		# Residue sequence number
				s_line[26].strip(),			# Code for insertion of residues
				s_line[30:38].strip(),		# Orthogonal coordinates for X in Angstroms
				s_line[38:46].strip(),		# Orthogonal coordinates for Y in Angstroms
				s_line[46:54].strip(),		# Orthogonal coordinates for Z in Angstroms
				s_line[54:60].strip(),		# Occupancy
				s_line[60:66].strip(),		# Temperature factor
				s_line[76:78].strip(),		# Element symbol
				s_line[78:80].strip(),		# Charge on the atom
			]
			b_valid = apply_parsing_filters(		# Applies the parsing filters to the line
				l_s_atom=l_s_atom_properties,		# The atom line to filter
				d_filters=d_filters			# The filters to apply
			)

			# If the line is validated by the filters
			if b_valid:

				# Try to save the line
				try:

					# Tries the conversion of each field
					d_atoms["element_type"].append(l_s_atom_properties[0].strip())						# ATOM or HETATM
					d_atoms["atom_serial"].append(int(l_s_atom_properties[1].strip()))					# Atom serial number
					d_atoms["atom_name"].append(l_s_atom_properties[2].strip())							# Atom name
					d_atoms["alternative_location"].append(l_s_atom_properties[3].strip())				# Alternate location indicator
					d_atoms["residue_name"].append(l_s_atom_properties[4].strip())						# Residue name
					d_atoms["chain_id"].append(l_s_atom_properties[5].strip())							# Chain identifier
					d_atoms["residue_serial"].append(int(l_s_atom_properties[6].strip()))				# Residue sequence number
					d_atoms["residue_insertion"].append(l_s_atom_properties[7].strip())					# Code for insertion of residues
					d_atoms["coord_x"].append(float(l_s_atom_properties[8].strip()))					# Orthogonal coordinates for X in Angstroms
					d_atoms["coord_y"].append(float(l_s_atom_properties[9].strip()))					# Orthogonal coordinates for Y in Angstroms
					d_atoms["coord_z"].append(float(l_s_atom_properties[10].strip()))					# Orthogonal coordinates for Z in Angstroms
					d_atoms["occupancy"].append(float(l_s_atom_properties[11].strip()))					# Occupancy
					d_atoms["temperature_factor"].append(float(l_s_atom_properties[12].strip()))		# Temperature factor
					d_atoms["element_symbol"].append(l_s_atom_properties[13].strip())					# Element symbol
					d_atoms["element_charge"].append(l_s_atom_properties[14].strip())					# Element symbol

				# If there is an error during the conversion
				except OSError:
					l_s_logs.append("ERROR : Incorrect value type in '{}' at line '{}'".format(p_file, i_line_count))		# Defines the error message
					terminate_program_process(		# Stops the program
						l_s_content=l_s_logs		# Content to save to logs
					)
				# End try
			# End if

		# If the line is does not contain atom data
		else:

			# If atomic data have been encountered
			if b_atom_encountered:
				l_s_trailing_pdb.append(s_line)		# Saves the trailing line

			# If atomic data have not been encountered
			else:
				l_s_leading_pdb.append(s_line)		# Saves the leading line
		# End if
	# End for
	# END STEP 3 ---------------------------------------- #


	# STEP 4 : ------------------------------------------ #
	s_name = p_file.split('/')[-1].split('.')[0]		# Extracts the name of the structure
	o_structure = PdbStructure()						# Creates a PDB structure object
	o_structure.load_structure(							# Loads the structure into the object
		s_name=s_name,									# Name of the structure
		l_s_leading_data=l_s_leading_pdb,				# Structure information
		l_s_trailing_data=l_s_trailing_pdb,				# Remarks on the structure
		d_atoms=d_atoms									# Dictionary of atom properties
	)
	# END STEP 4 ---------------------------------------- #


	# STEP 5 : Returns the object ----------------------- #
	return o_structure		# Returns the extracted structure
	# END STEP 5 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Auxiliary functions -------------------------------------------------------- #

def apply_parsing_filters(l_s_atom, d_filters):
	"""
	Applies filters to the PDB parsing, removing unwanted lines
	:param l_s_atom: The atom line to filter
	:param d_filters: Parsing parameters, filters to apply
	:return: A boolean, if the line pass all the filters or not
	"""

	b_valid_line = True					# If the line is valid and can be saved
	l_s_water_id = ["HOH", "OOW"]		# List of possible water molecules identifiers

	# If ATOM are discarded and the line contains a HETATM element
	if d_filters["b_discard_atom"] and l_s_atom[0] == "ATOM":
		return False		# The line cannot be saved

	# If HETATM are discarded and the line contains a HETATM element
	if d_filters["b_discard_hetatm"] and l_s_atom[0] == "HETATM":
		return False		# The line cannot be saved

	# If the Hydrogen must be discarded and the line contains a Hydrogen
	if d_filters["b_discard_hydrogen"] and l_s_atom[13] == " H":
		return False		# The line cannot be saved

	# If water is discarded and the line contains a water molecule
	if d_filters["b_discard_water"] and l_s_atom[4] in l_s_water_id:
		return False		# The line cannot be saved

	# If the alternative location for atom must be discarded and the line contains a alternative location
	if d_filters["b_discard_alternative"] and l_s_atom[3] not in " A":
		return False		# The line cannot be saved

	# If not all chain are considered and the element is from a discarded chain
	if len(d_filters["l_c_chain_white"]) > 0 and l_s_atom[5] not in d_filters["l_c_chain_white"]:
		return False		# The line cannot be saved

	# If the chain id is black listed
	if len(d_filters["l_c_chain_black"]) > 0 and l_s_atom[5] in d_filters["l_c_chain_black"]:
		return False		# The line cannot be saved

	# If not all residues are considered and the element is from a discarded residue
	if len(d_filters["l_s_residue_white"]) > 0 and l_s_atom[4] not in d_filters["l_s_residue_white"]:
		return False		# The line cannot be saved

	# If the residue is black listed
	if len(d_filters["l_s_residue_black"]) > 0 and l_s_atom[4] in d_filters["l_s_residue_black"]:
		return False		# The line cannot be saved

	# If not all residues are considered and the element is from a discarded residue
	if len(d_filters["l_i_residue_white"]) > 0 and l_s_atom[6] not in d_filters["l_i_residue_white"]:
		return False		# The line cannot be saved

	# If the residue id is black listed
	if len(d_filters["l_i_residue_black"]) > 0 and l_s_atom[6] in d_filters["l_i_residue_black"]:
		return False		# The line cannot be saved

	# If not all residues are considered and the element is from a discarded residue
	if len(d_filters["l_s_atom_white"]) > 0 and l_s_atom[13] not in d_filters["l_s_atom_white"]:
		return False		# The line cannot be saved

	# If the residues is black listed
	if len(d_filters["l_s_atom_black"]) > 0 and l_s_atom[13] in d_filters["l_s_atom_black"]:
		return False		# The line cannot be saved

	return True
# End function ------------------------------------------ #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from parse_pdb_file import parse_pdb_file
	# Extracts a PDB structure from a file and applies filters
	# In : (p) PDB file to extract, (d) parsing filters to apply
	# Out : (o) the object containing the structure

# Usage
# parse_pdb_file(			# Extracts a PDB structure into an object
# 	p_file=p_file,			# The PDB file to extract
# 	d_filters=d_filters		# The parsing filters to apply
# )

# ---------------------------------------------------------------------------- #
