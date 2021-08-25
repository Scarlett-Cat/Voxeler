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

# Program external resources
from resources import elem_config
	# Contains chemical elements custom properties

# Parameters
from config import global_parameters as gp
	# Contains the global variables

# General library
from lib.retrieve_element_mass import retrieve_element_mass
	# Retrieves the mass of each given element
	# In : (a/s) the element symbol, (a/s) element symbol in case of fail
	# Out : (a/f) the element mass
# ---------------------------------------------------------------------------- #



# Class ---------------------------------------------------------------------- #

class PdbStructure:
	"""
	A class representing the fields contained in a PDB file
	"""

	def __init__(self):
		"""
		Initializes the fields of the structure
		"""

		# PDB fields
		self.s_name = ""				# Name of the structure
		self.l_s_leading_data = []		# PDB information written above the atom properties
		self.l_s_trailing_data = []		# PDB information written under the atom properties

		# Structural fields
		self.i_atom_count = 0			# Number of atoms in the structure
		self.a_atoms = None				# Array of atoms properties
		self.a_max_coord = None			# Maximal coordinates for each axis
		self.a_min_coord = None			# Minimal coordinates for each axis

		# Grid fields
		self.a_grid = None				# 3D grid containing the structure
		self.l_l_elements = None		# Set of atoms contained in the structure

		# Solubilization fields
		self.o_tree = None						# A KDTree object representing the exact placement of atoms, used for distance determination

		# Comparison fields
		self.b_loaded = False		# Keeps tracks of the state of the structure

		# Ligands fields
		self.l_o_ligands = []		# A list for the ligands

		# Pocket fields
		self.l_i_pocket_residues = []		# List of the residues included in the pocket
		self.a_pocket_atoms = None			# Array of the pocket atoms properties
		self.a_pocket_grid = None			# 3D grid containing the pocket

		# Miscellaneous fields
		self.a_min_coord = None		# Minimum coordinates of the structure
		self.a_max_coord = None		# Maximum coordinates of the structure
		self.f_mass = 0.0			# Mass of the structure
	# End method


	def __repr__(self):
		"""
		Creates a human friendly representation of the information contained in the object
		"""

		# Preparing variables
		l_s_content = [		# List containing the content to print
			"> The structure object :"
		]

		# PDB fields
		l_s_content.append("s_name : {}".format(self.s_name))

		# Structural fields
		l_s_content.append("i_atom_count : {}".format(self.i_atom_count))
		l_s_content.append("a_atoms : {}".format(len(self.a_atoms)))

		# Grid fields
		l_s_content.append("b_loaded : {}".format(self.b_loaded))
		l_s_content.append("a_grid : {}".format(self.a_grid.size))

		return "\n".join(l_s_content)		# Returns the content to show
	# End method


	def load_structure(self, **kwargs):
		"""
		Loads in the object the base information about the structure
		"""

		# PDB fields
		self.s_name = kwargs["s_name"]								# Name of the structure
		self.l_s_leading_data = kwargs["l_s_leading_data"]			# PDB information written above the atom properties
		self.l_s_trailing_data = kwargs["l_s_trailing_data"]		# PDB information written under the atom properties

		# Structural fields
		self.i_atom_count = len(kwargs["d_atoms"]["element_type"])		# Retrieves the number of atoms
		self.a_atoms = np.arange(self.i_atom_count).astype(				# Array of atoms properties
			np.dtype([
				("element_type", np.str, 6),				# ATOM or HETATM
				("atom_serial", np.uint16, 1),				# Atom serial number
				("atom_name", np.str, 4),					# Atom name
				("alternative_location", np.str, 1),		# Alternate location indicator
				("residue_name", np.str, 3),				# Residue name
				("chain_id", np.str, 1),					# Chain identifier
				("residue_serial", np.int16, 1),			# Residue sequence number
				("residue_insertion", np.str, 1),			# Code for insertion of residues
				("coord_x", np.float32, 1),					# Orthogonal coordinates for X in Angstroms
				("coord_y", np.float32, 1),					# Orthogonal coordinates for Y in Angstroms
				("coord_z", np.float32, 1),					# Orthogonal coordinates for Z in Angstroms
				("occupancy", np.float16, 1),				# Occupancy
				("temperature_factor", np.float16, 1),		# Temperature factor
				("element_symbol", np.str, 2),				# Element symbol
				("element_charge", np.str, 2),				# Charge on the atom
				("element_mass", np.float16, 1),			# Mass of the atom
				("grid_x", np.int16, 1),					# X coordinates in the grid
				("grid_y", np.int16, 1),					# Y coordinates in the grid
				("grid_z", np.int16, 1),					# Z coordinates in the grid
				("custom_type", np.str, 3),					# A custom name for the element
			])
		)

		# For each field to save
		for s_key in kwargs["d_atoms"]:
			self.a_atoms[s_key] = kwargs["d_atoms"][s_key]		# Saves each field of the dictionary of atom properties

		self.a_atoms["element_mass"] = retrieve_element_mass(		# Retrieves the atomic mass of the given elements
			x_element_symbol=self.a_atoms["element_symbol"],		# Element symbol
			x_backup_symbol=self.a_atoms["atom_name"]				# Element symbol in case of fail
		)
		self.translate_custom_types()		# Translates to the custom element types

		self.l_l_elements = set(self.a_atoms["element_symbol"])		# List all the different elements contained in the structure
		l_s_elements = [None] * len(gp.D_ELEMENT_NUMBER)			# Creates an empty list with a slot for each possible element

		# For each chemical element
		for s_element in self.l_l_elements:

			i_element_number = gp.D_ELEMENT_NUMBER[s_element]			# Retrieves the atomic number of the element
			a_element_indexes = np.where(								# Retrieves the indexes of the elements
				self.a_atoms["element_symbol"] == s_element
			)

			l_s_elements[i_element_number] = [		# Orders each element by their atomic number
				s_element,							# Element symbol
				i_element_number,					# Atomic number of the element
				a_element_indexes,					# Indexes of the element in the structure
				None,								# Coordinates of the element in the grid
				None,								# VdW radius of the element
				None								# Sphere coordinates of the element
			]
		# End for

		self.l_l_elements = list(filter(None, l_s_elements))		# Removes empty elements in the list

		# Miscellaneous fields
		self.f_mass = sum(self.a_atoms["element_mass"])		# Sums the mass of each element
	# End method ---------------------------------------- #


	def actualize_properties(self):
		"""
		Actualizes the properties of the structure
		"""

		self.a_max_coord = np.array((			# Computes the maximal coordinates
			max(self.a_atoms["coord_x"]),		# For the x axis
			max(self.a_atoms["coord_y"]),		# For the y axis
			max(self.a_atoms["coord_z"])		# For the z axis
		))
		self.a_min_coord = np.array((			# Computes the minimal coordinates
			min(self.a_atoms["coord_x"]),		# For the x axis
			min(self.a_atoms["coord_y"]),		# For the y axis
			min(self.a_atoms["coord_z"])		# For the z axis
		))
	# End method ---------------------------------------- #


	def delete_grid(self):
		"""
		Deletes the grid and frees memory
		"""

		self.a_grid = None		# Deletes the object from memory
	# End method ---------------------------------------- #


	def translate_custom_types(self):
		"""
		Defines a new element type according to the element role within the structure
		"""

		# Preparing variables
		a_residue_names = self.a_atoms["residue_name"]		# Loads the names of residues
		a_atom_name = self.a_atoms["atom_name"]		# Loads the names of the atoms
		a_atom_symbol = self.a_atoms["element_symbol"]		# Loads the elements symbols
		l_s_custom_types = []		# Contains the list of converted types
		d_translate_custom = {		# Conversion dictionary for custom types
			"O": "OC",
			"H": "H",
			"N": "NAM",
			"C": "XOT",
			"CA": "XOT",
			"CB": "XOT",
			"OXT": "XOT"
		}

		# STEP 1 : Converting the atom types ---------------- #
		# For each element to convert
		for i_element in range(len(a_residue_names)):

			# If the residue is one of the main amino acids
			if a_residue_names[i_element] in elem_config.RES:

				# Hydrogen
				if a_atom_symbol[i_element] == "H":
					s_custom_type = "H"

				# If the atom is one of the main carbon chain
				elif a_atom_name[i_element] in d_translate_custom.keys():
					s_custom_type = d_translate_custom[a_atom_name[i_element]]

				# Nitrogen in Arginine
				elif a_residue_names[i_element] == "ARG" and a_atom_name[i_element] in elem_config.NARG[a_residue_names[i_element]]:
					s_custom_type = "NBAS"

				# Carbon SP2 in aromatic ring
				elif a_residue_names[i_element] in elem_config.CAR.keys() and a_atom_name[i_element] in elem_config.CAR[a_residue_names[i_element]]:
					s_custom_type = "CAR"

				# Oxygen in hydroxyl or phenol
				elif a_residue_names[i_element] in elem_config.OHY.keys() and a_atom_name[i_element] == elem_config.OHY[a_residue_names[i_element]]:
					s_custom_type = "OH"

				# Nitrogen in amide
				elif a_residue_names[i_element] in elem_config.NAM.keys() and a_atom_name[i_element] == elem_config.NAM[a_residue_names[i_element]]:
					s_custom_type = "NAM"

				# Nitrogen in Histidine
				elif a_residue_names[i_element] in elem_config.NHIS.keys() and a_atom_name[i_element] in elem_config.NHIS[a_residue_names[i_element]]:
					s_custom_type = "NBAS"

				# Central carbon from ARG, GLN, GLU, ASP, ASN
				elif a_residue_names[i_element] in elem_config.CE.keys() and elem_config.CE[a_residue_names[i_element]] == a_atom_name[i_element]:
					s_custom_type = "CAR"

				# Oxygen in carbonyl
				elif a_residue_names[i_element] in elem_config.OC.keys() and a_atom_name[i_element] == elem_config.OC[a_residue_names[i_element]]:
					s_custom_type = "OC"

				# Oxygen in carboxylate and oxygen in C-terminal
				elif a_residue_names[i_element] in elem_config.OOX.keys() and \
						(a_atom_name[i_element] == elem_config.OOX[a_residue_names[i_element]][0] or
						 a_atom_name[i_element] == elem_config.OOX[a_residue_names[i_element]][1]):
					s_custom_type = "OOX"

				# Nitrogen in Lysine
				elif a_residue_names[i_element] in elem_config.NLYS.keys() and a_atom_name[i_element] == elem_config.NLYS[a_residue_names[i_element]]:
					s_custom_type = "NBAS"

				# Unknown element within a amino acid
				else:
					s_custom_type = "XOT"
			# End if

			# If the element is a metallic atom
			elif a_atom_symbol[i_element] in elem_config.METAL:
				s_custom_type = "META"

			# If the element is a halogen
			elif a_atom_symbol[i_element] in elem_config.HALO:
				s_custom_type = "HALO"

			# If the element is a water molecule
			elif a_residue_names[i_element] == "HOH" and a_atom_name[i_element] == "O":
				s_custom_type = "OOW"

			# If the element is not known
			else:

				# If the element can be converted
				if a_atom_symbol[i_element] in d_translate_custom.keys():
					s_custom_type = d_translate_custom[a_atom_symbol[i_element]]

				# If it cannot
				else:
					s_custom_type = "HETATM"
			# End if

			l_s_custom_types.append(s_custom_type)		# Saves the new element type
		# End for
		# END STEP 1 ---------------------------------------- #

		# STEP 2 : Saving the list of custom types ---------- #
		self.a_atoms["custom_type"] = l_s_custom_types		# Saves the list of custom types
		# END STEP 2 ---------------------------------------- #
	# End method ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# from pdb_structure import PdbStructure
	# PDB structure and associated grid

# ---------------------------------------------------------------------------- #
