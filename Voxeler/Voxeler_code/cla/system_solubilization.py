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
import sys
	# Allows python to access the system commands

# Program external resources
from resources import elements
	# Contains chemical elements properties

# Parameters
# Classes
# General library
from lib.retrieve_specific_files import retrieve_specific_files
	# Retrieves files path, recursively or not, matching a specific pattern, or not
	# In : (p) directory to retrieve files from, (s) pattern to match,
	# In : (b) if the search needs to be recursive, (i) minimum number of match,
	# In : (i) maximum number of match
	# Out : (l(p)) a list of the file paths founds
from lib.read_file_content import read_file_content
	# Extracts the content of a file
	# In : (p) file's path
	# Out : (l(s)) the file's content
from lib.terminate_program_process import terminate_program_process
	# Stops the program and prints content
	# In : (l(s)) content to prompt
	# Out : None

# Specific modules
from lib.convert_element_symbol import convert_element_symbol
	# Converts atomic number into element symbols or element symbols into atomic numbers
	# In : (a/i/s) the atom symbol to convert
	# Out : (a/i/s) the converted atom data

# ---------------------------------------------------------------------------- #



# Class ---------------------------------------------------------------------- #

class SystemSolubilization:
	"""
	TODO
	"""

	def __init__(self):
		"""
		Initializes the fields
		"""

		# Structure fields
		self.l_o_structures = []				# A list for structure objects
		self.a_offset = np.array((0, 0, 0))		# Coordinates offset compared to the original placement
		self.a_max_coord = None					# Maximal real coordinates for each axis
		self.a_min_coord = None					# Minimal real coordinates for each axis

		# Grid fields
		self.i_grid_bleeding = 0		# Number of extra points around the grids
		self.f_grid_spacing = 5.0		# Distance between two grid points in Angstroms
		self.a_max_grid = None			# Maximal grid coordinates for each axis
		self.a_min_grid = None			# Minimal grid coordinates for each axis
		self.a_grid_size = None			# Size of the grid
		self.i_points_count = 0			# Number of points in the grid

		# Resources fields
		self.d_vdw_radius = {}				# Dictionary of VdW radius
		self.d_scaled_vdw = {}				# VdW radius scaled to the grid
		self.i_max_radius = 0				# The maximal VdW radius scaled to the grid
		self.d_d_distance_score = {}		# Dictionary containing the electronic densities

		# Progression fields
		self.i_current_structure = 0		# Keeps tacks of the structures being processed
		self.i_progress_length = 0			# The length of the progress bar
		self.l_l_tasks = []					# A list containing the tasks to process
		self.i_progress = 0					# The number of completed task
	# End method ---------------------------------------- #

	def __repr__(self):
		"""
		Creates a human friendly representation of the system and it's content
		"""

		# Preparing variables
		l_s_content = [		# List containing the content to print
			"> The solubilization system :"
		]

		# Structure fields
		l_s_content.append("l_o_structures : {}".format(len(self.l_o_structures)))
		l_s_content.append("a_offset : {}".format(self.a_offset))
		l_s_content.append("a_max_coord : {}".format(self.a_max_coord))
		l_s_content.append("a_min_coord : {}".format(self.a_min_coord))

		# Grid fields
		l_s_content.append("i_grid_bleeding : {}".format(self.i_grid_bleeding))
		l_s_content.append("f_grid_spacing : {}".format(self.f_grid_spacing))
		l_s_content.append("a_max_grid : {}".format(self.a_max_grid))
		l_s_content.append("a_min_grid : {}".format(self.a_min_grid))
		l_s_content.append("a_grid_size : {}".format(self.a_grid_size))

		# Resources fields
		l_s_content.append("d_d_distance_score : {} entries".format(len(self.d_d_distance_score.keys())))

		return "\n".join(l_s_content)		# Returns the content to show
	# End method

# System initialization

	def initialize_system(self, d_parameters):
		"""
		Initializes the system base fields
		:param d_parameters: Dictionary of the system parameters
		"""

		# Grid fields
		self.f_grid_spacing = d_parameters["f_grid_spacing"]		# The space between two grid points

		# Resources fields
		self.load_vdw_radius()				# Retrieves the VdW radius of chemical elements
		self.scale_vdw_radius()				# Scales the radius to the grid spacing
		self.load_electronic_densities(		# Retrieves the electronic densities of each atom type
			d_parameters=d_parameters
		)

		# Adapting the grid size
		self.i_grid_bleeding = np.rint(		# Defines the number of points to add around the grids
			self.i_max_radius * 2 + int(d_parameters["l_f_solubilization_radius"][1] / self.f_grid_spacing) * 2
		)

		self.a_offset = self.a_offset + self.i_grid_bleeding		# Records the scaled offset
	# End method ---------------------------------------- #


	def load_vdw_radius(self):
		"""
		Extracts the values of VdW radius
		"""

		# Preparing variables
		d_vdw_radius = {}						# A dictionary containing the VdW radius of each element
		d_o_elements = elements.ELEMENTS		# A dictionary of chemical elements properties

		# For each known element within the list of chemical elements
		for o_element in d_o_elements:
			d_vdw_radius[o_element.symbol] = d_o_elements[o_element.symbol].vdwrad		# Saves the VdW radius of each chemical element

		self.d_vdw_radius = d_vdw_radius		# Saves the dictionary of VdW radius
	# End method ---------------------------------------- #


	def scale_vdw_radius(self):
		"""
		Scales the VdW radius to the grid spacing, converting the radius from Angstroms to points
		"""

		# Preparing variables
		i_max_radius = 0		# Stores the maximal radius encountered
		d_scaled_vdw = {}		# Dictionary of VdW radius per element scaled with the grid spacing

		# For each element in the dictionary
		for s_element in self.d_vdw_radius:
			d_scaled_vdw[s_element] = np.rint(		# Defines the VdW radius in points
				self.d_vdw_radius[s_element] / self.f_grid_spacing #WHY ?
			).astype(np.int32)

			# If the radius is superior to the actual maximal
			if d_scaled_vdw[s_element] > i_max_radius:
				i_max_radius = d_scaled_vdw[s_element]		# Saves the actual radius
		# End for

		self.i_max_radius = i_max_radius		# Saves the maximal radius
		self.d_scaled_vdw = d_scaled_vdw		# Saves the scaled radius
	# End method ---------------------------------------- #


	def load_electronic_densities(self, d_parameters):
		"""
		Retrieves the electronic densities of each atom type
		:param d_parameters: Dictionary of parameters containing the path to the density files
		"""

		# Preparing variables
		d_d_distance_score = {}		# Dictionary of scoring dictionaries

		# Retrieving the density files
		l_p_density_files = retrieve_specific_files(				# Retrieves the path of files within a directory
			p_directory=d_parameters["p_electronic_densities"],		# Path to the directory to explore
			s_pattern="Oow*",										# Pattern to match within the directories
			b_recursive=True,										# If the function will search in subdirectories
			i_min_match=1,											# Minimal number of expected matches
			i_max_match=999999										# Maximal number of expected matches
		)

		# For each density file
		for p_file in l_p_density_files:

			d_score_rule_buffer = {}								# A temporary dictionary buffer
			a_distance = np.zeros(512).astype(np.float32)			# Array of electronic distances
			a_density = np.zeros(512).astype(np.float16)			# Array of electronic densities
			l_s_content = read_file_content(						# Retrieves the content of the file
				p_file=p_file,										# Path to the file to read
			)[1:]													# Do not take the first line
			s_file_name = p_file.split('/')[-1].split('.')[0]		# Extracts the last string before the extension

			# For each line in the density file
			for i_line in range(len(l_s_content)):
				a_distance[i_line] = l_s_content[i_line].split(' ')[1]		# Retrieves the distance
				a_density[i_line] = l_s_content[i_line].split(' ')[2]		# Retrieves the density

			# If the density needs to be normalized
			if d_parameters["b_normalize_densities"]:
				a_density = a_density / max(a_density)		# Normalizes the density value

			# For each distances in the file
			for i_distance in range(len(a_distance)):
				d_score_rule_buffer[a_distance[i_distance]] = a_density[i_distance]		# Saves the score corresponding to the distance

			d_d_distance_score[s_file_name.upper()] = d_score_rule_buffer		# Saves the temporary dict into the final one
			self.d_d_distance_score = d_d_distance_score						# Saves the electronic densities in the system object
		# End for
	# End method ---------------------------------------- #


	def actualize_properties(self):
		"""
		Actualizes the properties of the system
		"""

		# Preparing variables
		self.a_max_coord = np.array((np.NaN, np.NaN, np.NaN))		# Maximal coordinates of the structure
		self.a_min_coord = np.array((np.NaN, np.NaN, np.NaN))		# Minimal coordinates of the structure

		# For each registered structure
		for o_structure in self.l_o_structures:

			o_structure.actualize_properties()		# Actualizes the properties of the structure, such as coordinates
			self.a_max_coord = np.array((									# Computes the maximal coordinates
				max(o_structure.a_max_coord[0], self.a_max_coord[0]),		# For the x axis
				max(o_structure.a_max_coord[1], self.a_max_coord[1]),		# For the y axis
				max(o_structure.a_max_coord[2], self.a_max_coord[2])		# For the z axis
			))
			self.a_min_coord = np.array((									# Computes the minimal coordinates
				min(o_structure.a_min_coord[0], self.a_min_coord[0]),		# For the x axis
				min(o_structure.a_min_coord[1], self.a_min_coord[1]),		# For the y axis
				min(o_structure.a_min_coord[2], self.a_min_coord[2])		# For the z axis
			))
		# End for

		# Actualizing the grid properties
		self.a_max_grid = np.rint(self.a_max_coord / self.f_grid_spacing)		# Maximal coordinates in the grid
		self.a_min_grid = np.rint(self.a_min_coord / self.f_grid_spacing)		# Minimal coordinates in the grid

		# For each axis x, y and z
		for i_axis in range(3):
			self.a_offset[i_axis] -= self.a_min_grid[i_axis]		# Centers all the structures

		self.a_grid_size = np.array((		# Computes the number of points in each dimension of the grid
			(self.i_grid_bleeding * 2 + self.a_max_grid[0] - self.a_min_grid[0]),
			(self.i_grid_bleeding * 2 + self.a_max_grid[1] - self.a_min_grid[1]),
			(self.i_grid_bleeding * 2 + self.a_max_grid[2] - self.a_min_grid[2])
		)).astype(np.int32)
		self.i_points_count = (self.a_grid_size[0] * self.a_grid_size[1] * self.a_grid_size[2])		# Number of points in the grid
	# End method ---------------------------------------- #

# Grid management

	def generate_grid(self, o_structure, d_parameters):
		"""
		Generates a grid for a given PDB structure
		:param o_structure: The structure to be loaded into a grid
		:param d_parameters: Dictionary of the program parameters
		"""

		# Preparing variables
		l_l_elements = o_structure.l_l_elements		# Shortcut for the structure field containing sorted data for each element

		# Converting the structure coordinates
		o_structure.a_atoms["grid_x"] = np.rint(o_structure.a_atoms["coord_x"] / self.f_grid_spacing + self.a_offset[0])		# Converts the real X coordinates for the grid
		o_structure.a_atoms["grid_y"] = np.rint(o_structure.a_atoms["coord_y"] / self.f_grid_spacing + self.a_offset[1])		# Converts the real Y coordinates for the grid
		o_structure.a_atoms["grid_z"] = np.rint(o_structure.a_atoms["coord_z"] / self.f_grid_spacing + self.a_offset[2])		# Converts the real Z coordinates for the grid

		# For each chemical element in the structure
		for i_element in range(len(l_l_elements)):

			a_element_indexes = l_l_elements[i_element][2]				# Loads the indexes of the element
			l_l_elements[i_element][3] = (								# Retrieves the coordinates of the element
				o_structure.a_atoms["grid_x"][a_element_indexes],
				o_structure.a_atoms["grid_y"][a_element_indexes],
				o_structure.a_atoms["grid_z"][a_element_indexes]
			)
			o_structure.l_l_elements[i_element][3] = np.transpose(l_l_elements[i_element][3])		# Formats and saves the atom coordinates
		# End for

		# Creating the grid
		o_structure.a_grid = np.zeros(self.i_points_count).reshape(		# Initializes the grid
			self.a_grid_size[0],
			self.a_grid_size[1],
			self.a_grid_size[2]
		).astype(
			np.dtype([									# Defines the content of each point
				("element_symbol", np.uint8, 1),		# The atom symbol
				("atom_serial", np.uint16, 1),			# The atom serial number
				("score", np.float16, 1)				# The atom score
			])
		)

		# Loading the structure into the grid
		a_elements_code = convert_element_symbol(					# Switches between atom symbols and atomic numbers
			x_element=o_structure.a_atoms["element_symbol"],		# The elements to convert
		)
		o_structure.a_grid["element_symbol"][		# Loads the elements symbols in the grid
			o_structure.a_atoms["grid_x"],
			o_structure.a_atoms["grid_y"],
			o_structure.a_atoms["grid_z"]
		] = a_elements_code
		o_structure.a_grid["atom_serial"][		# Loads the elements symbols in the grid
			o_structure.a_atoms["grid_x"],
			o_structure.a_atoms["grid_y"],
			o_structure.a_atoms["grid_z"]
		] = o_structure.a_atoms["atom_serial"]

		# Generating VDW volumes
		self.generate_vdw_spheres(			# Generates the VdW volumes around each atom
			o_structure=o_structure,		# The structure to be incorporated
			d_parameters=d_parameters		# Dictionary of the program parameters
		)
	# End method ---------------------------------------- #


	def generate_vdw_spheres(self, o_structure, d_parameters, i_solubilization_radius=0, b_remove_volume=False):
		"""
		Generates the VdW volumes around each atom in a specific grid
		:param o_structure: The structure to be incorporated with its VdW volume
		:param d_parameters: Dictionary of the program parameters
		:param i_solubilization_radius: Radius used for the solubilization
		:param b_remove_volume: If the selected volume needs to be deleted
		"""

		# Preparing variables
		l_l_elements = o_structure.l_l_elements		# Shortcut for the structure field containing sorted data for each element

		# For each atom type in the structure
		for i_element in range(len(l_l_elements)):

			i_radius = self.d_scaled_vdw[l_l_elements[i_element][0]] + i_solubilization_radius		# Retrieves the VdW radius of the element
			l_i_radius_range = list(range(-i_radius, i_radius + 1))									# Builds a list of distances included in the sphere
			l_l_elements[i_element][4] = i_radius													# Saves the VdW radius
			l_l_elements[i_element][5] = self.create_vdw_sphere(
				d_parameters=d_parameters,				# Dictionary of the program parameters
				i_radius=i_radius,						# VdW radius of the element
				l_i_radius_range=l_i_radius_range		# Range of radius around the element
			)

		# For each chemical element present
		for i_element in range(len(l_l_elements)):

			# For each atom in the structure
			for a_atom in l_l_elements[i_element][3]:

				i_symbol = o_structure.a_grid["element_symbol"][a_atom[0], a_atom[1], a_atom[2]]
				i_serial = o_structure.a_grid["atom_serial"][a_atom[0], a_atom[1], a_atom[2]]

				a_sphere_coords = (									# Retrieves the coordinates of each point of the sphere
					l_l_elements[i_element][5][0] + a_atom[0],		# X coordinates
					l_l_elements[i_element][5][1] + a_atom[1],		# Y coordinates
					l_l_elements[i_element][5][2] + a_atom[2]		# Z coordinates
				)

				# If the volume needs to be removed
				if b_remove_volume:
					o_structure.a_grid["element_symbol"][a_sphere_coords] = 0		# Deletes the element symbol
					o_structure.a_grid["atom_serial"][a_sphere_coords] = 0			# Deletes the atom serial number

				# If the volume need to be filled with the atom properties
				else:
					o_structure.a_grid["element_symbol"][a_sphere_coords] = i_symbol		# Fills the sphere with the element
					o_structure.a_grid["atom_serial"][a_sphere_coords] = i_serial			# Fills the sphere with the atom serial number
		# End for
	# End method ---------------------------------------- #


	@staticmethod
	def create_vdw_sphere(d_parameters, i_radius, l_i_radius_range):
		"""
		Defines a spherical array of points for a specific geometry
		:param d_parameters: Dictionary of the program parameters
		:param i_radius: The raw VdW radius of the element, in grid points
		:param l_i_radius_range: The list of points in the radius range
		:return: A spherical array of relative coordinates
		"""

		# Preparing variables
		s_grid_geometry = d_parameters["s_grid_geometry"].upper()		# Converts to uppercase the grid geometry
		a_sphere = None													# Array of points coordinates for each grid point in the VdW range
		l_s_logs = []													# Creates an empty list for logs

		# If the grid geometry is a taxicab
		if s_grid_geometry == "TAXICAB":
			a_sphere = np.array(		# The taxicab sphere formula
				[
					(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range
					if abs(x) + abs(y) + abs(z) <= i_radius
				]
			).astype(np.int32)

		# If the grid geometry is uniform
		elif s_grid_geometry == "UNIFORM":
			a_sphere = np.array(		# The uniform sphere formula
				[
					(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range
				]
			).astype(np.int32)

		# If the grid geometry is a classic sphere
		elif s_grid_geometry == "SPHERE":
			a_sphere = np.array(		# The classic sphere formula
				[
					(x, y, z) for x in l_i_radius_range for y in l_i_radius_range for z in l_i_radius_range
					if x ** 2 + y ** 2 + z ** 2 <= i_radius ** 2
				]
			).astype(np.int32)

		# If the geometry is unknown
		else:
			l_s_logs.append("ERROR : Unknown sphere geometry '{}', known geometries are 'taxicab', 'uniform' and 'sphere'.".format(s_grid_geometry))		# Defines the error message
			terminate_program_process(		# Stops the program
				l_s_content=l_s_logs		# Content to save in the logs
			)
		# End if
		return a_sphere.T		# Returns the spherical array of points
	# End method ---------------------------------------- #


	def retrieve_nearest_score(self, s_interaction, f_distance):
		"""
		Retrieves the score corresponding to the closest distance to the query
		:param s_interaction: The atom_atom_neighbourindex interaction to analyse
		:param f_distance: The distance, in Angstroms, to compare
		:return: The score corresponding to the distance closest to the query
		"""

		# If the interaction is present within the density files
		if s_interaction in self.d_d_distance_score:
			l_f_distances = list(self.d_d_distance_score[s_interaction].keys())														# Loads every distances recorded for this interaction
			i_nearest = l_f_distances[min(range(len(l_f_distances)), key=lambda i: abs(l_f_distances[i] - f_distance).any())]		# Retrieves the index of the nearest recorded distance

			return self.d_d_distance_score[s_interaction][i_nearest]		# Returns the score corresponding to the distance closest to the query

		# If this interaction does not exist
		else:
			return 0.0		# Returns a null score
	# End method ---------------------------------------- #

# Progression

	def setup_progress(self, s_text=""):
		"""
		Initializes a progression bar
		"""

		# Preparing the toolbar

		# If the text field is not empty
		if s_text != "":
			s_header = "|{:^46}|".format(s_text)

		# If the text field is empty
		else:
			s_header = "|{:^46}|".format("Running structure comparison...")

		print(s_header)
		self.i_progress_length = 50		# Defines the length of the progression bar
		sys.stdout.write(				# Allocates space to the progression bar
			"|%s|" % (' ' * self.i_progress_length)
		)
		sys.stdout.flush()				# Do not resets the terminal
		sys.stdout.write(				# Writes the first char of the bar
			"\b" * (self.i_progress_length + 1)
		)
	# End method ---------------------------------------- #

	def update_progress(self):
		"""
		Updates the progression bar
		"""

		i_previous_percentage = int(		# Determines the "percentage" before the update
			self.i_progress * 50 / len(self.l_l_tasks)
		)
		self.i_progress += 1				# Updates the number of completed tasks
		i_current_percentage = int(			# Determines the percentage after the update
			self.i_progress * 50 / len(self.l_l_tasks)
		)

		# For each missing percentage since the last update
		for i in range((i_current_percentage - i_previous_percentage)): #TODO
			sys.stdout.write("-")		# Displays the progression
			sys.stdout.flush()			# Without resetting the terminal
	# End method ---------------------------------------- #

	def close_progress(self):
		"""
		Closes the progression bar
		"""
		sys.stdout.write("|\n")		# Ending the progression bar
		self.i_progress = 0			# Resets the progress advancement
	# End method ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# from system_solubilization import SystemSolubilization
	#

# ---------------------------------------------------------------------------- #
