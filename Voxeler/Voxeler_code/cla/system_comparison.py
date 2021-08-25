# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : June 2020
# Updated : July 2020
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import sys
	# Allows python to access the system commands
import numpy as np
	# Allows Numpy array manipulation

# Program external resources
from resources import elements
	# Contains chemical elements properties

# General library
from lib.convert_element_symbol import convert_element_symbol
	# Converts atomic number into element symbols or element symbols into atomic numbers
	# In : (a/i/s) the atom symbol to convert
	# Out : (a/i/s) the converted atom data
from lib.terminate_program_process import terminate_program_process
	# Stops the program and prints content
	# In : (l(s)) content to prompt
	# Out : None

# ---------------------------------------------------------------------------- #



# Class ---------------------------------------------------------------------- #

class SystemComparison:
	"""
	A class containing everything needed during the structure comparison
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
		self.f_grid_spacing = 0.0		# Distance between two grid points in Angstroms
		self.a_max_grid = None			# Maximal grid coordinates for each axis
		self.a_min_grid = None			# Minimal grid coordinates for each axis
		self.a_grid_size = None			# Size of the grid
		self.i_points_count = 0			# Number of points in the grid

		# Resources fields
		self.d_vdw_radius = {}		# Dictionary of VdW radius
		self.d_scaled_vdw = {}		# VdW radius scaled to the grid
		self.i_max_radius = 0		# The maximal VdW radius scaled to the grid

		# Comparison fields
		self.l_l_tasks = []				# List of tasks to perform
		self.i_progress_length = 0		# The length of the progress bar
		self.i_progress = 0				# The number of completed task
	# End method


	def __repr__(self):
		"""
		Creates a human friendly representation of the system and it's content
		"""

		# Preparing variables
		l_s_content = [		# List containing the content to print
			"> The comparison system :"
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
		self.load_vdw_radius()		# Retrieves the VdW radius of chemical elements
		self.scale_vdw_radius()		# Scales the radius to the grid spacing

		# Adapting the grid size
		self.i_grid_bleeding = np.rint(		# Defines the number of points to add around the grids
			self.i_max_radius + 1
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
				self.d_vdw_radius[s_element] / self.f_grid_spacing
			).astype(np.int32)

			# If the radius is superior to the actual maximal
			if d_scaled_vdw[s_element] > i_max_radius:
				i_max_radius = d_scaled_vdw[s_element]		# Saves the actual radius
		# End for

		self.i_max_radius = i_max_radius		# Saves the maximal radius
		self.d_scaled_vdw = d_scaled_vdw		# Saves the scaled radius
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

		# If the structure has never been loaded into a grid before
		if o_structure.b_loaded is not True:

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
		# End if

		# Creating the grid
		o_structure.a_grid = np.zeros(self.i_points_count).reshape(		# Initializes the grid
			self.a_grid_size[0],
			self.a_grid_size[1],
			self.a_grid_size[2]
		).astype(np.uint8)		# Stores values between 0 and 127

		# Loading the structure into the grid
		a_elements_code = convert_element_symbol(					# Switches between atom symbols and atomic numbers
			x_element=o_structure.a_atoms["element_symbol"],		# The elements to convert
		)
		o_structure.a_grid[
			o_structure.a_atoms["grid_x"],
			o_structure.a_atoms["grid_y"],
			o_structure.a_atoms["grid_z"]
		] = a_elements_code

		# Generating VDW volumes
		self.generate_vdw_spheres(			# Generates the VdW volumes around each atom
			o_structure=o_structure,		# The structure to be incorporated
			d_parameters=d_parameters		# Dictionary of the program parameters
		)
	# End method ---------------------------------------- #


	def generate_vdw_spheres(self, o_structure, d_parameters):
		"""
		Generates the VdW volumes around each atom in a specific grid
		:param o_structure: The structure to be incorporated with its VdW volume
		:param d_parameters: Dictionary of the program parameters
		"""

		# Preparing variables
		l_l_elements = o_structure.l_l_elements		# Shortcut for the structure field containing sorted data for each element

		# If the VdW radius by element has not been already retrieved
		if o_structure.b_loaded is not True:

			# For each atom type in the structure
			for i_element in range(len(l_l_elements)):

				i_radius = self.d_scaled_vdw[l_l_elements[i_element][0]]		# Retrieves the VdW radius of the element
				l_i_radius_range = list(range(-i_radius, i_radius + 1))			# Builds a list of distances included in the sphere
				l_l_elements[i_element][4] = i_radius							# Saves the VdW radius
				l_l_elements[i_element][5] = self.create_vdw_sphere(
					d_parameters=d_parameters,				# Dictionary of the program parameters
					i_radius=i_radius,						# VdW radius of the element
					l_i_radius_range=l_i_radius_range		# Range of radius around the element
				)
				o_structure.b_loaded = True		# Sets the structure as loaded

		# For each chemical element present
		for i_element in range(len(l_l_elements)):

			# For each atom in the structure
			for a_atom in l_l_elements[i_element][3]:

				a_sphere_coords = (									# Retrieves the coordinates of each point of the sphere
					l_l_elements[i_element][5][0] + a_atom[0],		# X coordinates
					l_l_elements[i_element][5][1] + a_atom[1],		# Y coordinates
					l_l_elements[i_element][5][2] + a_atom[2]		# Z coordinates
				)

				# If the comparison uses the type of elements
				if d_parameters["b_consider_elements"]:
					o_structure.a_grid[a_sphere_coords] = l_l_elements[i_element][1]		# Fills the sphere with the element

				# If only the volume is considered
				else:
					o_structure.a_grid[a_sphere_coords] = 1		# Fills the sphere with the same element
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
	# End method

# Progression

	def setup_progress(self):
		"""
		Initializes a progression bar
		"""

		# Preparing the toolbar
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

	# End method

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
		for i in range((i_current_percentage - i_previous_percentage)):
			sys.stdout.write("-")		# Displays the progression
			sys.stdout.flush()			# Without resetting the terminal
	# End method

	def close_progress(self):
		"""
		Closes the progression bar
		"""
		sys.stdout.write("|\n")		# Ending the progression bar
		self.i_progress = 0			# Resets the progress advancement
	# End method

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# from system_comparison import SystemComparison
	# Singleton object containing every structure to be compared

# ---------------------------------------------------------------------------- #
