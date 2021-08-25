# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import datetime
	# Time manipulation
from ete3 import Tree, TreeStyle
	# Phylogenetic tree tool
from copy import deepcopy
	# Allows the total copy of variables

# General library
from lib.write_file_content import write_file_content
	# Writes content to a file
	# In : (p) output file's path, (c) writing mode,
	# In : (l(s)) content to write
	# Out : None

# ---------------------------------------------------------------------------- #



# Class ---------------------------------------------------------------------- #

class TreePlot:
	"""
	A class converting a matrix of scores into a phylogenetic tree
	"""

	def __init__(self):
		"""
		Initializes the tree fields
		"""

		# Raw data fields
		self.d_nodes = {}			# The list of nodes registered
		self.l_s_nodes = []			# The list of nodes names
		self.l_l_matrix = []		# The matrix of scores
		self.i_node_index = 0		# The index of the last node registered

		# Computed fields
		self.s_tree = ""		# The tree as text, in newick format
		self.o_tree = None		# The tree object generated from the newick tree
		self.o_style = None		# The tree style to apply when rendered
	# End method ---------------------------------------- #

	def __repr__(self):
		"""
		Displays the fields of the tree
		"""

		l_s_content = [		# List containing the content to print
			"> The tree's matrix :"
		]
		s_line = " " * 8		# The first line indentation

		# For each node in the matrix
		for s_node in self.l_s_nodes:
			s_line += "{:<8}".format(s_node)		# Appends the name of each node

		l_s_content.append(s_line)		# Saves the first line

		# For each node in the matrix
		for i_node_index in range(len(self.l_l_matrix)):
			s_line = "{:<8}".format(self.l_s_nodes[i_node_index])		# Starts the line with the name of the node

			# For each node within
			for f_node in self.l_l_matrix[i_node_index]:

				# If the node value is None
				if f_node is None:
					s_line += "{:<8}".format(str(f_node))		# Appends None

				# If the node value is valid
				else:
					s_line += "{:<8.2f}".format(f_node)		# Appends the score of the two concerned nodes
			# End for

			l_s_content.append(s_line)		# Appends each line
		# End for

		return "\n".join(l_s_content)		# Returns the content to show
	# End method ---------------------------------------- #

# Data gathering

	def add_score(self, s_first_node, s_second_node, f_score):
		"""
		Adds a score and it's nodes, if necessary, to the matrix of scores
		:param s_first_node: The first element
		:param s_second_node: The second element
		:param f_score: The score shared by the elements
		"""

		# For each node to actualize
		for s_node in [s_first_node, s_second_node]:

			# If the node is not registered in the dictionary
			if s_node not in self.d_nodes.keys():

				self.d_nodes[s_node] = self.i_node_index				# Registers the node in the dictionary
				self.l_s_nodes.append(s_node)							# Saves the node name
				self.l_l_matrix.append([0.0] * self.i_node_index)		# Creates a place for each missing score
				self.i_node_index += 1									# Adds one to the number of registered nodes

				# For each node in the matrix
				for l_node in self.l_l_matrix:

					# If the list of scores is too short
					if len(l_node) < self.i_node_index:
						l_node.append(0.0)		# Appends an empty item to the list
				# End for
			# End if
		# End for

		i_first_index = self.d_nodes[s_first_node]						# Loads the index of the first node
		i_second_index = self.d_nodes[s_second_node]					# Loads the index of the second node
		self.l_l_matrix[i_first_index][i_second_index] = f_score		# Saves the score in the symmetrical matrix
		self.l_l_matrix[i_second_index][i_first_index] = f_score		# Saves the score in the symmetrical matrix
	# End method ---------------------------------------- #

# Tree generation

	def generate_tree(self, d_parameters):
		"""
		Generates a tree from the matrix of scores
		:param d_parameters: Parameters used for the tree generation and display
		"""

		# Computing the tree topology
		self.compute_tree(d_parameters=d_parameters)		# Computes the topology of the tree
		self.format_tree(d_parameters=d_parameters)			# Creates a tree object and format its style

	# End method ---------------------------------------- #

	def compute_tree(self, d_parameters):
		"""
		Computes the tree topology and saves it into a newick format
		:param d_parameters: Parameters used for the tree computation
		"""

		# Preparing variables
		l_l_matrix = deepcopy(self.l_l_matrix)		# Copies the matrix of scores
		l_s_nodes = deepcopy(self.l_s_nodes)		# Copies the list of leaf names
		f_max = 0.0									# Maximal score buffer
		i_max_first = 0								# X coordinates of the maximal score
		i_max_second = 0							# Y coordinates of the maximal score

		# Sorting the matrix
		# For each node merging needed
		for i_loop in range(len(self.l_s_nodes) - 1):

			# For each node in the matrix
			for i_first_node in range(len(l_l_matrix)):

				# For each node in the matrix
				for i_second_node in range(len(l_l_matrix)):

					# If the score is actually the best
					if l_l_matrix[i_first_node][i_second_node] > f_max:

						f_max = l_l_matrix[i_first_node][i_second_node]		# Saves the maximal score
						i_max_first = i_first_node							# Saves the x coordinates of the maximal score
						i_max_second = i_second_node						# Saves the y coordinates of the maximal score
				# End for
			# End for

			# If the tree only needs topology
			if d_parameters["b_only_topology"]:
				f_max = 1									# Defines the branch length as 1
				d_parameters["b_branch_length"] = False		# Do not display the branch lengths

			l_s_nodes[i_max_first] = "({}:{}, {}:{})".format(		# Formats the tree as text
				l_s_nodes[i_max_first],
				f_max,
				l_s_nodes[i_max_second],
				f_max,
			)

			# For each element, computes the mean scores
			for i_node in range(len(l_l_matrix)):
				l_l_matrix[i_max_first][i_node] = (l_l_matrix[i_max_first][i_node] + l_l_matrix[i_max_second][i_node]) / 2		# New average score
				l_l_matrix[i_node][i_max_first] = (l_l_matrix[i_node][i_max_first] + l_l_matrix[i_node][i_max_second]) / 2		# New average score

			l_l_matrix[i_max_first][i_max_first] = 0

			# For each element, deletes the merged elements
			for i_node in range(len(l_l_matrix)):
				l_l_matrix[i_node].pop(i_max_second)		# Deletes the row of the merged elements

			l_s_nodes.pop(i_max_second)			# Deletes the name of the merged element
			l_l_matrix.pop(i_max_second)		# Deletes the column of the merged element

			# For each element, resets the score of an element compared with itself
			for i_node in range(len(l_l_matrix)):
				l_l_matrix[i_node][i_node] = -1.0		# Sets to -1 the score of an element with itself

			f_max = 0.0				# Resets the maximum score comparator
			i_max_first = 0			# Resets the x index of the maximum score
			i_max_second = 0		# Resets the y index of the maximum score
		# End for

		self.s_tree = l_s_nodes[0] + ';'		# Saves the tree
	# End method ---------------------------------------- #

	def format_tree(self, d_parameters):
		"""
		Formats the newick tree as a tree object and applies style
		:param d_parameters: The parameters used for the tree format
		"""

		# Generating the tree
		self.o_tree = Tree(self.s_tree)		# Converts the newick tree into a tree object
		o_style = TreeStyle()			# Creates the tree style object

		# Defining the tree style
		# If the tree needs to be circular
		if d_parameters["s_tree_shape"].upper() == "CIRCULAR":
			o_style.mode = 'c'		# Defines the tree shape as fully circular

		# If the tree needs to be semicircular
		elif d_parameters["s_tree_shape"].upper() == "SEMICIRCULAR":
			o_style.mode = 'c'				# Defines the tree shape as fully circular
			o_style.arc_start = -180		# The starting point of the semicircle
			o_style.arc_span = 180			# The ending point of the semi circle

		# If the tree needs to be classic
		elif d_parameters["s_tree_shape"].upper() != "CLASSIC":

			l_s_content = [		# Defines the error message
				"ERROR : Unknown tree shape '{}', shape defined as classic".format(d_parameters["s_tree_shape"])
			]
			write_file_content(				# Writes the error to the logs
				p_file=None,				# Path to the file to be written
				s_writing_mode='a',			# Writing mode
				l_s_content=l_s_content		# Content to write
			)
		# End if

		o_style.show_leaf_name = d_parameters["b_leaf_name"]				# If the leaf names needs to be displayed
		o_style.show_branch_length = d_parameters["b_branch_length"]		# If the branch length needs to be displayed

		self.o_style = o_style		# Saves the tree style
	# End method ---------------------------------------- #

	def show_tree(self):
		"""
		Displays the tree in a OS window
		"""

		self.o_tree.show(
			tree_style=self.o_style
		)
	# End method ---------------------------------------- #

	def save_tree(self, d_parameters):
		"""
		Saves the tree in a file alongside its parameters
		:param d_parameters: The parameters used for saving the tree
		"""

		# Formatting the output files
		# If a file name has been given by the user
		if d_parameters["s_tree_name"] is not None:
			s_tree_name = d_parameters["s_tree_name"]		# Uses the user tree name

		# If the program needs to generate the tree name
		else:
			t_time = datetime.datetime.now()		# Retrieves the current time
			s_tree_name = "tree__{}_structures__{}_Ang__{}h{}m".format(
				self.i_node_index,
				str(d_parameters["f_grid_spacing"]).replace('.', '_'),
				t_time.hour,
				t_time.minute
			)

		# Creating the description file content
		l_s_content = [
			"Parameters of {}\n".format(s_tree_name)
		]

		# For each parameter used
		for s_key in d_parameters:

			# If the parameters are not files path
			if s_key[0:2] != "p_":

				l_s_content.append(		# Saves each parameters in the file
					"'{}' = '{}'".format(
						"_".join(s_key.split('_')[1:]),		# Removing variable prefixes
						d_parameters[s_key]
					)
				)

		# Creating the output paths
		p_tree = d_parameters["p_output_comparison"] + '/' + s_tree_name + ".png"
		p_text = d_parameters["p_output_comparison"] + '/' + s_tree_name + ".txt"

		# Saving the files
		i_width = 300 + 10 * self.i_node_index		# Computes the size of the output file
		self.o_tree.render(				# Saves the tree
			p_tree,						# The name of the output file
			w=i_width,					# The width of the rendered image
			units="px",					# The width unit
			tree_style=self.o_style		# The tree style to apply
		)
		write_file_content(				# Writes the error to the logs
			p_file=p_text,				# Path to the file to be written
			s_writing_mode='w',			# Writing mode
			l_s_content=l_s_content		# Content to write
		)
	# End method ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# from tree_plot import TreePlot
	# Generates, renders and saves trees

# ---------------------------------------------------------------------------- #
