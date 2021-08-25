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

# Parameters
from config import global_parameters as gp
	# Contains the global variables

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def convert_element_symbol(x_element):
	"""
	Converts atomic number into element symbols or element symbols into atomic numbers
	:param x_element: An array of element symbol or atomic numbers, or a single element symbol or atomic number
	:return: The corresponding atomic number or element symbol for each atom
	"""

	# If the variable is an array
	if isinstance(x_element, np.ndarray):
		x_element = x_element.tolist()  # Converts the array into a list

		# If the array contains strings
		if isinstance(x_element[0], str):
			x_converted = np.zeros(len(x_element)).astype(np.uint8)  # Creates an empty array of ints

			# For each strings in the array
			for i_symbol in range(len(x_element)):
				x_converted[i_symbol] = gp.D_ELEMENT_NUMBER[x_element[i_symbol]]  # Adds the atomic number

		# If the array contains ints
		elif isinstance(x_element[0], int):
			x_converted = np.zeros(len(x_element)).astype(np.str)  # Creates an empty array of strings

			# For each strings in the array
			for i_number in range(len(x_element)):
				x_converted[i_number] = gp.D_NUMBER_ELEMENT[x_element[i_number]]  # Adds the element symbol

	# If the variable is a string
	elif isinstance(x_element, str):
		x_converted = gp.D_ELEMENT_NUMBER[x_element]  # Gets the atom number of the element

	# If the variable is an int
	elif isinstance(x_element, int):
		x_converted = gp.D_NUMBER_ELEMENT[x_element]  # Gets the element symbol of the atom
	# End if

	return x_converted  # Returns the atomic number or element symbol

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from convert_element_symbol import convert_element_symbol
	# Converts atomic number into element symbols or element symbols into atomic numbers
	# In : (a/i/s) the atom symbol to convert
	# Out : (a/i/s) the converted atom data

# Usage
# convert_element_symbol(		# Switches between atom symbols and atomic numbers
# 	x_element=a_elements,		# The elements to convert
# )

# ---------------------------------------------------------------------------- #
