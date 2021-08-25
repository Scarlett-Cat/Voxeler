
# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import numpy as np
	# Optimized arrays

# Program external resources
from resources import elements
	# Contains chemical elements properties

# General library
from lib.terminate_program_process import terminate_program_process
	# Stops the program and prints content
	# In : (l(s)) content to prompt
	# Out : None

# ---------------------------------------------------------------------------- #



def retrieve_element_mass(x_element_symbol, x_backup_symbol):
	"""
	Retrieves the mass of each given element
	:param x_element_symbol: Array of element symbol or a single element symbol
	:param x_backup_symbol: Atom name in case of missing element symbol
	:return: The mass corresponding of the given element symbols
	"""

	# STEP 0 : Preparing variables ---------------------- #
	l_s_logs = []		# A list for log
	# END STEP 0 ---------------------------------------- #


	# STEP 1 : Retrieving the element mass -------------- #
	d_element_mass = {key.symbol: elements.ELEMENTS[key.symbol].mass for key in elements.ELEMENTS}
	# END STEP 1 ---------------------------------------- #


	# STEP 2 : Converting the element ------------------- #
	# If the element to convert is a numpy array
	if isinstance(x_element_symbol, np.ndarray):

		a_element_mass = np.zeros(len(x_element_symbol)).astype(np.float32)		# Creates an empty array for the elements mass
		x_element_symbol.tolist()												# Converts the numpy array into a Python list

		# For each element to process
		for i_element in range(len(x_element_symbol)):

			# Tries to retrieve the element mass
			try:
				a_element_mass[i_element] = d_element_mass[x_element_symbol[i_element]]		# Retrieves the mass of the element

			# If there is a symbol error
			except KeyError:
				a_element_mass[i_element] = d_element_mass[x_backup_symbol[i_element]]		# Retrieves the mass of the element

		return a_element_mass		# Returns the array of mass

	# If the element to convert is a string
	elif isinstance(x_element_symbol, str):

		# Tries to retrieve the element mass
		try:
			return d_element_mass[x_element_symbol]		# Retrieves the mass of the element

		# If there is a symbol error
		except KeyError:
			return d_element_mass[x_backup_symbol]  # Retrieves the mass of the element

	# If the argument type is wrong
	else:

		l_s_logs.append("ERROR : Wrong element symbol type '{}'".format(type(type(x_element_symbol))))		# Defines the error message
		terminate_program_process(		# Stops the program
			l_s_content=l_s_logs		# Content to save to logs
		)
	# END STEP 2 ---------------------------------------- #

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from retrieve_element_mass import retrieve_element_mass
	# Retrieves the mass of each given element
	# In : (a/s) the element symbol, (a/s) element symbol in case of fail
	# Out : (a/f) the element mass

# Usage
# retrieve_element_mass(		# Retrieves the atomic mass of the given elements
# 	x_element_symbol=,			# Element symbol
# 	x_backup_symbol=			# Element symbol in case of fail
# )

# ---------------------------------------------------------------------------- #
