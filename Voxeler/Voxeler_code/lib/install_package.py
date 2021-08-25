# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : July 2020
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

# Universal modules
import sys
	# Allows python to access the system commands
import subprocess
	# Allows Python to run external processes

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def install_package(s_package):
	"""
	Installs a package during the Python execution, useful when modules cannot be installed the classic way
	:param s_package: The name of the package to install
	"""

	subprocess.check_call([sys.executable, "-m", "pip", "install", s_package])		# Installs the package

# ---------------------------------------------------------------------------- #



# Reference ------------------------------------------------------------------ #

# Importation
# from install_package import install_package
	# Installs a package during the Python execution, useful when modules cannot be installed the classic way
	# In : (s) the name of the package to install
	# Out : None

# Usage
# install_package(			# Installs a package while Python is running
# s_package=s_package,		# The name of the package to install
# )

# ---------------------------------------------------------------------------- #
