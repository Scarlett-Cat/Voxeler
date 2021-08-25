# Information ---------------------------------------------------------------- #
# Author :	Anna Spampinato
# Github :	Scarlett-Cat
# Created : February 2021
# Updated :
# ---------------------------------------------------------------------------- #


# Importations --------------------------------------------------------------- #
import numpy as np

# Classes -------------------------------------------------------------------- #

# General library
from config import global_parameters as gp

# ---------------------------------------------------------------------------- #

# Main function -------------------------------------------------------------- #

def update_dic(a_atom_center, f_score, d_to_update):
    i_posi_x = a_atom_center[0]
    i_posi_y = a_atom_center[1]
    i_posi_z = a_atom_center[2]
    key_dic = str(i_posi_x) + "_" + str(i_posi_y) + "_" + str(i_posi_z)
    d_to_update[key_dic] = f_score
