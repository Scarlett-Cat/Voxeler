# Information ---------------------------------------------------------------- #
# Author :	Anna Spampinato
# Github :	Scarlett-Cat
# Created : February 2021
# Updated :
# ---------------------------------------------------------------------------- #



# Importations --------------------------------------------------------------- #

import numpy as np
# Allows Numpy array manipulation

# Classes -------------------------------------------------------------------- #

# General library

# ---------------------------------------------------------------------------- #



# Main function -------------------------------------------------------------- #

def update_dictionary(a_positions, d_water_density):
    """
    Updates the dictionary of occurrences
    :param a_positions: The array containing the positions of water molecules
    """
    for item in a_positions:
        i_xpos = round(item[0], 1)          # Saves the coordinates of the position
        i_ypos = round(item[1], 1)
        i_zpos = round(item[2], 1)
        s_pos = str(i_xpos) + "_" + str(i_ypos) + "_" + str(i_zpos)         # Convert the position in string
        if s_pos in d_water_density:            # Uses the position as a key
            d_water_density[s_pos] += 1         # Uses the occurrences as a value
        else:
            d_water_density[s_pos] = 1
    #return(d_water_density)