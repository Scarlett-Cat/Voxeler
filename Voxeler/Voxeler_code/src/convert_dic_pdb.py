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
from src.convert_grid_pdb import write_pdb_output


# ---------------------------------------------------------------------------- #

# Main function -------------------------------------------------------------- #

def thresh_positions(d_positions):
    """
    	Convert a dictionnary into an array containing only the desired value
    	:param d_positions: a dictionary which keys are positions and values occurences of each position
    	:return a_position_organized: array of the positions selected.
    """
    f_thresh = gp.D_PARAMETERS_SOLUBILIZATION[
        "f_occu_threshold"]  # The threshold of occurrences on wich a value is selected if its proportion of occurences is
    # strictly superior to this threshold
    i_number_run = gp.D_PARAMETERS_SOLUBILIZATION["i_launch_number"]
    i_len = 0
    l_position = []
    for key in d_positions:
        i_occurrence = float(d_positions[key])
        if i_occurrence / i_number_run > f_thresh:          # Only selects the positions with enough occurrences
            i_len += 1
            l_posi = key.split("_")
            i_pos_x = l_posi[0]
            i_pos_y = l_posi[1]
            i_pos_z = l_posi[2]
            l_position.append(i_pos_x)
            l_position.append(i_pos_y)
            l_position.append(i_pos_z)
            if gp.D_PARAMETERS_SOLUBILIZATION["s_type_b_factor"] is True:
                l_position.append(i_occurrence)     # Using the occurrence as B-factor
            else:
                f_score = gp.D_PDB_SCORING[key]     # Using the score as B-factor
                l_position.append(f_score)
    a_position = np.array(l_position)
    a_position_organized = a_position.reshape(i_len, 4)
    return a_position_organized
    # convert_dic_pdb(a_position_organized, d_parameters, o_system, o_structure)


def convert_dic_pdb(a_density, o_structure, d_parameters):
    """
	Convert an array into a PDB file
	:param a_density:
	:return:
	"""
    i_water_code = gp.D_ELEMENT_NUMBER["OOW"]

    # STEP 0 : Preparing variables ---------------------- #
    p_output = "{}/{}_{}.pdb".format(
        d_parameters["p_output_solubilization"],
        "randomax",
        o_structure.s_name
    )
    # END STEP 0 ---------------------------------------- #

    l_water_density = []

    for a_water in a_density:
        l_line_buffer = []

        l_line_buffer.append("HETATM")  # The element type
        l_line_buffer.append(65535 - i_water_code)  # The code of the water molecules
        l_line_buffer.append("OOW")  # The code of the water molecule
        l_line_buffer.append("")  # Possible alternative location
        l_line_buffer.append("HOH")  # Name of the residue
        l_line_buffer.append("")  # Water molecules are not part of a chain
        l_line_buffer.append(9999 - i_water_code)  # Serial number of the residue
        l_line_buffer.append("")  # Code for the insertion of a residue
        l_line_buffer.append(a_water[0])  # The occupancy of this atom at this position
        l_line_buffer.append(a_water[1])  # The occupancy of this atom at this position
        l_line_buffer.append(a_water[2])  # The occupancy of this atom at this position
        l_line_buffer.append(1.0)  # The occupancy of this atom at this position
        l_line_buffer.append(a_water[3])  # The occupancy of this atom at this position
        l_line_buffer.append("O")  # The symbol of the element
        l_line_buffer.append(0)  # No charge for water

        l_water_density.append(l_line_buffer)

    write_pdb_output(
        o_structure=o_structure,
        l_water_molecules=l_water_density,
        p_output=p_output
    )
