# Information ---------------------------------------------------------------- #
# Author :	Julien Lenoir
# Github :	Blackounet
# Created : June 2020
# Updated :
# ---------------------------------------------------------------------------- #


# Importations --------------------------------------------------------------- #

import os


# Allows file system operations

# ---------------------------------------------------------------------------- #


# Main function -------------------------------------------------------------- #

def init():
    """
	Initializes global variables used all across the software
	"""

    # STEP 0 : Preparing variables ---------------------- #
    global O_SYSTEM_COMPARISON  # Object containing the whole comparison system
    global O_SYSTEM_SOLUBILIZATION  # Object containing the whole solubilization system

    global D_PARAMETERS_GLOBAL  # Dictionary of the program main parameters
    global D_PARAMETERS_COMPARISON  # Dictionary of the parameters specific to the structure comparison
    global D_PARAMETERS_SOLUBILIZATION  # Dictionary of the parameters specific to the structure solubilization

    global D_EXPECTED_PARAMETERS_GLOBAL  # Dictionary of the expected parameters for the main program
    global D_EXPECTED_PARAMETERS_COMPARISON  # Dictionary of the expected parameters for the structure comparison
    global D_EXPECTED_PARAMETERS_SOLUBILIZATION  # Dictionary of the expected parameters for the structure solubilization

    global D_ELEMENT_NUMBER  # Dictionary of atomic number associated to the atom type
    global D_NUMBER_ELEMENT  # Dictionary of atom type associated to the atomic number

    global D_WATER_POSITION  # Dictionary of positions of the water
    global D_WATER_SCORING   # Dictionary of the association grid position - score
    global D_PDB_SCORING     # Dictionary of the association pdb position - score
    # END STEP 0 ---------------------------------------- #

    # STEP 1 : Initializing variables ------------------- #
    O_SYSTEM_COMPARISON = None
    O_SYSTEM_SOLUBILIZATION = None

    D_PARAMETERS_GLOBAL = {}
    D_PARAMETERS_COMPARISON = {}
    D_PARAMETERS_SOLUBILIZATION = {}

    D_WATER_POSITION = {}
    D_WATER_SCORING = {}
    D_PDB_SCORING = {}

    # Contains the parameter's name, it's key, it's type and it's default value
    D_EXPECTED_PARAMETERS_GLOBAL = {

        # Base parameters
        "path_to_logs": ["p_log", "path", "project.log"],  # Not shown in global_parameters.txt
        "cpu_allocated": ["i_cpu_allocated", "int", "None"],
        "memory_allocated": ["f_memory_allocated", "float", "4.0"],

        # Features requested
        "run_comparison": ["b_run_comparison", "bool", "True"],
        "run_solubilization": ["b_run_solubilization", "bool", "True"],

        # Parsing parameters
        "path_to_global_parameters": ["p_global_parameters", "path", "config/global_parameters.txt"],
        # Not shown in global_parameters.txt
        "path_to_comparison_parameters": ["p_comparison_parameters", "path", "config/comparison_parameters.txt"],
        "path_to_solubilization_parameters": ["p_solubilization_parameters", "path",
                                              "config/solubilization_parameters.txt"],
        "comment_delimiters_string": ["l_s_comment_delimiters", "list_str", ["# ", "\"\"\" ", "// "]],
        # Not shown in global_parameters.txt
        "comment_delimiters_char": ["l_c_comment_delimiters", "list_char", ["#", "\"", "/"]],
        # Not shown in global_parameters.txt
    }
    D_EXPECTED_PARAMETERS_COMPARISON = {

        # Input paths
        "path_to_PDB_directory": ["p_input_pdb", "path", "comparison/input/"],

        # Output paths
        "path_to_output_directory": ["p_output_comparison", "path", "output/comparison/"],

        # Comparison parameters
        "distance_between_points": ["f_grid_spacing", "float", "0.1"],
        "consider_elements": ["b_consider_elements", "bool", "True"],
        "comparison_normalisation": ["s_comparison_normalisation", "str", "Max"],
        "grid_geometry": ["s_grid_geometry", "str", "Sphere"],

        # Tree generation
        "tree_name": ["s_tree_name", "str", "None"],
        "show_tree": ["b_show_tree", "bool", "True"],
        "tree_shape": ["s_tree_shape", "str", "Classic"],
        "only_topology": ["b_only_topology", "bool", "True"],
        "leaf_name": ["b_leaf_name", "bool", "True"],
        "branch_length": ["b_branch_length", "bool", "True"],

        # Parsing parameters
        "discard_atom": ["b_discard_atom", "bool", "False"],
        "discard_hetatm": ["b_discard_hetatm", "bool", "False"],
        "discard_hydrogen": ["b_discard_hydrogen", "bool", "False"],
        "discard_water": ["b_discard_water", "bool", "False"],
        "discard_alternative": ["b_discard_alternative", "bool", "False"],
        "chain_white_list": ["l_c_chain_white", "list_char", ""],
        "chain_black_list": ["l_c_chain_black", "list_char", ""],
        "residue_white_list": ["l_s_residue_white", "list_str", ""],
        "residue_black_list": ["l_s_residue_black", "list_str", ""],
        "residue_id_white_list": ["l_i_residue_white", "list_int", ""],
        "residue_id_black_list": ["l_i_residue_black", "list_int", ""],
        "atom_white_list": ["l_s_atom_white", "list_str", ""],
        "atom_black_list": ["l_s_atom_black", "list_str", ""],
    }
    D_EXPECTED_PARAMETERS_SOLUBILIZATION = {

        # Input paths
        "path_to_PDB_directory": ["p_input_pdb", "path", "comparison/input/"],

        # Output paths
        "path_to_output_directory": ["p_output_solubilization", "path", "output/solubilization/"],

        # Solubilization parameters
        "distance_between_points": ["f_grid_spacing", "float", "0.1"],
        "grid_geometry": ["s_grid_geometry", "str", "Sphere"],
        "use_randomax": ["b_use_randomax", "bool", "False"],
        "launch_number": ["i_launch_number", "int", "1"],
        "occurrences_threshold": ["f_occu_threshold", "float", "O.5"],
        "solubilization_method": ["s_solubilization_method", "str", "Itermax"],
        "random_threshold": ["f_random_threshold", "float", "1.0"],
        "solubilization_radius": ["l_f_solubilization_radius", "list_float", "1.0, 2.0"],
        "max_neighbor_distance": ["f_max_neighbor_distance", "float", "5.0"],
        "max_neighbor_number": ["i_max_neighbor_number", "int", "10"],
        "scoring_metric": ["s_scoring_metric", "str", "euclidean"],
        "only_first_neighbor": ["b_only_first_neighbor", "bool", "False"],
        "scoring_per_residue": ["b_scoring_per_residue", "bool", "False"],
        "mean_score": ["b_mean_score", "bool", "True"],
        "min_water_score": ["f_min_water_score", "float", "0.5"],

        # Atom scoring
        "electronic_densities_folder": ["p_electronic_densities", "path", "resources/densities"],
        "normalize_electronic_densities": ["b_normalize_densities", "bool", "False"],
        "occurrence_b_factor": ["s_type_b_factor", "bool", "False"],

        # Parsing parameters
        "discard_atom": ["b_discard_atom", "bool", "False"],
        "discard_hetatm": ["b_discard_hetatm", "bool", "False"],
        "discard_hydrogen": ["b_discard_hydrogen", "bool", "False"],
        "discard_water": ["b_discard_water", "bool", "False"],
        "discard_alternative": ["b_discard_alternative", "bool", "False"],
        "chain_white_list": ["l_c_chain_white", "list_char", ""],
        "chain_black_list": ["l_c_chain_black", "list_char", ""],
        "residue_white_list": ["l_s_residue_white", "list_str", ""],
        "residue_black_list": ["l_s_residue_black", "list_str", ""],
        "residue_id_white_list": ["l_i_residue_white", "list_int", ""],
        "residue_id_black_list": ["l_i_residue_black", "list_int", ""],
        "atom_white_list": ["l_s_atom_white", "list_str", ""],
        "atom_black_list": ["l_s_atom_black", "list_str", ""],
    }

    D_ELEMENT_NUMBER = {  # Dictionary of atomic numbers
        "NX": 0,
        "H": 1,
        "HE": 2,
        "LI": 3,
        "BE": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "NE": 10,
        "NA": 11,
        "MG": 12,
        "AL": 13,
        "SI": 14,
        "P": 15,
        "S": 16,
        "CL": 17,
        "AR": 18,
        "K": 19,
        "CA": 20,
        "SC": 21,
        "TI": 22,
        "V": 23,
        "CR": 24,
        "MN": 25,
        "FE": 26,
        "CO": 27,
        "NI": 28,
        "CU": 29,
        "ZN": 30,
        "GA": 31,
        "GE": 32,
        "AS": 33,
        "SE": 34,
        "BR": 35,
        "KR": 36,
        "RB": 37,
        "SR": 38,
        "Y": 39,
        "ZR": 40,
        "NB": 41,
        "MO": 42,
        "TC": 43,
        "RU": 44,
        "RH": 45,
        "PD": 46,
        "AG": 47,
        "CD": 48,
        "IN": 49,
        "SN": 50,
        "SB": 51,
        "TE": 52,
        "I": 53,
        "XE": 54,
        "CS": 55,
        "BA": 56,
        "LA": 57,
        "CE": 58,
        "PR": 59,
        "ND": 60,
        "PM": 61,
        "SM": 62,
        "EU": 63,
        "GD": 64,
        "TB": 65,
        "DY": 66,
        "HO": 67,
        "ER": 68,
        "TM": 69,
        "YB": 70,
        "LU": 71,
        "HF": 72,
        "TA": 73,
        "W": 74,
        "RE": 75,
        "OS": 76,
        "IR": 77,
        "PT": 78,
        "AU": 79,
        "HG": 80,
        "TL": 81,
        "PB": 82,
        "BI": 83,
        "PO": 84,
        "AT": 85,
        "RN": 86,
        "FR": 87,
        "RA": 88,
        "AC": 89,
        "TH": 90,
        "PA": 91,
        "U": 92,
        "NP": 93,
        "PU": 94,
        "AM": 95,
        "CM": 96,
        "BK": 97,
        "CF": 98,
        "ES": 99,
        "FM": 100,
        "MD": 101,
        "NO": 102,
        "LR": 103,
        "RF": 104,
        "DB": 105,
        "SG": 106,
        "BH": 107,
        "HS": 108,
        "MT": 109,
        "META": 110,
        "OOW": 111
    }
    D_NUMBER_ELEMENT = {  # Dictionary of atomic symbols
        0: "NX",
        1: "H",
        2: "HE",
        3: "LI",
        4: "BE",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "NE",
        11: "NA",
        12: "MG",
        13: "AL",
        14: "SI",
        15: "P",
        16: "S",
        17: "CL",
        18: "AR",
        19: "K",
        20: "CA",
        21: "SC",
        22: "TI",
        23: "V",
        24: "CR",
        25: "MN",
        26: "FE",
        27: "CO",
        28: "NI",
        29: "CU",
        30: "ZN",
        31: "GA",
        32: "GE",
        33: "AS",
        34: "SE",
        35: "BR",
        36: "KR",
        37: "RB",
        38: "SR",
        39: "Y",
        40: "ZR",
        41: "NB",
        42: "MO",
        43: "TC",
        44: "RU",
        45: "RH",
        46: "PD",
        47: "AG",
        48: "CD",
        49: "IN",
        50: "SN",
        51: "SB",
        52: "TE",
        53: "I",
        54: "XE",
        55: "CS",
        56: "BA",
        57: "LA",
        58: "CE",
        59: "PR",
        60: "ND",
        61: "PM",
        62: "SM",
        63: "EU",
        64: "GD",
        65: "TB",
        66: "DY",
        67: "HO",
        68: "ER",
        69: "TM",
        70: "YB",
        71: "LU",
        72: "HF",
        73: "TA",
        74: "W",
        75: "RE",
        76: "OS",
        77: "IR",
        78: "PT",
        79: "AU",
        80: "HG",
        81: "TL",
        82: "PB",
        83: "BI",
        84: "PO",
        85: "AT",
        86: "RN",
        87: "FR",
        88: "RA",
        89: "AC",
        90: "TH",
        91: "PA",
        92: "U",
        93: "NP",
        94: "PU",
        95: "AM",
        96: "CM",
        97: "BK",
        98: "CF",
        99: "ES",
        100: "FM",
        101: "MD",
        102: "NO",
        103: "LR",
        104: "RF",
        105: "DB",
        106: "SG",
        107: "BH",
        108: "HS",
        109: "MT",
        110: "META",
        111: "OOW"
    }


# END STEP 1 ---------------------------------------- #


# STEP 2 : ------------------------------------------ #
# END STEP 2 ---------------------------------------- #

# ---------------------------------------------------------------------------- #


# Auxiliary functions -------------------------------------------------------- #

def loads_default_parameters():
    D_PARAMETERS_GLOBAL["p_log"] = os.path.abspath(D_EXPECTED_PARAMETERS_GLOBAL["path_to_logs"][2])
    D_PARAMETERS_GLOBAL["p_global_parameters"] = os.path.abspath(
        D_EXPECTED_PARAMETERS_GLOBAL["path_to_global_parameters"][2])
    D_PARAMETERS_GLOBAL["l_s_comment_delimiters"] = D_EXPECTED_PARAMETERS_GLOBAL["comment_delimiters_string"][2]
    D_PARAMETERS_GLOBAL["l_c_comment_delimiters"] = D_EXPECTED_PARAMETERS_GLOBAL["comment_delimiters_char"][2]

# ---------------------------------------------------------------------------- #


# Reference ------------------------------------------------------------------ #

# Importation
# from import global_parameters as gp
# Contains the global variables

# Usage
# init()		# Initializes the global variables

# ---------------------------------------------------------------------------- #
