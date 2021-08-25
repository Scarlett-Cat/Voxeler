# --- VOXELER --- #
Bioinformatic tool using 3D matrices to perform structure analysis
Developped by Julien Lenoir, based on the work of Cyril Bigot for the Computational Drug Discovery Group at the Faculty of Pharmacy at the University of Helsinki

# --- Foreword --- #
	Before installing or running Voxeler, you need to set in your terminal your working directory inside the "Voxeler/" folder

# --- Installation --- #
	Execute the "/install.sh" script (bash install.sh)
		It creates a Python virtual environement and install all the needed packages
		Some errors may appear, caused by missing packages on your system such as "python3-venv" or "unzip"
		In case of errors, here the complete list of the Python packages required by the programm, which you can install manually (the install script may help you)
			numpy
			sklearn
			ete3
			six
			PyQt5
			psutils
		Also, be sure that the density folder at "Voxeler_code/resources/densities.zip" is unziped (a "densities" folder should exist), if not, unzip it manually

# --- Running --- #
	0 - The first step is to configure the program
		Choose the tool to run, head to "/Voxeler_code/config/global_parameters.txt"
	1 - Configure the tool
		Choose the parameters of the tool contained within "/Voxeler_code/config/[tool]_parameters.txt"
	2 - Place your inputs
		Place your PDB files you need to process in the corresponding folder, in "Voxeler_code/inout_[tool]/input/"
	3 - Lauch the program
		Everything should be fine, still at the "Voxeler" folder root, executes the "run.sh" script (bash run.sh)
