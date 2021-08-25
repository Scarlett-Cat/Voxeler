#!/bin/bash

# Activating the virtual environment
source ./venv/bin/activate

# Changing the active directory
cd Voxeler_code

# Launching the programm
python3 ./main.py

# Exiting the virtual environment
deactivate
