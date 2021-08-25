#!/bin/bash


# Installing packages
echo "Installing packages"
{
	pip3 install numpy
	pip3 install psutil
	pip3 install ete3
	pip3 install PyQt5==5.11.3
	pip3 install six
	pip3 install sklearn
} || {
	echo "Some packages may not be correctly installed"
	echo "Here the packages currently installed :"
	pip list
	echo "numpy, psutil, ete3, PyQt5==5.11.3 and sklearn needs to be installed"
	echo "You can install the packages manually by running 'source venv/bin/activate' and 'pip3 install [package_name]'"
}
