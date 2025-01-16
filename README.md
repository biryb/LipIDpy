# lipyd
A script that produces an untargeted lipid library from one or multiple DDA files

# Scope
This script will only use MS2 information from mgf files. It's recommended to use the library produced by this script in conjunction with MS1 peak picking software, and the output file is directly compatible with mzmine4.
This is just a script, not a package.

# Usage
Download the script as well as LipidMatch rule-based lipid libraries. Create a virtual environment, and install dependencies via pip, conda, or alike.

# Dependencies
Python 3.10.9
pandas 1.5.3
polars 1.5.0
pyteomics 4.6
