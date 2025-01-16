# lipyd
A script that produces an untargeted lipid library from one or multiple DDA files

# Scope
This script will only use MS2 information from mgf files. It's recommended to use the library produced by this script in conjunction with MS1 peak picking software, and the output file is directly compatible with mzmine4.
The script is agnostic to the mass spectrometry platform as it purely uses the detected fragments. No parameter optimization is required.
This is just a script, not a package.

# Usage
Download the script as well as LipidMatch rule-based lipid libraries. Create a virtual environment, and install dependencies via pip, conda, or alike.

# Settings for MSConvert
You can convert the raw mass spec data to .mgf's via MSConvert:
![image](https://github.com/user-attachments/assets/41167f3e-1844-4fc6-afae-ad7680747e7d)

# Dependencies
LipidMatch lipid libraries.<br/> 
Python 3.10.9 <br/>
pandas 1.5.3 <br/>
polars 1.5.0 <br/>
pyteomics 4.6 <br/>

