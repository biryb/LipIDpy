# lipyd
A script that produces an untargeted lipid library from one or multiple DDA files

# Scope
This script will only use MS2 information from mgf files. It's recommended to use the library produced by this script in conjunction with MS1 peak picking software, and the output file is directly compatible with mzmine4. It only works for experimental conditions where mobile phase contained ammonium formate or ammonium acetate.
The script is agnostic to the mass spectrometry platform as it purely uses the detected fragments. No parameter optimization is required.
This is just a script, not a package.

# Usage
Download the script as well as LipidMatch rule-based lipid libraries. Create a virtual environment, and install dependencies via pip, conda, or alike.

# Setup
After downloading the files, set up a folder called lipyd somewhere in the file system.</br>
The folder is supposed to look like this:</br>

![image](https://github.com/user-attachments/assets/e1aafd89-6197-4527-9dde-8dc06c1b060f)


# Parameters
In rows 14-19 of the script, the user has to set a four parameters. If you're running an Orbitrap and your mobile phase contains ammonium formate, the only parameter you
need to change is the file path (dir_base).</br>
MP_additive - ammonium salt used in the analysis. This can be either "Formate" or "Acetate". A string value.</br>
polarity - POS or NEG</br>
dir_base - the absolute path to the folder where the script and data files are wrapped in r'path' (e.g. r'/Volumes/name/lipyd/')</br>
mztol - this is to tolerance for mz in MS2 spectra, and dictates how much a detected fragment is allowed to deviate from the known exact value. I usually set this to 3.5 mmu (my instrument is QE HF-X at 60,000 resolution for MS2), but anything between 1.5-10 might be appropriate

# Settings for MSConvert
You can convert the raw mass spec data to .mgf's via MSConvert:
![image](https://github.com/user-attachments/assets/41167f3e-1844-4fc6-afae-ad7680747e7d)

# Dependencies
LipidMatch lipid libraries (supplied as a part of this repository; credit to https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1744-3)<br/> 
Python 3.10.9 <br/>
polars 1.20.0 <br/>
pyteomics 4.6 <br/>

