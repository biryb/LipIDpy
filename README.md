# LipIDpy
A command line tool that produces an untargeted lipid library from one or multiple DDA files.

# Scope
This script will only use MS2 information from mgf files. It's recommended to use the library produced by this script in conjunction with MS1 peak picking software, and the output file is directly compatible with mzmine4. It only works for experimental conditions where mobile phase contained ammonium formate or ammonium acetate.
The script is agnostic to the mass spectrometry platform as it purely uses the detected fragments. No parameter optimization is required.
This is just a script, not a package or a command line tool. The script can be run in an IDE. Download it and open in VSCode or Spyder or alike.

# Usage

First, ensure you have Poetry installed. If not, install it:<br>
Windows:
```bash
winget install Poetry.Poetry
```
Mac/Linux:
```bash
curl -sSL https://install.python-poetry.org | python3 -
```

Then, clone the repository:
```bash
poetry add git+https://github.com/biryb/lipidpy.git
```
```bash
poetry run lipidpy <path_mgf_files> <path_library_files>
```


# Parameters
In rows 14-19 of the script, the user has to set four parameters. If you're running an Orbitrap and your mobile phase contains ammonium formate, the only parameter you
need to change is the file path (dir_base).</br>
MP_additive - ammonium salt used in the analysis. This can be either "Formate" or "Acetate". A string value.</br>
polarity - POS or NEG</br>
dir_base - the absolute path to the folder where the script and data files are wrapped in r'path' (e.g. r'/Volumes/name/LipIDpy/')</br>
mztol - this is the tolerance for mz in MS2 spectra, and dictates how much a detected fragment is allowed to deviate from the known exact value. I usually set this to 3.5 mmu (my instrument is QE HF-X at 60,000 resolution for MS2), but anything between 1.5-10 might be appropriate

# Settings for MSConvert
You can convert the raw mass spec data to .mgf's via MSConvert:
![image](https://github.com/user-attachments/assets/41167f3e-1844-4fc6-afae-ad7680747e7d)

# Dependencies
LipidMatch lipid libraries (supplied as a part of this repository; credit to https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1744-3)<br> 
Python 3.11.7 <br>
Poetry 2.0.1<br>

# Contact
If you like this project, have questions/ideas, or would to chat, reach out to birgittaryback@gmail.com

