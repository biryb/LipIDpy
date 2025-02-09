# LipIDpy
A command line tool that produces an untargeted lipid library from one or multiple DDA files.

# Scope
This script will only use MS2 information from mgf files. It's recommended to use the library produced by this script in conjunction with MS1 peak picking software, and the output file is directly compatible with mzmine4. It only works for experimental conditions where mobile phase contained ammonium formate or ammonium acetate.
The script is agnostic to the mass spectrometry platform as it purely uses the detected fragments. No parameter optimization is required.

# Usage
## Install from GitHub

To install the package directly from GitHub, run the following command:
```bash
pip install git+https://github.com/biryb/LipIDpy.git
```

## Install from a Local Directory

If you'd like to install from your local repository, download the repository, navigate to the folder containing the `setup.py` file and run:
```bash
pip install .
```

# Usage

It's recommended to use a virtual environment to run LipIDpy.<br>
If using Conda, open a terminal (Mac) or command line/powershell (Windows) and run the following commands:<br>

```bash
conda create --name lipidpy python=3.11
conda activate lipidpy
```

In Mac terminal, the prompt should appear as
```bash
(lipidpy) id@macname lipidpy %                                                                         
```

Then run
```bash
pip install git+https://github.com/biryb/LipIDpy.git
```
to install LipIDpy.<br>

That's it - you're ready to process samples. The basic usage is like this:

```bash
lipidpy <path_mgf_files> <path_lipid_library>
```
path_lipid_library is a path to the lipid library (formate or acetate) downloaded via this repo!

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

