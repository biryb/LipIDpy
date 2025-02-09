import os
import argparse
import polars as pl
from pyteomics import mgf
import utils


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process MGF files in a given folder and save results.")
    parser.add_argument("mgf_folder", type=str, help="Path to the folder containing .mgf files")
    parser.add_argument("library_folder", type=str, help="Path to the folder containing lipid library files")
    args = parser.parse_args()

    # Validate input folder
    if not os.path.isdir(args.mgf_folder):
        print(f"Error: The folder '{args.mgf_folder}' does not exist.")
        return
    print(f"Looking for mgf in {args.mgf_folder}")
    # Define output directory inside the provided folder
    dir_output = os.path.join(args.mgf_folder, "output")
    os.makedirs(dir_output, exist_ok=True)

    # Path to parameters file
    params_file = os.path.join(dir_output, "parameters.txt")

    # Read or create parameters
    parameters = utils.read_or_create_parameters(params_file)

    # Extract parameter values
    MP_additive = parameters["MP_additive"]
    polarity = parameters["polarity"]
    mztol = parameters["mztol"]

    print(f"Using parameters: MP_additive={MP_additive}, polarity={polarity}, mztol={mztol}")

    # Set up library directory based on MP_additive
    
    dir_library = args.library_folder
    '''
    dir_acetate_libraries = os.path.join(args.library_folder, "LipidMatch_Libraries_Acetate")
    dir_formate_libraries = os.path.join(args.library_folder, "LipidMatch_Libraries_Formate")

    if MP_additive.upper() == 'ACETATE':
        dir_library = dir_acetate_libraries
    elif MP_additive.upper() == 'FORMATE':
        dir_library = dir_formate_libraries
    else:
        print("Mobile phase additive ('MP_additive') must be formate or acetate")
        return
    '''
    # Get relevant library files
    library_files = os.listdir(args.library_folder)
    if polarity == 'NEG':
        files_sel = [x for x in library_files if "NEG" in x.upper() and x[0] != '.']
    else:
        files_sel = [x for x in library_files if "NEG" not in x.upper() and x[0] != '.']
        files_sel = [x for x in files_sel if x != 'PEG_Pos.csv']  # Exclude PEG contaminants

    # Merge LipidMatch libraries into one dataframe
    df_alllipids = utils.merge_lipid_libraries(files_sel, dir_library)

    # Load rules for lipid matching
    df_rules = utils.define_rules(os.path.join(args.library_folder, "LIPID_ID_CRITERIA.xlsx"))

    # Process MGF files
    mgf_files = [x for x in os.listdir(args.mgf_folder) if x.endswith(".mgf")]
    if not mgf_files:
        print("No .mgf files found in the given folder.")
        return

    df_foundlipids = utils.match_dda_to_alllipids(mgf_files, args.mgf_folder, df_alllipids, mztol)
    df_foundlipids = utils.prepare_lipid_data_with_rules(df_foundlipids,df_rules)

    # Aggregate and filter lipids based on rules
    print(df_foundlipids.head())
    print(df_foundlipids.columns)
    df_foundlipids = utils.aggregate_and_filter_lipids_by_rules(df_foundlipids)

    # Save results
    utils.save_results(df_foundlipids,dir_output)

if __name__ == "__main__":
    main()
