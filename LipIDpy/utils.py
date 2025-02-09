import os
import polars as pl
from pyteomics import mgf
import numpy as np

def read_or_create_parameters(params_file):
    """Reads parameters from file if it exists; otherwise, creates the file with default values."""
    default_params = {
        "MP_additive": "Formate",  # Formate or Acetate
        "polarity": "POS",  # POS or NEG
        "mztol": 0.0035  # Tolerance in mmu
    }
    
    if os.path.exists(params_file):
        parameters = {}
        with open(params_file, "r") as f:
            for line in f:
                key, value = line.strip().split("=")
                if key == "mztol":
                    parameters[key] = float(value)
                else:
                    parameters[key] = value
        return parameters
    else:
        with open(params_file, "w") as f:
            for key, value in default_params.items():
                f.write(f"{key}={value}\n")
        return default_params

def merge_lipid_libraries(list_of_library_files,dir_library):
    exclude=['Precursor_Library_POS.csv','LIPID_ID_CRITERIA.csv','Precursor_Library_NEG.csv']
    files_sel=  [x for x in list_of_library_files if x not in exclude]
    df_alllipids = pl.DataFrame()
    for file in files_sel:
        if file not in exclude and file[0]!='~' and file[0]!='.' and '.csv' in file:
            df_tmp = pl.read_csv(os.path.join(dir_library, file))
            df_tmp = df_tmp.with_columns(pl.lit( df_tmp[df_tmp.columns[0]].str.strip_chars('+Na').str.strip_chars('+H').str.strip_chars('-H').str.strip_chars('+NH4')).alias('lipid'))
            df_tmp = df_tmp.with_columns(pl.lit( df_tmp.columns[1]).alias('adduct'))
            df_tmp = df_tmp.with_columns(pl.lit( df_tmp[df_tmp.columns[1]].alias('precursor_mz')))
            fragment_names = df_tmp.columns[2:-3]
            if len(fragment_names)>0:
                dict_og_cols = dict(zip(fragment_names,[x+3 for x in range(0,len(fragment_names))]))
                df_tmp = df_tmp.drop(df_tmp.columns[:2])
                df_tmp = df_tmp.unpivot(index=['lipid','adduct','precursor_mz'],on=fragment_names,variable_name='Fragment',value_name='fragment_mz')
                df_tmp = df_tmp.with_columns(pl.lit( file).alias('source_file'))
                df_tmp = df_tmp.with_columns(df_tmp['Fragment'].replace(dict_og_cols).alias('original_fragment_column').cast(int))
                df_tmp = df_tmp.with_columns(df_tmp['fragment_mz'].cast(float).alias('fragment_mz'))
                df_alllipids = pl.concat([df_alllipids,df_tmp])
    df_alllipids = df_alllipids.filter(df_alllipids['adduct'].str.contains('\+'))
    return df_alllipids

def define_rules(rule_file):
    df_rules = pl.read_excel(rule_file)
    return df_rules.rename({'CSV File': 'source_file'})

def initial_matching(file_name, df_lipid_library, mztol=0.0035):
    dda_file = mgf.read(file_name)
    df_foundlipids = pl.DataFrame()

    for elem in dda_file:
        rt_s = elem['params']['rtinseconds']
        precursor_mz = elem['params']['pepmass'][0]
        scan_number = int(elem['params']['title'].split('scan=')[1].strip('"'))
        df_possible_lipids = df_lipid_library.filter(
            pl.col("precursor_mz").is_between(precursor_mz - mztol, precursor_mz + mztol)
        )

        if df_possible_lipids.shape[0] > 0:
            df_possible_lipids = df_possible_lipids.with_columns(pl.lit(scan_number).alias('scan_number'))
            df_possible_lipids = df_possible_lipids.with_columns(pl.lit(rt_s).alias('rt_s'))
            df_foundlipids = pl.concat([df_foundlipids, df_possible_lipids])

    return df_foundlipids

def match_dda_to_alllipids(dda_file_names, dir_data, df_alllipids, mztol=0.0035):
    df_foundlipids = pl.DataFrame()

    for dda_file_name in dda_file_names:
        file_name = os.path.join(dir_data, dda_file_name)
        df_tmp = initial_matching(file_name, df_alllipids, mztol)
        df_foundlipids = pl.concat([df_foundlipids, df_tmp])

    return df_foundlipids

def transform_fragment_inds_to_int(df):
    df = df.with_columns(
        pl.col("Columns containing all necessary fragments for ID").str.split(";").cast(pl.List(pl.Int64))
        .alias("All_essential_fragments"))
    df = df.with_columns(pl.col("Columns containing fragments where atleast one fragment must be observed for confirmation")
        .str.split(";").cast(pl.List(pl.Int64)).alias("Qualifying_fragments"))
    return df

def mark_lipid_match_rows(df):
    df = df.with_columns(
        pl.when(pl.col('original_fragment_column').is_in(pl.col('All_essential_fragments')))
        .then(1).otherwise(0).alias('Essential_fragment'))
    df = df.with_columns(pl.when(pl.col('original_fragment_column').is_in(pl.col('Qualifying_fragments')))
        .then(1).otherwise(0).alias('Qualifying_fragment'))
    return df

def prepare_lipid_data_with_rules(df, df_rules):
    df = df.unique()
    df = df.with_columns(pl.concat_str([df['scan_number'], df["lipid"], df["adduct"], df["Fragment"], df["rt_s"]], separator="_").alias("row_ID"))
    df = df.join(df_rules, on='source_file', how='left')
    df = transform_fragment_inds_to_int(df)
    df = mark_lipid_match_rows(df)
    df = df.with_columns(n_essential_fragments=pl.col("All_essential_fragments").list.len())
    df = df.with_columns(n_qualifying_fragments=pl.col("Qualifying_fragments").list.len())
    return df

def aggregate_and_filter_lipids_by_rules(df):
    df = df.unique().group_by(['scan_number', 'lipid', 'adduct', 'rt_s']).agg([
        pl.sum('Essential_fragment').alias('n_essential_fragments_detected'),
        pl.sum('Qualifying_fragment').alias('n_qualifying_fragments_detected'),
        pl.mean('precursor_mz'),
        pl.mean('n_essential_fragments'),
        pl.col("Fragment").str.concat("|").alias('Fragments_detected')
    ])
    df = df.with_columns((df['n_essential_fragments_detected'] / df['n_essential_fragments']).alias('Fraction_essential_detected'))
    return df.filter(pl.col('Fraction_essential_detected') == 1)

def save_results(df, dir_output):
    os.makedirs(dir_output, exist_ok=True)
    df.write_csv(os.path.join(dir_output, 'rulebased_lipids.csv'))
    df_mzmine = df.rename({'lipid': 'name', 'rt_s': 'rt', 'precursor_mz': 'mz'})[['name', 'rt', 'mz', 'adduct']]
    df_mzmine.write_csv(os.path.join(dir_output, 'rulebased_lipids_mzmine.csv'))
