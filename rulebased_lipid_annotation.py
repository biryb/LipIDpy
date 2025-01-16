#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 10:32:31 2024

@author: imac9
"""

import os
import polars as pl
import pyteomics
from pyteomics import mgf

#%% Set up parameters

MP_additive = 'Formate' # Formate or Acetate
polarity = 'POS' # POS or NEG
dir_base = r'/Volumes/bryback/_research_and_development/lipyd/' # replace this with the absolute path to where the script is located, e.g. r'/Volumes/name/lipyd/'
mztol=0.0035 # Tolerance for each MS2 fragment to deviate from the predicted fragment m/z. The unit is mmu, not ppm. So this is 3.5mmu. Normal range is probably 1.5-10 mmu.

#%% Setting up all the relative paths needed (you shouldn't have to touch this part)
dir_analysis = os.path.join(dir_base, "analysis")
dir_data = os.path.join(dir_base, "data")
dir_acetate_libraries  =  os.path.join(dir_base, "LipidMatch_Libraries_Acetate")
dir_formate_libraries  =  os.path.join(dir_base, "LipidMatch_Libraries_Formate")

if MP_additive.upper() == 'ACETATE':
    dir_library = dir_acetate_libraries
elif MP_additive.upper() == 'FORMATE':
    dir_library = dir_formate_libraries

else:
    print("Mobile phase additive ('MP_additive') must be formate or acetate")
    
library_files = os.listdir(dir_library)


#%% Just run this; figures out what library files to use
if polarity=='NEG':
    files_sel = [x for x in library_files if "NEG"  in x.upper() and x[0]!='.']

else:
    files_sel = [x for x in library_files if "NEG"  not in x.upper() and x[0]!='.']
    files_sel = [x for x in files_sel if x!='PEG_Pos.csv'] # Exclude PEG, these are contaminants and the script always finds a lot of these
    
#%% Just run this; functions that will be used below


def merge_lipid_libraries(list_of_library_files,dir_library):
    exclude=['Precursor_Library_POS.csv','LIPID_ID_CRITERIA.csv','Precursor_Library_NEG.csv']

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
                df_tmp = df_tmp.with_columns(df_tmp['Fragment'].replace(dict_og_cols).alias('original_fragment_column'))
                df_tmp = df_tmp.with_columns(df_tmp['fragment_mz'].cast(float).alias('fragment_mz'))
                df_alllipids = pl.concat([df_alllipids,df_tmp])
    df_alllipids = df_alllipids.filter(df_alllipids['adduct'].str.contains('\+'))
    return df_alllipids

def define_rules(rule_file):
    df_rules = pl.read_excel(rule_file)
    df_rules = df_rules.rename({'CSV File':'source_file'})
    return df_rules

def initial_matching(file_name, df_lipid_library,mztol = 0.0035):
    dda_file = mgf.read(file_name)
    counter_match = 0
    counter_elems = 0
    
    df_foundlipids = pl.DataFrame()
        
    for elem in dda_file:
        precursor_mz=elem['params']['pepmass'][0]
        scan_number = elem['params']['title'].split('scan=')[1].strip('"')
        counter_elems+=1
        df_possible_lipids = df_alllipids.filter(pl.col("precursor_mz").is_between(precursor_mz-mztol,precursor_mz+mztol))
        df_possible_lipids = df_possible_lipids.with_columns(pl.lit(0).alias('fragment_ID_found'))
        df_possible_lipids_save = pl.DataFrame()
        if df_possible_lipids.shape[0]>0:
            for mz_dda in elem['m/z array']:
                df_match = df_possible_lipids.filter(pl.col("fragment_mz").is_between(mz_dda-mztol,mz_dda+mztol))
                if df_match.shape[0]>0 and df_match.shape[0]<100:
                    df_match_tmpsave=df_match
                    df_possible_lipids_tmpsave=df_possible_lipids
                    mz_dda_tmpsave=mz_dda
                    df_possible_lipids = df_possible_lipids.update(
                                               pl.select(fragment_mz = df_match['fragment_mz'], fragment_ID_found = 1),
                                               on = "fragment_mz",
                                               how = "outer"
                                            )
                    counter_match+=1
                    df_possible_lipids_save=df_possible_lipids.with_columns(pl.lit(precursor_mz).alias('precursor_mz'))
                    df_possible_lipids_save=df_possible_lipids_save.with_columns(pl.lit(elem['params']['rtinseconds']).alias('rt_s'))
        
        df_possible_lipids_save=df_possible_lipids_save.with_columns(pl.lit(int(scan_number)).alias('scan_number'))
        if df_possible_lipids_save.shape[0]>1:
            df_foundlipids = pl.concat([df_foundlipids,df_possible_lipids_save])
    return df_foundlipids

def match_dda_to_alllipids(dda_file_names, dir_data, df_alllipids, mztol=mztol):
    df_foundlipids = pl.DataFrame()
    if len(dda_file_names)==0:
        print("No valid DDA files found")
    for i,dda_file_name in enumerate(dda_file_names):
        print('Matching %s'%dda_file_name)
        file_name = os.path.join(dir_data, dda_file_name)
        df_foundlipids_tmp = initial_matching(file_name, df_alllipids,mztol=mztol)
        df_foundlipids_tmp = df_foundlipids_tmp.with_columns(df_foundlipids_tmp['scan_number']+'_%s'%(str(i)))
        if df_foundlipids.shape[0]==0:
            df_foundlipids = df_foundlipids_tmp
        else:
            df_foundlipids = pl.concat([df_foundlipids,df_foundlipids_tmp])
    return df_foundlipids


def transform_fragment_inds_to_int(df):
    df = df.with_columns(
        pl.col("Columns containing all necessary fragments for ID")
        .str.split(";")  # Split by semicolon into a list of strings
        .cast(pl.List(pl.Int64))
        .alias("All_essential_fragments")  # Name the new column
        )
        
    df = df.with_columns(
        pl.col("Columns containing fragments where atleast one fragment must be observed for confirmation")
        .str.split(";")  # Split by semicolon into a list of strings
        .cast(pl.List(pl.Int64))
        .alias("Qualifying_fragments")  # Name the new column
        )
    return df

def mark_lipid_match_rows(df):
    
    df = df.with_columns(
        pl.when(pl.col('original_fragment_column').is_in(
        pl.col('All_essential_fragments'))).
        then(1)
        .otherwise(0)
        .alias('Essential_fragment'))
    
    df = df.with_columns(
        pl.when(pl.col('original_fragment_column').is_in(
        pl.col('Qualifying_fragments'))).
        then(1)
        .otherwise(0)
        .alias('Qualifying_fragment'))

    return df

def prepare_lipid_data_with_rules(df):
    # Filter the results; remove duplicate entries and create a unique ID for each remaining lipid
    df = df.unique()
    df = df.with_columns(pl.lit(df['scan_number']+"_"+df["lipid"]+
                                                        "_"+df["adduct"]+"_"+df["Fragment"]+"_"+df["rt_s"]).alias("row_ID"))
    
    # add the rules to the found lipids so that in the future we can filter out rows that don't fulfill the rules
    df = df.join(df_rules,on='source_file')
    df = df.drop('source_file')
    
    # Transform the column indices from strings to lists of integers
    #format_fragment_info_columns(df_foundlipids_rules)
    
    # change the datatype of the original fragment column to integer
    df = df.with_columns(pl.lit(df['original_fragment_column'].cast(int).alias('original_fragment_column')))
    
    # Create a column which gets value 1 if a fragment was detected and it was on the list of essential fragments for that lipid
    df = transform_fragment_inds_to_int(df)
    df = mark_lipid_match_rows(df)
    return df

def aggregate_and_filter_lipids_by_rules(df):
    # Remove columns that are irrelevant; remove duplicates
    
    df = df[['scan_number','precursor_mz','rt_s','lipid','adduct',
                                                 'Fragment','Qualifying_fragment','Essential_fragment',
                                                 'n_essential_fragments','n_qualifying_fragments']].unique()
    df = df.unique()
    
    # Aggregate the data frame so that for each scan and lipid, we no longer have every MS2 fragment
    # but the sum of essential and qualifying fragments that were detected
    
    df_agg = df.group_by(['scan_number','lipid','adduct','rt_s']).agg([pl.sum(
        'Essential_fragment').alias('n_essential_fragments_detected'),
        pl.sum('Qualifying_fragment').alias('n_qualifying_fragments_detected'), pl.mean('n_essential_fragments'),
        pl.first('n_qualifying_fragments'),
        pl.mean('precursor_mz'), pl.col("Fragment").str.concat("|").alias('Fragments_detected')])
    
    
    # Calculate the fraction of essential fragments that were detected
    df_agg = df_agg.with_columns(pl.lit(df_agg['n_essential_fragments_detected']/df_agg['n_essential_fragments'])
                                                                     .alias('Fraction_essential_detected'))
    
    
    # Only keep lipids for which 100% of essential fragments were detected
    df_agg = df_agg.filter(pl.col('Fraction_essential_detected')==1)
    
    # For each scan, keep the lipid(s) that maximize the fraction of qualifying fragments detected
    df_agg = df_agg.filter(pl.col('n_qualifying_fragments_detected') == 
                                                                         pl.col('n_qualifying_fragments_detected').max().over('scan_number'))
    # Give lipids unique names by appending the lipid names with the adduct name
    df_agg = df_agg.with_columns(pl.lit(df_agg['lipid']+" "+df_agg['adduct']).alias('lipid-adduct'))
    # Keep only one annotation per scan
    df_agg = df_agg.group_by(['scan_number']).agg([
        pl.col('lipid').str.concat("|"),
        pl.col('adduct').str.concat("|"),
        pl.first('rt_s'),
        pl.mean('precursor_mz'),
        pl.first('Fraction_essential_detected'),
        pl.first('n_qualifying_fragments_detected'),
        pl.first('n_qualifying_fragments'),
        pl.col('Fragments_detected').str.concat("|"),
        pl.col('lipid-adduct').str.concat("|")])
    # Add a "representative" lipid
    df_agg = df_agg.with_columns(pl.col("lipid")
                .str.split("|")
                .list.get(0).alias("Representative_lipid"))

    return df_agg

#%% The following is the equivalent of main(), just not incorporated into one. You should be able to just run the code section (ctrl+enter)

# Merge LipidMatch libraries into one dataframe and create a rules dataframe for linking theoretical lipids to detected MS2 spectra.

df_alllipids = merge_lipid_libraries(files_sel,dir_library) # this produces a large dataframe with all lipids from LipidMatch
df_rules = define_rules(os.path.join(dir_library, "LIPID_ID_CRITERIA.xlsx")) # Creates a dataframe for judging whether an MS2 spectrum could represent a given lipid. The rules are class-dependent.


#%%  This snippet will read in each DDA file, and match it to the theoretical lipids
# then it'll figure out which fragments correspond to a known lipid fragment
# and only keep fragments/data that correspond to fragments that occur in lipids

dda_file_names = [x for x in os.listdir(dir_data) if x[0]!='.']
df_foundlipids = match_dda_to_alllipids(dda_file_names, dir_data, df_alllipids, mztol=mztol)
df_foundlipids = prepare_lipid_data_with_rules(df_foundlipids)
# Add the number of total qualifying and essential fragments
# (then later we'll set thresholds to what percentage of each fragment type must be present)
df_foundlipids = df_foundlipids.with_columns(n_essential_fragments = pl.col("All_essential_fragments").list.len())
df_foundlipids = df_foundlipids.with_columns(n_qualifying_fragments = pl.col("Qualifying_fragments").list.len())
df_foundlipids = aggregate_and_filter_lipids_by_rules(df_foundlipids)


#%% Write the lipid library into an excel file as well as mzmine compatible csv file (with less information)
dir_output = os.path.join(dir_base, "output")


try:  
    os.mkdir(dir_output)  
except OSError as error:  
    print(error)  


df_foundlipids.write_excel(dir_output+r'/rulebased_lipids.xlsx')

df_rulebased_mzmine = df_foundlipids.rename({'lipid':'name','rt_s':'rt','precursor_mz':'mz'})
df_rulebased_mzmine = df_rulebased_mzmine[['name','rt','mz','adduct']]
df_rulebased_mzmine = df_rulebased_mzmine.with_columns(pl.lit(df_rulebased_mzmine['rt']/60).alias('rt'))
df_rulebased_mzmine.write_csv(dir_output+r'/rulebased_lipids_mzmine.xlsx')
