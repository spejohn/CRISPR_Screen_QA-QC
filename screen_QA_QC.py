
#Import packages

#pip install gspread oauth2client

import argparse
import os
from pathlib import Path, PurePosixPath


from gooey import Gooey, GooeyParser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from plotly import graph_objects as go
from scipy.stats import pearsonr
from sklearn.metrics import auc, precision_recall_curve, roc_curve
from sklearn.linear_model import LinearRegression


#Access a Google Sheet tab and search for "First" and "Last" columns within the Sheet "Test"

def calc_stats(df1: pd.DataFrame, df2: pd.DataFrame, R_list: list, slope_list: list):
    # Determine the column to match rows by
    match_column = 'GENE' if 'DZ' in df1 or 'DZ' in df2 else 'id' if 'MGK' in df1 or 'MGK' in df2 else None

    if match_column is not None:
        # Match rows by the determined column
        df1 = df1.set_index(match_column)
        df2 = df2.set_index(match_column)
        df1, df2 = df1.align(df2, join='inner')
        df1 = df1.reset_index()
        df2 = df2.reset_index()

    data_columns = ['normZ', 'fdr_synth', 'fdr_supp', 'rank_synth', 'neg|lfc', 'neg|score', 'pos|score', 'neg|fdr', 'pos|fdr', 'neg|rank']
    for column in data_columns:
        if column in df1.columns and column in df2.columns:
            X = df1[column].values.reshape(-1, 1)
            y = df2[column].values
            reg = LinearRegression().fit(X, y)
            r_squared = reg.score(X, y)
            slope = reg.coef_[0]
        
            if column in R_list:
                        R_list[column] = [r_squared]  # Replace the list with a new list containing one value
                        slope_list[column] = [slope]  # Replace the list with a new list containing one value
        else:
            print(f"Column '{column}' does not exist in one or both dataframes.")

def create_fill_QC_df(authors_file_path:str, data_folder: str, output_dir: str):
    df = pd.read_csv(authors_file_path)
    df.columns = ['total_authors']
    # Remove rows that begin or end with "_"
    df = df[~df['total_authors'].str.startswith("_")]
    df = df[~df['total_authors'].str.endswith("_")]
    
    unique_authors = df['total_authors'].unique()
    authors_df = pd.DataFrame(unique_authors, columns = ['author_First_Last'])

    
    QC_folders = []  # List to store data before creating DataFrame
    
    for author in authors_df['author_First_Last']:
        author_folder = Path(data_folder) / str(author)
        #print(f"Searching in folder: {author_folder}")
        if author_folder.exists() and author_folder.is_dir():
            for folder in author_folder.iterdir():
                if folder.is_dir():
                    print(f"Found folder: {folder.name}")
                    QC_folders.append({"First_Last": str(author), "Screen_ID": str(folder.name)})

    column_names=["First_Last", "Screen_ID", 
                        "MGK_R2_neg|lfc", 
                         "MGK_R2_neg|score", 
                         "MGK_R2_pos|score", 
                         "MGK_R2_neg|fdr", 
                         "MGK_R2_pos|fdr", 
                         "MGK_R2_neg|rank", 
                         "DZ_R2_normZ", 
                         "DZ_R2_fdr_synth",
                         "DZ_R2_fdr_supp", 
                         "DZ_R2_rank_synth", 
                         "MGK_slope_neg|lfc", 
                         "MGK_slope_neg|score", 
                         "MGK_slope_neg|fdr",
                         "MGK_slope_pos|fdr",
                         "MGK_slope_neg|rank",
                         "DZ_slope_normZ",
                         "DZ_slope_fdr_synth",
                         "DZ_slope_fdr_supp",
                         "DZ_slope_rank_synth"]
    
    # Convert list of dictionaries to DataFrame
    QC_df = pd.DataFrame(QC_folders, columns=column_names)

    # Loop through each row of QC_df
    for index, row in QC_df.iterrows():
        # Construct directory path for each row
        dir_path = str(data_folder + f"\\" + row['First_Last'] + f"\\" + row['Screen_ID'] + f"\\" + row['Screen_ID'])
        print(f"Searching in directory: {dir_path}")
        
        # Define file suffixes for different types of files
        suffixes = ['_pDZ', '_gDZ', '_pMGK', '_gMGK']
        
        # Initialize variables to store statistics
        MGK_R2 = {'neg|lfc': [], 'neg|score': [], 'pos|score': [], 'neg|fdr': [], 'pos|fdr': [], 'neg|rank': []}
        DZ_R2 = {'normZ': [], 'fdr_synth': [], 'fdr_supp':[], 'rank_synth': []}
        MGK_slope = {'neg|lfc': [], 'neg|score': [], 'pos|score':[], 'neg|fdr': [], 'pos|fdr': [], 'neg|rank': []}
        DZ_slope = {'normZ': [], 'fdr_synth': [], 'fdr_supp': [], 'rank_synth': []}

        # Reset the dataframes for each iteration
        pDZ_df = None
        gDZ_df = None
        pMGK_df = None
        gMGK_df = None

        # Loop through the suffixes
        for suffix in suffixes:
            # Create the directory path based on the suffix
            file_str = dir_path + suffix
            file_path = Path(file_str).with_suffix(".csv")
                
            # Check if the file exists
            if file_path.exists() and file_path.is_file():
                print(f"Created path to {file_path}.")
                # Read the file into a DataFrame
                if f'_pDZ' in file_path.name:
                    pDZ_df = pd.read_csv(file_path)
                    print(f"Reading {row['Screen_ID']}{suffix} file into dataframe.")
                if f'_gDZ' in file_path.name:
                    gDZ_df = pd.read_csv(file_path)
                    print(f"Reading {row['Screen_ID']}{suffix} file into dataframe.")
                if f'_pMGK' in file_path.name:
                    pMGK_df = pd.read_csv(file_path)
                    print(f"Reading {row['Screen_ID']}{suffix} file into dataframe.")
                if f'_gMGK' in file_path.name:
                    gMGK_df = pd.read_csv(file_path)
                    print(f"Reading {row['Screen_ID']}{suffix} file into dataframe.")
                
            if pDZ_df is not None and gDZ_df is not None:
                #Use calc_stats() with pDZ_df and gDZ_df with DZ_R2 and DZ_slope lists.
                calc_stats(pDZ_df, gDZ_df, DZ_R2, DZ_slope)
                print(f"Calculated stats for pDZ and gDZ tables.")
            if pMGK_df is not None and gMGK_df is not None:
                calc_stats(pMGK_df, gMGK_df, MGK_R2, MGK_slope)
                print(f"Calculated stats for pMGK and gMGK tables.")

        # Fill the statistics into the corresponding columns in QC_df
        print(f"Filling calculated stats from {dir_path}{suffix}")
        if pDZ_df is not None and gDZ_df is not None:
            for key, value in DZ_R2.items():
                QC_df.at[index, f"DZ_R2_{key}"] = value
            for key, value in DZ_slope.items():
                QC_df.at[index, f"DZ_slope_{key}"] = value
        if pMGK_df is not None and gMGK_df is not None:
            for key, value in MGK_R2.items():
                QC_df.at[index, f"MGK_R2_{key}"] = value
            for key, value in MGK_slope.items():
                QC_df.at[index, f"MGK_slope_{key}"] = value
        else:
            print(f"No values were calculated from dataframe pairs.")


    output_path = Path(output_dir) / "final_QC_df_output.csv"
    QC_df.to_csv(output_path)
    return QC_df

#author_path = "C:\\Users\\spejo\\Documents\\BI_Pipeline\\QC_tests\\authors.csv"
#data_path = "C:\\Users\\spejo\\Documents\\BI_Pipeline\\QC_tests"
#output_folder = "C:\\Users\\spejo\\Documents\\BI_Pipeline\\QC_tests\\output_practice"
#create_fill_QC_df(author_path, data_path, output_folder)

'''
#Access a Google Sheet tab and search for "First" and "Last" columns within the Sheet "Test"
import gspread
from oauth2client.service_account import ServiceAccountCredentials

def find_columns(sheet_url):
    # Use the service account credentials to authenticate
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    credentials = ServiceAccountCredentials.from_json_keyfile_name('path_to_your_credentials_file.json', scope)
    gc = gspread.authorize(credentials)

    # Open the Google Spreadsheet using its URL
    spreadsheet = gc.open_by_url(sheet_url)

    # Select the "Test" worksheet from the spreadsheet
    worksheet = spreadsheet.worksheet('Test')

    # Get all the values in the worksheet
    values = worksheet.get_all_values()

    # Get the first row (which should contain the column names)
    column_names = values[0]

    # Check if "First" and "Last" are in the column names
    if 'First' in column_names and 'Last' in column_names:
        print('Columns "First" and "Last" found.')
    else:
        print('Columns "First" and "Last" not found.')

# Call the function with the URL of your Google Sheet
find_columns('your_google_sheet_url')

#TODO: create a conditional that removes if beginning or ending with "_"

#TODO add arguments to run all functions
def run_pipeline(args):
    create_authors_df(args.authors_file_path)
    create_QC_df()
'''

@Gooey(optional_cols=1)
def main():
    parser = GooeyParser(description="Run QA/QC on provided & generated DrugZ/Mageck files.")
    #input_group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument(
        'authors_file_path', 
        widget='FileChooser', 
        help="Select path to the authors.csv file"
    )
    parser.add_argument(
        "data_folder",
        type=str,
        help="Select path to Data folder",
        widget='DirChooser'
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Select path to directory for storing QA/QC file. If not specified, an 'output' folder in the current directory will be created.",
        default=".",
        widget="DirChooser"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="If set, overwrite existing output files.",
        default=False,
    )
    args = parser.parse_args()
    create_fill_QC_df(args.authors_file_path, args.data_folder, args.output_dir)

if __name__ == "__main__":
    main()
