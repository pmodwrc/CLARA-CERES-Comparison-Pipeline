###############################################################################
""" Import the necessary Packages """
import os
import gc
import pandas as pd

from datetime import datetime as dt
from dateutil.relativedelta import relativedelta

import yaml
import csv
from scipy.io import readsav
import xarray as xr

# Load functions from other files
from . import math_and_conversion as mc
###############################################################################
#%%
###############################################################################

""" File Metadata & Naming """

###############################################################################
#%%
###############################################################################
def extract_file_dates(directory, file_extension):
    """
    Extracts start and end dates from filenames in a given directory.

    This function scans a directory for files matching the specified prefix 
    and extension. It extracts timestamps from filenames, where the timestamp 
    is either split by a separator or just a single timestamp which is 
    formatted according to the provided `date_format`. 
    The extracted data is returned as a DataFrame.

    Parameters:
        directory (dict): A dict contsining the file path information
        - folder (str): Path to the directory containing the files.
        - prefix (str): Prefix that identifies relevant files.
        - suffix (str, optional): Suffix that identifies rlevant files
        - date_format (str): Format in which the date is represented 
                             in the filename (e.g. '%Y%m%d%H').
        - date_separator (str, optional): If the date range is given, this is the
                                          separator between the start and end date.
                                          Otherwise the date range is assumed to 
                                          be the full Year/Month/Day 
        file_extension (str): File extension (e.g., '.nc', '.pkl').
        
    Returns:
        pd.DataFrame: A DataFrame with columns:
            - 'start_time' (datetime): Start timestamp parsed from filename.
            - 'end_time' (datetime): End timestamp parsed from filename.
            - 'filename' (str): Original filename.
    """
    files = {'start_time': [], 'end_time': [], 'filename': []}
    for fname in os.listdir(directory['folder']):
        name, ext = os.path.splitext(fname)
        if ext.lower() != f".{file_extension.lower()}":
            continue
        if directory['prefix'] and not name.startswith(directory['prefix']):
            continue
        if directory['suffix'] and not name.endswith(directory['suffix']):
            continue
        files['filename'].append(fname)
        
        middle = name
        if directory['prefix']:
            middle = middle.replace(directory['prefix'],'')
        if directory['suffix']:
            middle = middle.replace(directory['suffix'],'')
        
        if not directory['date_separator']:
            date_range = relativedelta(days=1)
            if '%d' not in directory['date_format']:
                date_range = relativedelta(months=1)
            if '%m' not in directory['date_format']:
                date_range = relativedelta(years=1)
            start_dt = dt.strptime(middle, directory['date_format'])
            end_dt = start_dt+date_range
        else:
            start, end = middle.split(directory['date_separator'])
            start_dt = dt.strptime(start, directory['date_format'])
            end_dt = dt.strptime(end, directory['date_format'])
        files['start_time'].append(start_dt)
        files['end_time'].append(end_dt)
    return pd.DataFrame(files).sort_values('start_time').reset_index(drop=True)

###############################################################################
#%%
###############################################################################
def filter_files_by_date_range(files_df, start, end):
    """
    Filters a DataFrame of files for those with a date range overlap.
    """
    return files_df[(files_df['start_time'] <= end) & (files_df['end_time'] >= start)]

###############################################################################
#%%
###############################################################################
def create_filename(file_dict, start_time, end_time=''):
    """
    Creates a filename based on the specified file structure.
    """
    date_string = start_time.strftime(file_dict['date_format'])
    if file_dict['date_separator']:
        date_string += file_dict['date_separator']+end_time.strftime(file_dict['date_format'])
    filename = file_dict['prefix']+date_string
    if file_dict['suffix']:
        filename += file_dict['suffix']
    return filename

###############################################################################
#%%
###############################################################################

""" Loading Data """

###############################################################################
#%%
###############################################################################
def load_config(config_path="config.yaml"):
    """
    Loads the configuration from a YAML file and returns the relevant information.
    """
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    module = config['module']
    output = config['output']
    directory = config['directory']
    parameters = config['parameters']
    input_file_structure = config['input_file_structure']
    return module, output, directory, parameters, input_file_structure
###############################################################################
#%%
###############################################################################
def read_norsat_position(time, filename, input_file_structure):
    """
    Reads a CLARA satellite position and orientation CSV file and extracts relevant data.

    Parameters:
        time (datetime.datetime): The start time for filtering the data.
        filename (str): Path to the CSV file.
        input_file_structure (dict): A dict on how the file is structured with
            - 'time_format'
            - 'column_names'

    Returns:
        pd.DataFrame: DataFrame containing the extracted and filtered data.
    """
    # The key of the data columns used in input_file_structure
    # and now for the extracted data
    data_columns=["qx", "qy", "qz", "qw", "sx", "sy", "sz"]
    
    
    # Define the time range for extraction (entire day of the given time)
    start_time = time
    end_time = start_time.replace(hour=23, minute=59, second=59)

    extracted_data = {col: [] for col in data_columns}
    extracted_data['CLARA_time_utc'] = []

    with open(filename, newline='') as csvfile:
        csv_reader = csv.DictReader(csvfile, delimiter=",", lineterminator="\n")

        for row in csv_reader:
            # Parse timestamp from the CSV
            row_time = dt.strptime(row[input_file_structure['column_names']['time']], input_file_structure['time_format'])
            
            # Skip rows outside the time range
            if row_time < start_time:
                continue
            if row_time > end_time:
                break

            # Validate and extract data
            row_valid = all(mc.safe_float_conversion(row[input_file_structure['column_names'][col]]) != "No Float" for col in data_columns)

            if row_valid and row[input_file_structure['column_names']['time']] != 'nan':
                extracted_data['CLARA_time_utc'].append(row_time)
                for col in data_columns:
                    extracted_data[col].append(float(row[input_file_structure['column_names'][col]]))

    return pd.DataFrame(extracted_data)
###############################################################################
#%%
###############################################################################
def get_clara_olr(OLR_files_df,OLR_file_path,fov_angle, input_file_structure):
    """
    Reads and processes CLARA Outgoing Longwave Radiation (OLR) data from multiple .save files.

    Parameters:
    ----------
    OLR_files_df : pd.DataFrame
        A DataFrame containing a column 'filename' with the list of OLR .sav files to be processed.
    OLR_file_path : str
        Path to the directory containing the OLR .save files.
    fov_angle : float
        The opening angle of the satellite's field of view, used for radiance conversion.
    input_file_structure : dict
        A dictionary defining the structure of the .sav files, including:
        - 'key': The main key in the .save file that contains data.
        - 'column_names': A dictionary mapping required fields to their corresponding keys in the .save file.

    Returns:
    -------
    CLARA_OLR : pd.DataFrame
        A DataFrame containing the processed OLR data with the following columns:
        - 'OLR_index': Index of each measurement in the dataset.
        - 'JDAY': Julian day timestamp of the measurement.
        - 'CLARA_radiance': Converted radiance values.
        - 'CLARA_OLR_time': Timestamp in datetime format.
        - 'sav filename': The name of the source .sav file for each row.
    """
    
    OLR_dfs = []  # List to store individual DataFrames from each file

    # Iterate over each OLR file in the DataFrame
    for filename in OLR_files_df['filename']:
        # Load the .save file using IDL's readsav function
        OLR_sav = readsav(os.path.join(OLR_file_path, filename))

        # Extract time and radiance values from the file using input_file_structure
        OLR_df = pd.DataFrame({
            'JDAY': OLR_sav[input_file_structure['key']][input_file_structure['column_names']['time']][0].tolist(),
            'CLARA_radiance': OLR_sav[input_file_structure['key']][input_file_structure['column_names']['radiance']][0].tolist()
        })

        # Convert Julian day timestamps to standard datetime format
        OLR_df['CLARA_OLR_time'] = pd.to_datetime(OLR_df['JDAY'], origin='julian', unit='D')

        # Reset index, rename it to 'OLR_index', and sort by time
        OLR_df = OLR_df.reset_index().rename(columns={'index': 'OLR_index'}).sort_values('CLARA_OLR_time')

        # Store the filename in the DataFrame for traceability
        OLR_df['sav filename'] = filename

        # Append processed DataFrame to the list
        OLR_dfs.append(OLR_df)
    # Concatenate all individual DataFrames into one and sort by time
    CLARA_OLR = pd.concat(OLR_dfs).reset_index(drop=True).sort_values('CLARA_OLR_time')

    # Convert radiance values using the specified field-of-view angle
    CLARA_OLR['CLARA_radiance'] = mc.clara_radiance_conversion(fov_angle)*CLARA_OLR['CLARA_radiance']
    return CLARA_OLR
###############################################################################
#%%
###############################################################################
def load_and_merge_ceres(CERES_files_df, CERES_file_path, input_file_structure):
    """
    Optimized function to load and merge CERES .nc files efficiently.
    """

    # Convert filenames to full paths
    file_paths = [os.path.join(CERES_file_path, file) for file in CERES_files_df['filename']]

    # Use open_mfdataset for efficient loading (lazy loading enabled)
    time_ds = xr.open_mfdataset(file_paths, combine="nested", concat_dim="time", chunks={})
    position_ds = xr.open_mfdataset(file_paths, group=f'/{input_file_structure["group_names"]["position"]}', decode_times=False, combine="nested", concat_dim="time", chunks={})
    radiance_ds = xr.open_mfdataset(file_paths, group=f'/{input_file_structure["group_names"]["radiance"]}', combine="nested", concat_dim="time", chunks={})
    angles_ds = xr.open_mfdataset(file_paths, group=f'/{input_file_structure["group_names"]["angles"]}', combine="nested", concat_dim="time", chunks={})
    surface_ds = xr.open_mfdataset(file_paths, group=f'/{input_file_structure["group_names"]["surface"]}', combine="nested", concat_dim="time", chunks={})

    # Convert to DataFrames
    time_df = time_ds[input_file_structure['variable_names']['time']].to_dataframe()
    position_df = position_ds.to_dataframe()
    radiance_df = radiance_ds.to_dataframe()
    angles_df = angles_ds.to_dataframe()
    surface_df = surface_ds.to_dataframe()

    # Create time mapping
    time_map = (time_df.reset_index(drop=True)
                      .rename(columns={input_file_structure['variable_names']['time']: 'julian_observation_time'})
                      .rename_axis(input_file_structure['variable_names']['time']))

    # Map indices to time mapping
    position_df.index = position_df.index.map(time_map['julian_observation_time'])
    radiance_df.index = radiance_df.index.map(time_map['julian_observation_time'])
    angles_df.index = angles_df.index.map(time_map['julian_observation_time'])
    surface_df.index = surface_df.index.set_levels(
        surface_df.index.levels[0].map(time_map['julian_observation_time']), 
        level=input_file_structure['variable_names']['time']
    )

    return position_df, radiance_df, angles_df, surface_df
###############################################################################
#%%
###############################################################################
