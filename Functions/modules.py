###############################################################################
""" Import the necessary Packages """
import os
import gc
import pandas as pd
import geopandas as gpd
import dask_geopandas as dgp

import calendar
from datetime import datetime as dt

from pysolar.solar import get_altitude

import pickle

from tqdm import tqdm

# Load functions from other files
from . import file_handling as fh
from . import coordinates as coord
from . import matching_and_merging as mm
from . import output as out
###############################################################################
#%%
###############################################################################
def create_coordinates_file(config_path):
    """
    Creates monthly coordinate files for the CLARA satellite from the provided CSV files.

    This function processes satellite position and orientation data, computes various coordinate values,
    and stores them in pickle files grouped by month.

    Parameters:
    ----------
    config_path (str): path to the config file which includes:
        - module: Determine if this module is active or not
        - directory['norsat_pos']: Input file handling
        - input_file_structure['norsat_pos']: Input file structure
        - directory['clara_pos']: Output file handling
        - output['save_intermediate']['CLARA_position_csv']: If the output 
            is also saved as a .csv file
        - parameters: Parameters for get_coordinates
            - fov_angle: Field of view angle for satellite in degrees
            - num_footprint_vertices: Number of points to sample on the fov cone
            - max_view_zenith_angle: Limiting View Zenith Angle
    Returns:
    -------
    None: The function saves processed data to pickle files in the specified output directory.
    
    Functions:
    ---------
    file_handling.py:
    - load_config           - create_filename
    - read_norsat_position
    coordinates.py:
    - get_coordinates       -calculate_local_time_after_midnight
    """

    ###########################################################################
    """ 1: Load in the configuration settings """
    module, output, directory, parameters, input_file_structure = fh.load_config(config_path)

    # Ensure the module is enabled for combination processing
    if module not in ('Both','Coordinates'):
        return

    """ 2: Directory """
    ### a) Check the Directory and Create Folders if Necessary ###
    for key in directory:
        if key in ['norsat_pos']:
            if not os.path.exists(directory[key]['folder']):
                print(f"{directory[key]['folder']} does not exist.")
                return
        elif key in ['clara_pos']:
            os.makedirs(directory[key]['folder'], exist_ok=True)
            if output['save_intermediate']['CLARA_position_csv']:
                os.makedirs(directory[key]['folder']+'/csv', exist_ok=True)

    
    ### b) Create a DataFrame containing input file information ###
    files = {'filename': [], 'date': []}
    for fname in os.listdir(directory['norsat_pos']['folder']):
        name, ext = os.path.splitext(fname)
        if ext.lower() != '.csv':
            continue
        if directory['norsat_pos']['prefix'] and not name.startswith(directory['norsat_pos']['prefix']):
            continue
        if directory['norsat_pos']['suffix'] and not name.startswith(directory['norsat_pos']['suffix']):
            continue
        files['filename'].append(fname)

        date = name
        if directory['norsat_pos']['prefix']:
            date = date.replace(directory['norsat_pos']['prefix'],'')
        if directory['norsat_pos']['suffix']:
            date = date.replace(directory['norsat_pos']['suffix'],'')
        date = dt.strptime(date,directory['norsat_pos']['date_format'])
        files['date'].append(date)
    if not files['filename']:
        print("No files found matching the given criteria.")
        return
    files_df = pd.DataFrame(files)
    
    files_df['Year']  = files_df['date'].dt.year
    files_df['Month']  = files_df['date'].dt.month
    files_df['Day']  = files_df['date'].dt.day
    
    files_df = files_df.set_index(['Year', 'Month', 'Day']).sort_index()

    
    """ 3: Process NorSat files and group them into months """
    for year in files_df.index.get_level_values(0).unique():
        for month in files_df.loc[year].index.get_level_values(0).unique():
            # Initialize a list to store data for each day in the month
            df_list = []
            for day in tqdm(files_df.loc[(year, month)].index.get_level_values(0).unique(),
                            total=files_df.loc[(year, month)].shape[0],
                            desc=f"Processing {calendar.month_name[month]} {year}"):
                # Get the filename and information for the Day
                filename = files_df.loc[(year, month, day)]['filename']
                full_filename = os.path.join(directory['norsat_pos']['folder'], filename)
                time_dt = files_df.loc[(year, month, day)]['date']

                ### a)  Extract the data from the CSV file ###
                extracted_data = fh.read_norsat_position(time_dt, full_filename, input_file_structure['norsat_pos'])
                
                ### b) Calculate coordinates ###
                if extracted_data.empty:
                    coordinates = pd.DataFrame()
                else:
                    coordinates = coord.get_coordinates(
                        extracted_data, parameters['fov_angle'], 
                        parameters['num_footprint_vertices'], 
                        parameters['max_view_zenith_angle']
                    )

                if not coordinates.empty:
                    df_list.append(coordinates)

            """ 4: Combine daily DataFrames into one monthly DataFrame """
            if df_list:
                coordinates_month = pd.concat(df_list, ignore_index=True)

                ### a) Add local time and solar zenith angle to the dataframe ###
                coordinates_month['CLARA_local_time'] = coordinates_month.apply(
                    lambda row: coord.calculate_local_time_after_midnight(
                        row['CLARA_fov_longitude'], 
                        row['CLARA_fov_latitude'], 
                        row['CLARA_time_utc']
                    ), axis=1)
                
                coordinates_month['CLARA_solar_zenith_angle'] = coordinates_month.apply(
                    lambda row: float(90) - get_altitude(row['CLARA_fov_latitude'], row['CLARA_fov_longitude'], row['CLARA_time_utc'].tz_localize('UTC')), axis=1)

                # Rearrange columns in the desired order
                coordinates_month = coordinates_month[[
                    'CLARA_time_utc', 'CLARA_satellite_eci_coordinates',
                    'CLARA_satellite_pointing_vector_eci', 'CLARA_view_zenith_angle',
                    'CLARA_fov_longitude', 'CLARA_fov_latitude',
                    'CLARA_fov_longitude_positive', 'CLARA_fov_colatitude',
                    'CLARA_fov_footprint', 'CLARA_fov_footprint_match', 
                    'CLARA_subsatellite_longitude', 'CLARA_subsatellite_latitude', 
                    'CLARA_subsatellite_longitude_positive', 
                    'CLARA_subsatellite_colatitude', 'CLARA_local_time',
                    'CLARA_solar_zenith_angle'
                ]]

                ### b) Save the Data ###
                # Generate the filename
                save_filename = fh.create_filename(
                    directory['clara_pos'], 
                    min(coordinates_month['CLARA_time_utc']).to_pydatetime(), 
                    max(coordinates_month['CLARA_time_utc']).to_pydatetime()
                )
                
                # Save .pkl file
                save_file_path_pkl = os.path.join(directory['clara_pos']['folder'], 
                                                  save_filename+'.pkl')
                coordinates_month.to_pickle(save_file_path_pkl)
                print(f"Data saved to {save_file_path_pkl}")
                
                # save .csv file 
                if output['save_intermediate']['CLARA_position_csv']:
                    save_file_path_csv = os.path.join(directory['clara_pos']['folder']+'/csv', 
                                                      save_filename+'.csv')
                    coordinates_month.to_csv(save_file_path_csv)
                    print("--------------------------------------------------")
                    print(f"Data saved to {save_file_path_csv}")
    
    return
###############################################################################
#%%
###############################################################################
def combine_CLARA_CERES(config_path):
    """
    Combines CLARA and CERES satellite datasets by matching spatial and temporal data.

    Parameters:
    ----------
    config_path : str, optional
        Path to the configuration file (default is "config.yaml").

    Returns:
    -------
    None
        The function processes and saves the combined CLARA-CERES data in various formats.

    Functions:
    ----------
    file_handling.py:
    - load_config                  - create_filename, 
    - extract_file_dates           - filter_files_by_date_range
    - get_clara_olr                - load_and_merge_ceres
    matching_and_merging.py:
    - match_ceres_clara_positions  - filter_ceres_clara_by_time
    - add_ceres_data_to_clara      - get_clara_surface
    output.py:
    - save_netCDF
    
    """

    """ 1: Load in the configuration settings """
    module, output, directory, parameters, input_file_structure = fh.load_config(config_path)

    # Ensure the module is enabled for combination processing
    if module not in ('Both','Combination'):
        return

    """ 2: Directory """
    ### a) Check the Directory and Create Folders if Necessary ###
    for key in directory:
        # Verify that required input directories exist
        if key in ['clara_olr', 'ceres','clara_pos']:
            if not os.path.exists(directory[key]['folder']):
                print(f"{directory[key]['folder']} does not exist.")
                return
        # Ensure directories for combined output exist
        elif key in ['combined']:
            os.makedirs(directory[key]['folder'], exist_ok=True)
            # Create separate folders for saving intermediate files
            if output['save_intermediate'][key] in '.pkl':
                os.makedirs(directory[key]['folder']+'/pkl', exist_ok=True)
            elif output['save_intermediate'][key] in '.csv':
                os.makedirs(directory[key]['folder']+'/csv', exist_ok=True)
            elif output['save_intermediate'][key] in 'both':
                os.makedirs(directory[key]['folder']+'/pkl', exist_ok=True)
                os.makedirs(directory[key]['folder']+'/csv', exist_ok=True)
        # Ensure directories for the combination key exist if necessary
        elif key in ['combination_key']:
            # Create separate folders for saving intermediate files
            if output['save_intermediate'][key] in '.pkl':
                os.makedirs(directory[key]['folder'], exist_ok=True)
            elif output['save_intermediate'][key] in '.csv':
                os.makedirs(directory[key]['folder'], exist_ok=True)
            elif output['save_intermediate'][key] in 'both':
                os.makedirs(directory[key]['folder']+'/pkl', exist_ok=True)
                os.makedirs(directory[key]['folder']+'/csv', exist_ok=True)
        elif key in ['spatial_match']:
            if output['save_intermediate'][key]:
                os.makedirs(directory[key]['folder'], exist_ok=True)
                
    ### b) Create DataFrames containing file information ###
    CLARA_position_files = fh.extract_file_dates(directory['clara_pos'],'pkl')
    CLARA_OLR_files = fh.extract_file_dates(directory['clara_olr'],'save')
    CERES_files = fh.extract_file_dates(directory['ceres'],'nc')

    ### c) Check Time ranges are matching ###
    if not (
        (CLARA_OLR_files['start_time'].min() 
             <= CLARA_position_files['start_time'].min()) 
        and 
        (CLARA_OLR_files['end_time'].max() 
             >= CLARA_position_files['end_time'].max())
    ):
        print('WARNING: There are no CLARA OLR files for some date ranges')
        return
    if not (
        (CERES_files['start_time'].min() 
             <= CLARA_position_files['start_time'].min()) 
        and 
        (CERES_files['end_time'].max() 
             >= CLARA_position_files['end_time'].max())
    ):
        print('WARNING: There are no CERES files for some date ranges')
        return

    """ 3: Process each CLARA position file """
    for i, CLARA_file in CLARA_position_files.iterrows():
        function_timing = [dt.now()]           # Track function execution times
        ### a) Load CLARA Position Data ###
        CLARA_filename = CLARA_file['filename']

        with open(os.path.join(directory['clara_pos']['folder'], CLARA_filename), 'rb') as f:
            CLARA_position = pickle.load(f)
        
        print(f"({i+1}/{len(CLARA_position_files)}) Working on {CLARA_filename} - ({function_timing[0]})")

        ### b) Find matching CLARA OLR and CERES data ###
        date_filter_start = CLARA_file['start_time'].to_pydatetime()
        date_filter_end = CLARA_file['end_time'].to_pydatetime()
        
        CLARA_OLR_files_filtered = fh.filter_files_by_date_range(CLARA_OLR_files, 
                                                                 date_filter_start, 
                                                                 date_filter_end)
        CERES_files_filtered = fh.filter_files_by_date_range(CERES_files, 
                                                             date_filter_start, 
                                                             date_filter_end)
        
        print(f"  Working with {len(CERES_files_filtered)} CERES files")

        ### c) Load the CLARA OLR data ###
        CLARA_OLR = fh.get_clara_olr(CLARA_OLR_files_filtered, 
                                     directory['clara_olr']['folder'], 
                                     parameters['fov_angle'], 
                                     input_file_structure['clara_olr'])
        
        """ 4: Load the CERES data and match it with the CLARA data """
        ### a) Load the CERES data ###
        CERES_position, CERES_radiance, CERES_angles, CERES_surface = fh.load_and_merge_ceres(
            CERES_files_filtered, directory['ceres']['folder'], 
            input_file_structure['ceres']
        )
        function_timing.append(dt.now())
        print(f"    (1/6) Load CERES Data: {function_timing[-1]-function_timing[-2]}")
        
        ### b) Spatially combine CERES and CLARA data ###
        spatially_matched_df = mm.match_ceres_clara_positions(
            CERES_position, CLARA_position, input_file_structure['ceres']
        )
        del CERES_position
        
        # Save in .pkl or .csv format
        if output['save_intermediate']['spatial_match']:
            spatial_match_filename = fh.create_filename(directory['spatial_match'], 
                                                        date_filter_start, 
                                                        date_filter_end)
            spatial_match_full_filename_gdf = directory['spatial_match']['folder']+'/'+spatial_match_filename+'.parquet'
            
            spatially_matched_df.to_parquet(spatial_match_full_filename_gdf, engine='pyarrow', compression='snappy')
        
        # Print Step Completion
        function_timing.append(dt.now())
        print(f"    (2/6) Spatially Match CLARA and CERES: {function_timing[-1]-function_timing[-2]}")
        if output['save_intermediate']['spatial_match']:
            print('    ------------------------------------------------------')
            print(f"      Saved to {spatial_match_full_filename_gdf}")
        
        ### c) Filter the CERES&CLARA combination based on time ###
        combination_key = mm.filter_ceres_clara_by_time(
            spatially_matched_df, parameters['time_bin_size']
        )
        del spatially_matched_df
        # Save in .pkl or .csv format
        key_filename = fh.create_filename(directory['combination_key'], date_filter_start, date_filter_end)
        
        if output['save_intermediate']['combination_key'] in '.pkl':
            key_full_filename_pkl = directory['combination_key']['folder']+'/'+key_filename+'.pkl'
            combination_key.to_pickle(key_full_filename_pkl)
        if output['save_intermediate']['combination_key'] in '.csv':
            key_full_filename_csv = directory['combination_key']['folder']+'/'+key_filename+'.csv'
            combination_key.to_csv(key_full_filename_csv)
        if output['save_intermediate']['combination_key'] in 'both':
            key_full_filename_pkl = directory['combination_key']['folder']+'/pkl/'+key_filename+'.pkl'
            key_full_filename_csv = directory['combination_key']['folder']+'/csv/'+key_filename+'.csv'
            combination_key.to_pickle(key_full_filename_pkl)
            combination_key.to_csv(key_full_filename_csv)
        
        # Print Step Completion
        function_timing.append(dt.now())
        print(f"    (3/6) Temporal Match: {function_timing[-1]-function_timing[-2]}")
        if output['save_intermediate']['combination_key'] in '.pkl':
            print('    ------------------------------------------------------')
            print(f"      Saved to {key_full_filename_pkl}")
        if output['save_intermediate']['combination_key'] in '.csv':
            print('    ------------------------------------------------------')
            print(f"      Saved to {key_full_filename_csv}")
        if output['save_intermediate']['combination_key'] in 'both':
            print('    ------------------------------------------------------')
            print(f"      Saved to {key_full_filename_pkl}")
            print(f"      Saved to {key_full_filename_csv}")
        
        ### d) Combine CLARA and CERES data using the combination key ###
        CLARA_CERES_combined = mm.add_ceres_data_to_clara(
            CLARA_position, CLARA_OLR, combination_key, CERES_radiance, 
            CERES_angles, input_file_structure['ceres']
        )
        del CLARA_position, CLARA_OLR, CERES_radiance, CERES_angles
        # Print Step Completion
        function_timing.append(dt.now())
        print(f"    (4/6) Combine CLARA with CERES: {function_timing[-1]-function_timing[-2]}")
        
        ### e) Create the CLARA_surface DataFrame using the combination key ###
        CLARA_surface = mm.get_clara_surface(combination_key,CERES_surface, input_file_structure['ceres'])
        del CERES_surface
        # Print Step Completion
        function_timing.append(dt.now())
        print(f"    (5/6) CLARA_surface: {function_timing[-1]-function_timing[-2]}")
        
        """ 5: Save the Data """
        ### a) Apply time difference filters ###
        if output['filters']['pos_olr_time_diff']:
            filter_olr = CLARA_CERES_combined['CLARA_position_olr_time_diff'] < pd.Timedelta(seconds=output['filters']['pos_olr_time_diff'])
        else:
            filter_olr = True
        if output['filters']['clara_ceres_time_diff']:
            filter_ceres = CLARA_CERES_combined['mean_time_diff'] < pd.Timedelta(seconds=output['filters']['clara_ceres_time_diff'])
        else:
            filter_ceres = True
        if output['filters']['clara_ceres_time_diff'] or output['filters']['pos_olr_time_diff']:
            combined_filter = filter_olr & filter_ceres
            CLARA_CERES_combined = CLARA_CERES_combined[combined_filter]
            CLARA_surface = CLARA_surface.loc[combined_filter[combined_filter].index]
            combination_key = combination_key[combined_filter]
            
        ### b) Save combined data ###
        # Generate filenames
        combined_filename = fh.create_filename(directory['combined'], date_filter_start, date_filter_end)
        netCDF_full_filename = directory['combined']['folder']+'/'+combined_filename+'.nc'
        pkl_full_filename = directory['combined']['folder']+'/pkl/'+combined_filename+'.pkl'
        csv_full_filename = directory['combined']['folder']+'/csv/'+combined_filename+'.csv'

        # Save NetCDF file 
        out.save_netCDF(CLARA_CERES_combined, CLARA_surface, netCDF_full_filename, parameters, output['filters'])
        del CLARA_surface

        # Save in .pkl or .csv format
        if output['save_intermediate']['combined'] in '.pkl':
            CLARA_CERES_combined.to_pickle(pkl_full_filename)
        if output['save_intermediate']['combined'] in '.csv':
            CLARA_CERES_combined.to_csv(csv_full_filename)
        if output['save_intermediate']['combined'] in 'both':
            CLARA_CERES_combined.to_pickle(pkl_full_filename)
            CLARA_CERES_combined.to_csv(csv_full_filename)
        # Print Step Completion
        function_timing.append(dt.now())
        print(f"    (6/6) Save the Data: {function_timing[-1]-function_timing[-2]}")
        print(f"    Saved to {netCDF_full_filename}")
        if output['save_intermediate']['combined'] in '.pkl':
            print('    ------------------------------------------------------')
            print(f"      Saved to {pkl_full_filename}")
        if output['save_intermediate']['combined'] in '.csv':
            print('    ------------------------------------------------------')
            print(f"      Saved to {csv_full_filename}")
        if output['save_intermediate']['combined'] in 'both':
            print('    ------------------------------------------------------')
            print(f"      Saved to {pkl_full_filename}")
            print(f"      Saved to {csv_full_filename}")

        """ 6: Finish the Loop """
        ### a) Free memory ###
        del CLARA_CERES_combined, combination_key
        gc.collect()
        ### b) Print total time for the iteration ###
        function_timing.append(dt.now())
        print(f"TOTAL: {function_timing[-1]-function_timing[0]}")
        print('##############################################################\n')
    return
