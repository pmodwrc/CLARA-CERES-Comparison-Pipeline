###############################################################################
""" Import the necessary Packages """
import os
import gc

import numpy as np
import pandas as pd
import geopandas as gpd
import dask_geopandas as dgp

from datetime import datetime as dt
from datetime import timedelta
from shapely.geometry import Point, Polygon, MultiPolygon


###############################################################################
#%%
###############################################################################
""" Matching & Filtering """
###############################################################################
#%%
###############################################################################
def match_ceres_clara_positions(CERES_positions_df, CLARA_positions_df, input_file_structure):
    """
    Performs a spatial join between CERES and CLARA position data.

    This function converts longitude and latitude from the CERES dataset into 
    Point geometries and transforms both CERES and CLARA data into GeoDataFrames. 
    It then performs a spatial join to find CERES observation points that fall 
    within the CLARA footprint polygons.

    Parameters:
        CERES_positions_df (pd.DataFrame): DataFrame containing CERES position 
                                           data with longitude, colatitude, and time.
        CLARA_positions_df (pd.DataFrame): DataFrame containing CLARA footprint 
                                           polygons and observation times.
        input_file_structure: dict
            A dictionary defining the structure of the .nc files, including:
            - 'variable_names'

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame with CERES points matched to CLARA 
                          footprints where spatial overlap exists.
    """
    # Convert CERES longitude and latitude into a GeoPoint
    CERES_positions_df['instrument_fov_point'] = gpd.points_from_xy(
        CERES_positions_df[input_file_structure['variable_names']['fov_longitude']], 
        CERES_positions_df[input_file_structure['variable_names']['fov_colatitude']]
    )
    
    # Convert CERES DataFrame into a GeoDataFrame with a geometry column
    CERES_gdf = gpd.GeoDataFrame(
        CERES_positions_df[['instrument_fov_point']].reset_index().rename(
            columns={input_file_structure['variable_names']['time']: 'CERES_time'}
        ), 
        geometry='instrument_fov_point'
    )
    
    # Convert CLARA DataFrame into a GeoDataFrame (already contains a Polygon)
    CLARA_gdf = gpd.GeoDataFrame(
        CLARA_positions_df[['CLARA_fov_footprint_match', 'CLARA_time_utc']], 
        geometry='CLARA_fov_footprint_match'
    )
    CLARA_gdf = CLARA_gdf.explode(index_parts=False)
    
    # Bounding Box Filtering
    CLARA_gdf_bounds = CLARA_gdf.total_bounds
    CERES_gdf_filtered = CERES_gdf.cx[
    CLARA_gdf_bounds[0]:CLARA_gdf_bounds[2], 
    CLARA_gdf_bounds[1]:CLARA_gdf_bounds[3]
    ]
    
    # Convert to Dask-GeoPandas
    dask_CERES = dgp.from_geopandas(CERES_gdf_filtered, npartitions=10)
    dask_CLARA = dgp.from_geopandas(CLARA_gdf, npartitions=10)
    
    # Perform a spatial join: Keep only CERES points within CLARA footprints
    dask_spatial = dask_CERES.sjoin(dask_CLARA, predicate='within', how='inner')
    
    # Convert back to GeoPandas
    spatially_matched = dask_spatial.compute()
    
    spatially_matched = spatially_matched.drop(columns='instrument_fov_point')
    
    return spatially_matched
###############################################################################
#%%
###############################################################################
def get_nearest_time_deltas(time_deltas_series, time_bin_size):
    """
    Selects the lowest one or two unique time differences from a series.

    If the second smallest time difference is exactly 'time_bin_size' minutes 
    away from the smallest, both are kept. Otherwise, only the smallest is returned.
    This accounts for the rounding of the time differences

    Parameters:
        time_deltas_series (pd.Series): Series of unique time differences.
        time_bin_size (int): Time difference threshold in minutes.

    Returns:
        list: List containing one or two closest unique time differences.
    """
    # Sort unique time differences and convert to Timedelta
    unique_sorted_deltas = np.sort(time_deltas_series.unique())
    unique_sorted_deltas = list(map(pd.Timedelta, unique_sorted_deltas))

    # If the second smallest time difference is exactly the threshold apart, keep both
    if len(unique_sorted_deltas) > 1 and (
        unique_sorted_deltas[1] - unique_sorted_deltas[0] == pd.Timedelta(minutes=time_bin_size)
    ):
        return unique_sorted_deltas[:2]
    
    return [unique_sorted_deltas[0]]
###############################################################################
#%%
###############################################################################
def filter_ceres_clara_by_time(spatially_matched_df, time_bin_size):
    """
    Filters CERES points within CLARA footprints based on time difference.

    This function calculates the absolute time difference between CERES 
    and CLARA observation times, rounds it, and retains the closest matches. 
    It ensures that points that might have been excluded due to rounding 
    are included in the final filtered dataset.

    The function returns a DataFrame (`combination_key`) that maps CLARA 
    footprints to their corresponding CERES observations, including the 
    observation times, the number of CERES matches per footprint, and the 
    mean time difference. This mapping can be used for further analysis or 
    to merge with other datasets.

    Other Functions used:
        get_nearest_time_deltas()

    Parameters:
        spatially_matched_df (gpd.GeoDataFrame): GeoDataFrame containing 
                                                 spatially matched CERES and CLARA data.
        time_bin_size (int): Time difference threshold in minutes.

    Returns:
        pd.DataFrame: Filtered dataset containing CERES observations matched 
                      to CLARA footprints, grouped by CLARA index.
    """
    
    # Compute absolute time difference between CERES and CLARA observations
    spatially_matched_df['abs_time_diff'] = (
        spatially_matched_df['CERES_time'] - spatially_matched_df['CLARA_time_utc']
    ).abs()

    # Round time differences to nearest 'time_bin_size' minutes
    spatially_matched_df['rounded_time_diff'] = spatially_matched_df['abs_time_diff'].dt.round(f"{time_bin_size}min")

    # Identify closest unique time deltas for each CLARA footprint
    time_filter = (
        spatially_matched_df.groupby('index_right')['rounded_time_diff']
        .apply(lambda x: get_nearest_time_deltas(x, time_bin_size))
        .explode()
        .to_frame()
        .reset_index()
    )

    # Utilize the time_filter to filter the original dataset
    filtered_matches_df = spatially_matched_df.merge(
        time_filter, 
        on=['index_right', 'rounded_time_diff'], 
        how='inner'
    )
    # Delete Variables not used anymore to free up Memory
    del time_filter
    # Group by CLARA index to aggregate CERES times and compute mean time difference.
    combination_key = (
        filtered_matches_df.groupby('index_right')[['CERES_time', 'abs_time_diff']]
        .agg(list)
    )
    # Delete Variables not used anymore to free up Memory
    del filtered_matches_df
    # Compute the mean time difference per CLARA footprint
    combination_key['mean_time_diff'] = combination_key['abs_time_diff'].apply(np.mean)

    # Drop the time differences list, keeping only the mean
    combination_key = combination_key.drop('abs_time_diff', axis=1)

    # Count how many CERES points matched each CLARA footprint
    combination_key['num_CERES_matches'] = combination_key['CERES_time'].apply(len)

    # Rename index for clarity
    combination_key.index.name = 'CLARA_index'
    
    # Delete Variables not used anymore to free up Memory
    gc.collect()
    return combination_key
###############################################################################
#%%
###############################################################################
""" DataFrame Merging & Transformation """
###############################################################################
#%%
###############################################################################
def add_ceres_data_to_clara(CLARA_df, OLR_df, combination_key, CERES_radiance, CERES_angles, input_file_structure):
    """
    Enhances CLARA DataFrame by incorporating CERES-derived radiance and angle data.

    This function uses the combination key, which links CLARA footprints to 
    corresponding CERES observations, to calculate the mean CERES values for 
    each footprint. The computed values are then added to the CLARA dataset.

    Parameters:
        CLARA_df (pd.DataFrame): DataFrame containing CLARA position data.
        OLR_df (pd.DataFrame: DataFrame containing CLARA OLR data.
        combination_key (pd.DataFrame): DataFrame mapping CLARA footprints to 
                                        CERES observation times.
        CERES_radiance (pd.DataFrame): DataFrame containing CERES longwave radiance values.
        CERES_angles (pd.DataFrame): DataFrame containing CERES solar and view zenith angles.
        input_file_structure: dict
            A dictionary defining the structure of the .nc files, including:
            - 'variable_names'

    Returns:
        pd.DataFrame: The original CLARA DataFrame augmented with mean CERES values.
    """

    # Add the CLARA_OLR data
    enhanced_CLARA = pd.merge_asof(CLARA_df, OLR_df, left_on='CLARA_time_utc', right_on='CLARA_OLR_time', direction='nearest')
    enhanced_CLARA['CLARA_position_olr_time_diff'] = abs(enhanced_CLARA['CLARA_time_utc']-enhanced_CLARA['CLARA_OLR_time'])
    enhanced_CLARA = enhanced_CLARA.drop(columns=['OLR_index','JDAY','sav filename'])
    
    # Compute mean CERES longwave radiance for each CLARA footprint
    enhanced_CLARA['CERES_longwave_radiance'] = combination_key['CERES_time'].apply(
        lambda timestamps: CERES_radiance[input_file_structure['variable_names']['radiance']].loc[timestamps].mean() if timestamps else None
    )

    # Compute mean CERES solar zenith angle for each CLARA footprint
    enhanced_CLARA['CERES_solar_zenith_angle'] = combination_key['CERES_time'].apply(
        lambda timestamps: CERES_angles[input_file_structure['variable_names']['solar_zenith_angle']].loc[timestamps].mean() if timestamps else None
    )

    # Compute mean CERES view zenith angle for each CLARA footprint
    enhanced_CLARA['CERES_view_zenith_angle'] = combination_key['CERES_time'].apply(
        lambda timestamps: CERES_angles[input_file_structure['variable_names']['view_zenith_angle']].loc[timestamps].mean() if timestamps else None
    )

    # Compute mean CERES relative azimuth angle for each CLARA footprint
    enhanced_CLARA['CERES_relative_azimuth_angle'] = combination_key['CERES_time'].apply(
        lambda timestamps: CERES_angles[input_file_structure['variable_names']['relative_azimuth_angle']].loc[timestamps].mean() if timestamps else None
    )
    
    # Add mean_time_diff and num_CERES_matches
    enhanced_CLARA['mean_time_diff'] = combination_key['mean_time_diff']
    enhanced_CLARA['num_CERES_matches'] = combination_key['num_CERES_matches']
    enhanced_CLARA['CERES_time'] = combination_key['CERES_time']

    return enhanced_CLARA
###############################################################################
#%%
###############################################################################
def get_clara_surface(combination_key, CERES_surface, input_file_structure):
    """
    Generates a multi-index DataFrame for CLARA surface data using matched CERES surface data.

    This function takes the combination key (which links CLARA and CERES observations) and 
    aggregates CERES surface types for each CLARA footprint. It then normalizes the 
    surface type coverage percentages and ensures a consistent multi-index format.

    Parameters:
        combination_key (pd.DataFrame): DataFrame linking CLARA and CERES observations, 
                                        indexed by 'CLARA_index'.
        CERES_surface (pd.DataFrame): DataFrame containing CERES surface data.
        input_file_structure: dict
            A dictionary defining the structure of the .nc files, including:
            - 'variable_names'

    Returns:
        pd.DataFrame: Multi-index DataFrame where each CLARA footprint (CLARA_index) 
                      has associated surface type coverage percentages.
    """

    # Remove unnecessary columns from the combination key
    key_cleaned = combination_key.reset_index().drop(columns=['mean_time_diff', 'num_CERES_matches'])
    
    # Drop NaN values from CERES surface data and reset index for merging
    CERES_surface_cleaned = CERES_surface.dropna().reset_index()

    # Expand the combination key to have one row per CERES_time
    key_exploded = key_cleaned.explode('CERES_time')

    # Merge CERES surface data with the expanded combination key
    merged_surface = pd.merge(
        CERES_surface_cleaned, key_exploded, 
        left_on='time', right_on='CERES_time', 
        how='right'
    )
    
    # Remove rows without CLARA index matches and unnecessary columns
    filtered_surface = merged_surface.dropna(subset=['CLARA_index']).drop(columns=[input_file_structure['variable_names']['time'], input_file_structure['variable_names']['surface_index'], 'CERES_time'])
    
    # Group by CLARA footprint and surface type, summing up the coverage
    surface_coverage = filtered_surface.groupby(['CLARA_index', input_file_structure['variable_names']['surface_type']])[
        input_file_structure['variable_names']['surface_coverage']
    ].sum().to_frame()
    
    # Sort surface types within each CLARA footprint by coverage percentage (descending order)
    sorted_surface = surface_coverage.groupby('CLARA_index', group_keys=False).apply(
        lambda x: x.sort_values(by=input_file_structure['variable_names']['surface_coverage'], ascending=False)
    )
    
    # Normalize surface type coverage to percentages
    sorted_surface[input_file_structure['variable_names']['surface_coverage']] = (
        sorted_surface.groupby('CLARA_index')[input_file_structure['variable_names']['surface_coverage']]
        .transform(lambda x: 100 * x / x.sum())
    )

    # Reset index while keeping 'surface_igbp_type'
    sorted_surface_reset = sorted_surface.reset_index(level=input_file_structure['variable_names']['surface_type'], drop=False)
    
    # Assign a 'SurfType' rank within each CLARA footprint
    sorted_surface_reset[input_file_structure['variable_names']['surface_index']] = sorted_surface_reset.groupby('CLARA_index').cumcount()

    # Convert to a MultiIndex DataFrame indexed by ('CLARA_index', 'SurfType')
    CLARA_surface = sorted_surface_reset.reset_index().set_index(['CLARA_index', input_file_structure['variable_names']['surface_index']])
    
    # Ensure every CLARA footprint has a complete SurfType range, even if missing data
    all_surf_types = np.arange(CLARA_surface.index.get_level_values(level=input_file_structure['variable_names']['surface_index']).max() + 1)
    all_clara_indices = CLARA_surface.index.get_level_values(level='CLARA_index').unique()
    
    new_index = pd.MultiIndex.from_product([all_clara_indices, all_surf_types], names=['CLARA_index', input_file_structure['variable_names']['surface_index']])
    
    # Reindex to include missing (CLARA_index, SurfType) pairs
    CLARA_surface = CLARA_surface.reindex(new_index)
    CLARA_surface = CLARA_surface.rename(columns={
        input_file_structure['variable_names']['surface_index']: 'surface-igbp_type', 
        input_file_structure['variable_names']['surface_coverage']: 'surface_igbp_type_coverage'
    })
    
    return CLARA_surface
###############################################################################
#%%
###############################################################################
