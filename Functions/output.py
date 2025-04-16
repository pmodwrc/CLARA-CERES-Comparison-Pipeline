###############################################################################
""" Import the necessary Packages """
import os
import numpy as np
import pandas as pd

import netCDF4 as nc

# Load functions from other files
from . import math_and_conversion as mc
###############################################################################
#%%
###############################################################################
def save_netCDF(CLARA_CERES_combined, CLARA_surface, save_file_path, variables_dict, filters_dict):
    """
    Creates and saves a NetCDF file containing satellite position, radiance, 
    and surface data from CLARA and CERES satellites.

    Parameters:
    - CLARA_CERES_combined (DataFrame): Merged dataset of CLARA and CERES observations.
    - CLARA_surface (DataFrame): Surface classification data for different observation points.
    - save_file_path (str): Path where the NetCDF file will be saved.
    - variables_dic (dict): contains the variables that can be adjusted. Important here are:
        - max_view_zenith_angle (float): Limiting view zenith angle for the CLARA data analysis.
        - num_footprint_vertices (int): Number of vertices in the CLARA footprint polygon.
        - fov_angle (float): Opening angle of CLARA's field of view.
        - time_bin_size (int): Time bin size (e.g., 5 minutes) used for grouping the CERES observations.
    - filters_dict (dict): Filters applied to the data.
        - pos_olr_time_diff (float): the maximum time difference between the CLARA position and OLR time (in seconds)
        - clara_ceres_time_diff (float): the maximum time difference between the CLARA and CERES measurement (in minutes)

    Returns:
    - None (saves NetCDF file to specified path)
    """
    
    # Open a new NetCDF file
    ncfile = nc.Dataset(save_file_path, 'w', format='NETCDF4')

    """ Create the Dimensions """
    dim_time_len = len(CLARA_CERES_combined)
    if 'time' not in ncfile.dimensions:
        ncfile.createDimension('time', size=dim_time_len)
    
    dim_SurfType_len = CLARA_surface.index.get_level_values('SurfType').max()+1
    if 'SurfType' not in ncfile.dimensions:
        ncfile.createDimension('SurfType', size=dim_SurfType_len)
    
    dim_VectroXYZ_len = 3
    if 'VectorXYZ' not in ncfile.dimensions:
        ncfile.createDimension('VectorXYZ', size=dim_VectroXYZ_len)
    
    # Footprint vertices (+1 to close the polygon)
    dim_footprint_vertices_len = variables_dict['num_footprint_vertices']+1
    if 'footprint_vertices' not in ncfile.dimensions:
        ncfile.createDimension('footprint_vertices', size=dim_footprint_vertices_len)
    
    dim_lon_lat_len = 2
    if 'lon_lat' not in ncfile.dimensions:
        ncfile.createDimension('lon_lat', size=dim_lon_lat_len)
    
    dim_match_len = CLARA_CERES_combined['num_CERES_matches'].max()
    if 'match' not in ncfile.dimensions:
        ncfile.createDimension('match', size=dim_match_len)
    """ Global Variables """
    ########################### Create the Variables ##########################
    if 'time' not in ncfile.variables:
        time_var = ncfile.createVariable('time', 'f8', ('time',), zlib=True, 
                                         fill_value=np.finfo('f8').max)
        time_var.units = "days since 1970-01-01 00:00:00"
        time_var.calendar = "gregorian"
        time_var.long_name = "Time of observation (UTC)"
    if 'SurfType' not in ncfile.variables:
        SurfType = ncfile.createVariable('SurfType', 'i4', ('SurfType',), 
                                         zlib=True, fill_value=np.iinfo('i4').max)
        SurfType.long_name = "Index of the most prevalent surface types"
    ############################# Load in the Data ############################
    time_var[:] = (
        (CLARA_CERES_combined['CLARA_time_utc'] 
         - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)
    ).values
    SurfType[:] = np.arange(dim_SurfType_len)
    ###########################################################################
    """ Header """
    ############################# Create the Group ############################
    if 'Header' not in ncfile.groups:
        header_group = ncfile.createGroup('Header')
    else:
        header_group = ncfile.groups['Header']
    ########################### Create the Variables ##########################
    if 'number_of_measurements' not in header_group.variables:
        var_number_of_measurements = header_group.createVariable(
            'number_of_measurements', 'i8', zlib=True, fill_value=np.iinfo('i8').max
        )
        var_number_of_measurements.long_name = 'Total number of CLARA measurements'
        var_number_of_measurements.description = "The total count of CLARA measurements included in this netCDF file."
        var_number_of_measurements.coverage_content_type = 'referenceInformation'
    if 'clara_max_view_zenith_angle' not in header_group.variables:
        var_clara_max_view_zenith_angle = header_group.createVariable(
            'clara_max_view_zenith_angle', 'f4', zlib=True, fill_value=np.finfo('f4').max
        )
        var_clara_max_view_zenith_angle.long_name = 'Maximum allowed view zenith angle for CLARA observations'
        var_clara_max_view_zenith_angle.description = "The maximum view zenith angle, in degrees, that is allowed for CLARA observations in this dataset."
        var_clara_max_view_zenith_angle.units = 'degrees'
        var_clara_max_view_zenith_angle.coverage_content_type = 'referenceInformation'
    if 'num_footprint_vertices' not in header_group.variables:
        var_num_footprint_vertices = header_group.createVariable(
            'num_footprint_vertices', 'i2', zlib=True, fill_value=np.iinfo('i2').max
        )
        var_num_footprint_vertices.long_name = 'Number of vertices defining the CLARA footprint area'
        var_num_footprint_vertices.description = "The number of vertices used to define the CLARA footprint area."
        var_num_footprint_vertices.coverage_content_type = 'referenceInformation'
    if 'clara_fov_angle' not in header_group.variables:
        var_clara_fov_angle = header_group.createVariable(
            'clara_fov_angle', 'f4', zlib=True, fill_value=np.finfo('f4').max
        )
        var_clara_fov_angle.long_name = 'Half-angle of the CLARA field of view (FOV)'
        var_clara_fov_angle.description = "The half-angle, in degrees, of the CLARA field of view (FOV) cone. This defines the angular extent of the satellite's observation region."
        var_clara_fov_angle.units = 'degrees'
        var_clara_fov_angle.coverage_content_type = 'referenceInformation'
    if 'ceres_matching_time_bin_size' not in header_group.variables:
        var_ceres_matching_time_bin_size = header_group.createVariable(
            'ceres_matching_time_bin_size', 'f4', zlib=True, fill_value=np.finfo('f4').max
        )
        var_ceres_matching_time_bin_size.long_name = 'Time threshold (in minutes) for binning CERES matches'
        var_ceres_matching_time_bin_size.description = "The time threshold, in minutes, used to group CERES matches into bins."
        var_ceres_matching_time_bin_size.units = 'minutes'
        var_ceres_matching_time_bin_size.coverage_content_type = 'referenceInformation'
    if filters_dict['pos_olr_time_diff']:
        if 'filter_pos_olr_time_diff' not in header_group.variables:
            var_filter_pos_olr_time_diff = header_group.createVariable(
                'filter_pos_olr_time_diff', 'f4', zlib=True, fill_value=np.finfo('f4').max
            )
            var_filter_pos_olr_time_diff.long_name = 'Maximum allowed time difference between CLARA position time and OLR time'
            var_filter_pos_olr_time_diff.description = "The maximum allowed time difference, in seconds, between the CLARA position time and the corresponding OLR measurement time."
            var_filter_pos_olr_time_diff.units = 'seconds'
            var_filter_pos_olr_time_diff.coverage_content_type = 'referenceInformation'
        var_filter_pos_olr_time_diff[:] = filters_dict['pos_olr_time_diff']
    if filters_dict['clara_ceres_time_diff']:
        if 'filter_clara_ceres_time_diff' not in header_group.variables:
            var_filter_clara_ceres_time_diff = header_group.createVariable(
                'filter_clara_ceres_time_diff', 'f4', zlib=True, fill_value=np.finfo('f4').max
            )
            var_filter_clara_ceres_time_diff.long_name = 'Maximum allowed mean time difference between CLARA position and CERES time'
            var_filter_clara_ceres_time_diff.description = "The maximum allowed mean time difference, in minutes, between the CLARA position time and the corresponding CERES observation time."
            var_filter_clara_ceres_time_diff.units = 'minutes'
            var_filter_clara_ceres_time_diff.coverage_content_type = 'referenceInformation'
        var_filter_clara_ceres_time_diff[:] = filters_dict['clara_ceres_time_diff']
    ############################# Load in the Data ############################
    var_number_of_measurements[:] = dim_time_len
    var_clara_max_view_zenith_angle[:] = variables_dict['max_view_zenith_angle']
    var_num_footprint_vertices[:] = variables_dict['num_footprint_vertices']
    var_clara_fov_angle[:] = variables_dict['fov_angle']
    var_ceres_matching_time_bin_size[:] = variables_dict['time_bin_size']
    ###########################################################################
    """ Time and Position """
    ############################# Create the Group ############################
    if 'Time_and_Position' not in ncfile.groups:
        time_and_position_group = ncfile.createGroup('Time_and_Position')
    else:
        time_and_position_group = ncfile.groups['Time_and_Position']
    ########################### Create the Variables ##########################
    if 'instrument_fov_longitude' not in time_and_position_group.variables:
        var_instrument_fov_longitude = time_and_position_group.createVariable(
            'instrument_fov_longitude', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_instrument_fov_longitude.long_name = 'Longitude of instrument field of view at surface'
        var_instrument_fov_longitude.description = "Geodetic longitude of the instrument's surface field of view, expressed in degrees east from 0 to 360. This represents the surface intercept location observed at each time step."
        var_instrument_fov_longitude.units = 'degrees_east'
        var_instrument_fov_longitude.coverage_content_type = 'coordinate'
        var_instrument_fov_longitude.valid_range = np.array([0.0, 360.0], dtype='f4')
        var_instrument_fov_longitude.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'instrument_fov_latitude' not in time_and_position_group.variables:
        var_instrument_fov_latitude = time_and_position_group.createVariable(
            'instrument_fov_latitude', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_instrument_fov_latitude.long_name = 'Latitude of instrument field of view at surface'
        var_instrument_fov_latitude.description = "Geodetic latitude of the instrument's surface field of view, expressed in degrees north from -90 to 90. This location corresponds to the surface intercept point observed at each time step."
        var_instrument_fov_latitude.units = 'degrees_north'
        var_instrument_fov_latitude.coverage_content_type = 'coordinate'
        var_instrument_fov_latitude.valid_range = np.array([-90.0, 90.0], dtype='f4')
        var_instrument_fov_latitude.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'instrument_fov_colatitude' not in time_and_position_group.variables:
        var_instrument_fov_colatitude = time_and_position_group.createVariable(
            'instrument_fov_colatitude', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_instrument_fov_colatitude.long_name = 'Colatitude of instrument field of view at surface'
        var_instrument_fov_colatitude.description = "Geodetic colatitude of the instrument's surface field of view, expressed in degrees from 0° at the North Pole to 180° at the South Pole. This is complementary to latitude and may be used in spherical coordinate systems."
        var_instrument_fov_colatitude.units = 'degrees'
        var_instrument_fov_colatitude.coverage_content_type = 'physicalMeasurement'
        var_instrument_fov_colatitude.valid_range = np.array([0.0, 180.0], dtype='f4')
        var_instrument_fov_colatitude.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'instrument_fov_footprint' not in time_and_position_group.variables:
        var_instrument_fov_footprint = time_and_position_group.createVariable(
            'instrument_fov_footprint', 'f4', ('time', 'footprint_vertices', 'lon_lat'), zlib=True, fill_value=np.finfo('f4').max
        )
        var_instrument_fov_footprint.long_name = 'Surface footprint polygon of instrument field of view'
        var_instrument_fov_footprint.description = "Geodetic coordinates (longitude, latitude) of the surface footprint vertices defining the field of view at each time step."
        var_instrument_fov_footprint.units = 'degrees'
        var_instrument_fov_footprint.coverage_content_type = 'coordinate'
        var_instrument_fov_footprint.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'UTC_time' not in time_and_position_group.variables:
        var_UTC_time = time_and_position_group.createVariable(
            'UTC_time', 'f8', ('time',), zlib=True, fill_value=np.finfo('f8').max
        )
        var_UTC_time.long_name = 'Time of observation (UTC) in number of days since the beginning of the Julian period'
        var_UTC_time.description = "Time of each observation expressed in UTC, as days since 1970-01-01 00:00:00. This variable follows CF conventions and can be interpreted using standard calendar handling."
        var_UTC_time.units = 'days since 1970-01-01 00:00:00'
        var_UTC_time.coverage_content_type = 'physicalMeasurement'
        var_UTC_time.valid_range = np.array([10000.0, 400000.0], dtype='f8')
        var_UTC_time.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'local_time' not in time_and_position_group.variables:
        var_local_time = time_and_position_group.createVariable(
            'local_time', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_local_time.long_name = 'Local solar time at observation point (in hours after midnight)'
        var_local_time.description = "Local time at the surface observation point, expressed in decimal hours after local midnight (0–24). This value is derived from UTC time and the observation longitude."
        var_local_time.units = 'hours'
        var_local_time.coverage_content_type = 'auxiliaryInformation'
        var_local_time.valid_range = np.array([0.0, 24.0], dtype='f4')
        var_local_time.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'solar_zenith_angle' not in time_and_position_group.variables:
        var_solar_zenith_angle = time_and_position_group.createVariable(
            'solar_zenith_angle', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_solar_zenith_angle.long_name = 'Solar zenith angle at observed location'
        var_solar_zenith_angle.description = "Solar zenith angle at the geodetic location and time of observation, expressed in degrees."
        var_solar_zenith_angle.units = 'degrees'
        var_solar_zenith_angle.coverage_content_type = 'auxiliaryInformation'
        var_solar_zenith_angle.valid_range = np.array([0.0, 180.0], dtype='f4')
        var_solar_zenith_angle.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    var_instrument_fov_longitude[:] = CLARA_CERES_combined['CLARA_fov_longitude_positive'].values
    var_instrument_fov_latitude[:] = CLARA_CERES_combined['CLARA_fov_latitude'].values
    var_instrument_fov_colatitude[:] = CLARA_CERES_combined['CLARA_fov_colatitude'].values
    var_instrument_fov_footprint[:] = mc.extract_polygon_coords(CLARA_CERES_combined['CLARA_fov_footprint'].values, variables_dict['num_footprint_vertices']+1)
    var_UTC_time[:] = ((CLARA_CERES_combined['CLARA_time_utc'] - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)).values
    var_local_time[:] = CLARA_CERES_combined['CLARA_local_time'].values
    var_solar_zenith_angle[:] = CLARA_CERES_combined['CLARA_solar_zenith_angle'].values
    ###########################################################################
    """ Radiance """
    ############################# Create the Group ############################
    if 'Radiance' not in ncfile.groups:
        radiance_group = ncfile.createGroup('Radiance')
    else:
        radiance_group = ncfile.groups['Radiance']
    ########################### Create the Variables ##########################
    if 'longwave_radiance_clara' not in radiance_group.variables:
        var_longwave_radiance_clara = radiance_group.createVariable(
            'longwave_radiance_clara', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_longwave_radiance_clara.long_name = 'Outgoing longwave radiance measured by CLARA'
        var_longwave_radiance_clara.description = "Top-of-atmosphere outgoing longwave radiance measured by the CLARA instrument, expressed in watts per square meter per steradian. This represents the upwelling thermal emission in the longwave spectral range at the time and location of observation."
        var_longwave_radiance_clara.units = 'W/m^2/sr'
        var_longwave_radiance_clara.coverage_content_type = 'physicalMeasurement'
        var_longwave_radiance_clara.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'longwave_radiance_ceres' not in radiance_group.variables:
        var_longwave_radiance_ceres = radiance_group.createVariable(
            'longwave_radiance_ceres', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_longwave_radiance_ceres.long_name = 'Matched outgoing longwave radiance form CERES'
        var_longwave_radiance_ceres.description = "Mean top-of-atmosphere outgoing longwave radiance from CERES observations, spatially and temporally matched to the CLARA observation. Values are expressed in watts per square meter per steradian."
        var_longwave_radiance_ceres.units = 'W/m^2/sr'
        var_longwave_radiance_ceres.coverage_content_type = 'physicalMeasurement'
        var_longwave_radiance_ceres.valid_range = np.array([0.0, 200.0], dtype='f4')
        var_longwave_radiance_ceres.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    var_longwave_radiance_clara[:] = CLARA_CERES_combined['CLARA_radiance'].values
    var_longwave_radiance_ceres[:] = CLARA_CERES_combined['CERES_longwave_radiance'].values
    ###########################################################################
    """ Surface Map """
    ############################# Create the Group ############################
    if 'Surface_Map' not in ncfile.groups:
        surface_map_group = ncfile.createGroup('Surface_Map')
    else:
        surface_map_group = ncfile.groups['Surface_Map']
    ########################### Create the Variables ##########################
    if 'surface_igbp_type' not in surface_map_group.variables:
        var_surface_igbp_type = surface_map_group.createVariable(
            'surface_igbp_type', 'i2', ('time', 'SurfType'), zlib=True, fill_value=np.iinfo('i2').max
        )
        var_surface_igbp_type.long_name = 'Surface type index (IGBP type)'
        var_surface_igbp_type.description = "International Geosphere-Biosphere Programme (IGBP) land cover type indices present within the instrument's field of view at each observation time. Each index corresponds to a specific land cover class (e.g., evergreen needleleaf forest, grassland, urban)."
        var_surface_igbp_type.coverage_content_type = 'physicalMeasurement'
        var_surface_igbp_type.valid_range = np.array([1,20], dtype='i2')
        var_surface_igbp_type.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
        var_surface_igbp_type.flag_value = np.array([1,20], dtype='i2')
        var_surface_igbp_type.flag_meanings = 'evergreen_needleleaf_forest evergreen_broadleaf_forest deciduous_needleleaf_forest deciduous_broadleaf_forest mixed_forest closed_shrublands open_shrublands woody_savannas savannas grasslands permanent_wetlands croplands urban_and_built-up cropland_mosaics snow_and_ice_permanent bare_soil_and_rocks water_bodies tundra fresh_snow sea_ice'
    if 'surface_igbp_type_coverage' not in surface_map_group.variables:
        var_surface_igbp_type_coverage = surface_map_group.createVariable(
            'surface_igbp_type_coverage', 'i2', ('time', 'SurfType'), zlib=True, fill_value=np.iinfo('i2').max
        )
        var_surface_igbp_type_coverage.long_name = 'Surface type percent coverage'
        var_surface_igbp_type_coverage.description = "Percentage of each International Geosphere-Biosphere Programme (IGBP) land cover type present within the instrument's field of view at each observation time. Coverage values are expressed as percentages (0 to 100%) of the observed area corresponding to each surface type index."
        var_surface_igbp_type_coverage.units = 'percent'
        var_surface_igbp_type_coverage.coverage_content_type = 'physicalMeasurement'
        var_surface_igbp_type_coverage.valid_range = np.array([0, 100], dtype='i2')
        var_surface_igbp_type_coverage.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    # Replace NaN values with the fill value
    CLARA_surface = CLARA_surface.fillna(np.iinfo('i2').max)
    var_surface_igbp_type[:] = CLARA_surface['surface_igbp_type'].unstack().values
    var_surface_igbp_type_coverage[:] = CLARA_surface['surface_igbp_type_coverage'].unstack().values
    ###########################################################################
    """ CLARA Satellite """
    ############################# Create the Group ############################
    if 'CLARA_Satellite' not in ncfile.groups:
        clara_satellite_group = ncfile.createGroup('CLARA_Satellite')
    else:
        clara_satellite_group = ncfile.groups['CLARA_Satellite']
    ########################### Create the Variables ##########################
    if 'subsatellite_longitude' not in clara_satellite_group.variables:
        var_subsatellite_longitude = clara_satellite_group.createVariable(
            'subsatellite_longitude', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_subsatellite_longitude.long_name = 'Longitude of the subsatellite point'
        var_subsatellite_longitude.description = "Longitude of the surface point directly beneath the CLARA satellite, expressed in degrees east (0–360). This value represents the location on the Earth's surface where the satellite is directly overhead at each time step."
        var_subsatellite_longitude.units = 'degrees_east'
        var_subsatellite_longitude.coverage_content_type = 'auxiliaryInformation'
        var_subsatellite_longitude.valid_range = np.array([0.0, 360.0], dtype='f4')
        var_subsatellite_longitude.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'subsatellite_colatitude' not in clara_satellite_group.variables:
        var_subsatellite_colatitude = clara_satellite_group.createVariable(
            'subsatellite_colatitude', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_subsatellite_colatitude.long_name = 'Colatitude of the subsatellite point'
        var_subsatellite_colatitude.description = "Colatitude of the surface point directly beneath the CLARA satellite, expressed in degrees (0–180), where 0° corresponds to the North Pole and 180° corresponds to the South Pole. This value represents the location on the Earth's surface where the satellite is directly overhead at each time step."
        var_subsatellite_colatitude.units = 'degrees'
        var_subsatellite_colatitude.coverage_content_type = 'auxiliaryInformation'
        var_subsatellite_colatitude.valid_range = np.array([0.0, 180.0], dtype='f4')
        var_subsatellite_colatitude.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'satellite_position_eci' not in clara_satellite_group.variables:
        var_satellite_position_eci = clara_satellite_group.createVariable(
            'satellite_position_eci', 'f8', ('time', 'VectorXYZ'), zlib=True, fill_value=np.finfo('f8').max
        )
        var_satellite_position_eci.long_name = 'Satellite position in Earth-Centered Inertial (ECI) coordinates'
        var_satellite_position_eci.description = "Position of the satellite in ECI (Earth-Centered Inertial) coordinates at each time step. The position is provided as a 3-component vector, with values for the X, Y, and Z coordinates, expressed in kilometers (km). These components define the satellite's location relative to the Earth's center."
        var_satellite_position_eci.units = 'km'
        var_satellite_position_eci.coverage_content_type = 'auxiliaryInformation'
        var_satellite_position_eci.coordinates = 'Time_and_Position/UTC_time'
    if 'satellite_pointing_vector_eci' not in clara_satellite_group.variables:
        var_satellite_pointing_vector_eci = clara_satellite_group.createVariable(
            'satellite_pointing_vector_eci', 'f8', ('time', 'VectorXYZ'), zlib=True, fill_value=np.finfo('f8').max
        )
        var_satellite_pointing_vector_eci.long_name = 'Normalized satellite pointing vector in Earth-Centered Inertial (ECI) coodrinates'
        var_satellite_pointing_vector_eci.description = "Direction in which the satellite is pointing, represented as a unit vector in ECI coordinates at each time step. The vector defines the orientation of the satellite relative to the Earth's center, with components in the X, Y, and Z directions. The vector is normalized to unit length."
        var_satellite_pointing_vector_eci.coverage_content_type = 'auxiliaryInformation'
        var_satellite_pointing_vector_eci.valid_range = np.array([-1.0,1.0], dtype='f8')
        var_satellite_pointing_vector_eci.coordinates = 'Time_and_Position/UTC_time'
    if 'view_zenith_angle' not in clara_satellite_group.variables:
        var_view_zenith_angle = clara_satellite_group.createVariable(
            'view_zenith_angle', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_view_zenith_angle.long_name = 'CLARA viewing zenith angle at surface'
        var_view_zenith_angle.description = "The angle between the satellite's viewing vector (or pointing vector) and the nadir vector, expressed in degrees. This represents the zenith angle of the instrument's view relative to the Earth's surface, where 0° corresponds to nadir (directly below the satellite), and larger values indicate a more oblique view."
        var_view_zenith_angle.units = 'degrees'
        var_view_zenith_angle.coverage_content_type = 'physicalMeasurement'
        var_view_zenith_angle.valid_range = np.array([0.0, variables_dict['max_view_zenith_angle']], dtype='f4')
        var_view_zenith_angle.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    var_subsatellite_longitude[:] = CLARA_CERES_combined['CLARA_subsatellite_longitude_positive'].values
    var_subsatellite_colatitude[:] = CLARA_CERES_combined['CLARA_subsatellite_colatitude'].values
    var_satellite_position_eci[:] = np.vstack(CLARA_CERES_combined['CLARA_satellite_eci_coordinates'].values)
    var_satellite_pointing_vector_eci[:] = np.vstack(CLARA_CERES_combined['CLARA_satellite_pointing_vector_eci'].values)
    var_view_zenith_angle[:] =CLARA_CERES_combined['CLARA_view_zenith_angle'].values
    ###########################################################################
    """ CERES Satellite """
    ############################# Create the Group ############################
    if 'CERES_Satellite' not in ncfile.groups:
        ceres_satellite_group = ncfile.createGroup('CERES_Satellite')
    else:
        ceres_satellite_group = ncfile.groups['CERES_Satellite']
    ########################### Create the Variables ##########################
    if 'mean_solar_zenith_angle' not in ceres_satellite_group.variables:
        var_mean_solar_zenith_angle = ceres_satellite_group.createVariable(
            'mean_solar_zenith_angle', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_mean_solar_zenith_angle.long_name = 'Mean solar zenith angle from CERES observations'
        var_mean_solar_zenith_angle.description = "The average solar zenith angle at the surface during matched CERES observations, expressed in degrees."
        var_mean_solar_zenith_angle.units = 'degrees'
        var_mean_solar_zenith_angle.coverage_content_type = 'auxiliaryInformation'
        var_mean_solar_zenith_angle.valid_range = np.array([0.0, 180.0], dtype='f4')
        var_mean_solar_zenith_angle.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'mean_relative_azimuth_angle' not in ceres_satellite_group.variables:
        var_mean_relative_azimuth_angle = ceres_satellite_group.createVariable(
            'mean_relative_azimuth_angle', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_mean_relative_azimuth_angle.long_name = 'Mean relative azimuth angle from CERES observations'
        var_mean_relative_azimuth_angle.description = "The average relative azimuth angle at the surface during matched CERES observations, expressed in degrees. The relative azimuth angle represents the angle between the Sun-satellite plane and the satellite-observation direction."
        var_mean_relative_azimuth_angle.units = 'degrees'
        var_mean_relative_azimuth_angle.coverage_content_type = 'auxiliaryInformation'
        var_mean_relative_azimuth_angle.valid_range = np.array([0.0, 360.0], dtype='f4')
        var_mean_relative_azimuth_angle.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'mean_view_zenith_angle' not in ceres_satellite_group.variables:
        var_mean_view_zenith_angle = ceres_satellite_group.createVariable(
            'mean_view_zenith_angle', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_mean_view_zenith_angle.long_name = 'Mean view zenith angle from CERES observations'
        var_mean_view_zenith_angle.description = "The average view zenith angle at the surface during matched CERES observations, expressed in degrees."
        var_mean_view_zenith_angle.units = 'degrees'
        var_mean_view_zenith_angle.coverage_content_type = 'auxiliaryInformation'
        var_mean_view_zenith_angle.valid_range = np.array([0.0, 90.0], dtype='f4')
        var_mean_view_zenith_angle.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    var_mean_solar_zenith_angle[:] = CLARA_CERES_combined['CERES_solar_zenith_angle'].values
    var_mean_relative_azimuth_angle[:] = CLARA_CERES_combined['CERES_relative_azimuth_angle'].values
    var_mean_view_zenith_angle[:] = CLARA_CERES_combined['CERES_view_zenith_angle'].values
    ###########################################################################
    """ Matching Information """
    ############################# Create the Group ############################
    if 'Matching_Information' not in ncfile.groups:
        matching_information_group = ncfile.createGroup('Matching_Information')
    else:
        matching_information_group = ncfile.groups['Matching_Information']
    ########################### Create the Variables ##########################
    if 'ceres_match_times' not in matching_information_group.variables:
        var_ceres_match_times = matching_information_group.createVariable(
            'ceres_match_times', 'f8', ('time', 'match'), zlib=True, fill_value=np.finfo('f8').max
        )
        var_ceres_match_times.long_name = 'Timestamps of matched CERES observations'
        var_ceres_match_times.description = "Timestamps of the CERES instrument observations that correspond to each CLARA observation, expressed as days since 1970-01-01 00:00:00 (UNIX epoch). These times are aligned to the matched CLARA observations."
        var_ceres_match_times.units = 'days since 1970-01-01 00:00:00'
        var_ceres_match_times.coverage_content_type = 'auxiliaryInformation'
        var_ceres_match_times.valid_range = np.array([10000.0, 400000.0], dtype='f8')
        var_ceres_match_times.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'num_ceres_matches' not in matching_information_group.variables:
        var_num_ceres_matches = matching_information_group.createVariable(
            'num_ceres_matches', 'i8', ('time',), zlib=True, fill_value=np.iinfo('i8').max
        )
        var_num_ceres_matches.long_name = 'Number of matched CERES observations'
        var_num_ceres_matches.description = "The total number of CERES observations that correspond to each CLARA observation. This value indicates how many CERES data points are spatially and temporally matched to the CLARA observation."
        var_num_ceres_matches.coverage_content_type = 'auxiliaryInformation'
        var_num_ceres_matches.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'mean_ceres_match_time_diff' not in matching_information_group.variables:
        var_mean_ceres_match_time_diff = matching_information_group.createVariable(
            'mean_ceres_match_time_diff', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_mean_ceres_match_time_diff.long_name = 'Mean time difference between CLARA and CERES observations'
        var_mean_ceres_match_time_diff.description = "The average time difference, in seconds, between the CERES observations and the corresponding CLARA observation."
        var_mean_ceres_match_time_diff.units = 'seconds'
        var_mean_ceres_match_time_diff.coverage_content_type = 'auxiliaryInformation'
        var_mean_ceres_match_time_diff.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'olr_measurement_time' not in matching_information_group.variables:
        var_olr_measurement_time = matching_information_group.createVariable(
            'olr_measurement_time', 'f8', ('time',), zlib=True, fill_value=np.finfo('f8').max
        )
        var_olr_measurement_time.long_name = 'Time of CLARA OLR measurement'
        var_olr_measurement_time.description = "The time, in days since 1970-01-01 00:00:00, of the CLARA outgoing longwave radiation (OLR) measurement that corresponds to the matched CLARA position time."
        var_olr_measurement_time.units = 'days since 1970-01-01 00:00:00'
        var_olr_measurement_time.coverage_content_type = 'auxiliaryInformation'
        var_olr_measurement_time.valid_range = np.array([10000.0, 400000.0], dtype='f8')
        var_olr_measurement_time.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    if 'olr_position_time_diff' not in matching_information_group.variables:
        var_olr_position_time_diff = matching_information_group.createVariable(
            'olr_position_time_diff', 'f4', ('time',), zlib=True, fill_value=np.finfo('f4').max
        )
        var_olr_position_time_diff.long_name = 'Time difference between the UTC_time of the CLARA position and the OLR measurement time'
        var_olr_position_time_diff.description = "The time difference, in seconds, between the CLARA outgoing longwave radiation (OLR) measurement time and the corresponding CLARA position time."
        var_olr_position_time_diff.units = 'seconds'
        var_olr_position_time_diff.coverage_content_type = 'auxiliaryInformation'
        var_olr_position_time_diff.coordinates = 'Time_and_Position/UTC_time  Time_and_Position/instrument_fov_latitude  Time_and_Position/instrument_fov_longitude'
    ############################# Load in the Data ############################
    CERES_times = CLARA_CERES_combined[['CERES_time']].explode('CERES_time').reset_index()
    CERES_times = CERES_times.set_index(['index', CERES_times.groupby('index').cumcount()])
    var_ceres_match_times[:] = ((CERES_times['CERES_time'] - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)).unstack().fillna(np.finfo('f8').max).values
    var_num_ceres_matches[:] = CLARA_CERES_combined['num_CERES_matches'].values
    var_mean_ceres_match_time_diff[:] = CLARA_CERES_combined['mean_time_diff'].dt.total_seconds().values
    var_olr_measurement_time[:] = ((CLARA_CERES_combined['CLARA_OLR_time'] - pd.Timestamp("1970-01-01")) / pd.Timedelta(days=1)).values
    var_olr_position_time_diff[:] = CLARA_CERES_combined['CLARA_position_olr_time_diff'].dt.total_seconds().values
    
    ncfile.close()
    return