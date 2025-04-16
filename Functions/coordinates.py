###############################################################################
""" Import the necessary Packages """
import numpy as np
import pandas as pd

from datetime import datetime as dt
from shapely.geometry import Point, Polygon, MultiPolygon

import pymap3d.ecef as ecef

from timezonefinder import TimezoneFinder
from pytz import timezone

# Load functions from other files
from . import math_and_conversion as mc
###############################################################################
#%%
###############################################################################
""" Satellite Viewing Geometry """
###############################################################################
#%%
###############################################################################
def get_clara_fov_lonlat(view_vector_eci, sat_position_eci, timestamp, earth_radius=6371000):
    """
    Calculates the latitude and longitude of the point on Earth's surface 
    where the CLARA satellite is looking.

    Parameters:
        view_vector_eci (np.array): Unit vector representing the satellite's viewing direction (ECI coordinates).
        sat_position_eci (np.array): Satellite's position in ECI coordinates.
        timestamp (datetime): Time of observation, required for ECI to geodetic conversion.
        earth_radius (float, optional): Earth's radius in meters. Default is 6,371,000 m.

    Returns:
        tuple: (latitude, longitude) of the intersection point if it exists.
               (None, None) if no intersection occurs.
    """

    # Normalize the viewing direction vector
    view_vector_eci = view_vector_eci / np.linalg.norm(view_vector_eci)

    # Compute the discriminant (delta) to check for intersection
    delta = (
        (np.dot(view_vector_eci, sat_position_eci))**2 
        - (np.linalg.norm(sat_position_eci)**2 - earth_radius**2)
    )

    # If delta is negative, no intersection occurs (satellite is pointing away from Earth)
    if delta < 0:
        return (None, None)

    # Compute the distance from the satellite to the Earth's surface along the view vector
    d = -np.dot(view_vector_eci, sat_position_eci) + np.sqrt(delta)

    # Compute the intersection point in ECI coordinates
    surface_intersection_eci = sat_position_eci + d * view_vector_eci

    # Convert the intersection point from ECI to geodetic coordinates (lat, lon, alt)
    latitude, longitude, _ = ecef.eci2geodetic(
        -surface_intersection_eci[0], 
        -surface_intersection_eci[1], 
        -surface_intersection_eci[2], 
        timestamp
    )

    return longitude, latitude
###############################################################################
#%%
###############################################################################
def generate_fov_cone_vectors(view_vector, cone_angle_deg, num_vectors):
    """
    Generates a set of unit vectors forming a cone around a given direction.

    Parameters:
        view_vector (np.array): The main direction vector (unit vector).
        cone_angle_deg (float): The cone opening angle in degrees.
        num_vectors (int): The number of vectors to generate around the cone.

    Returns:
        np.array: An array of shape (num_vectors, 3) containing unit vectors on the cone.
    """

    # Normalize the main direction vector
    view_vector = view_vector / np.linalg.norm(view_vector)

    # Define an arbitrary perpendicular vector to create a basis
    basis_vector = np.array([1, 0, 0]) if abs(view_vector[0]) < 0.9 else np.array([0, 1, 0])

    # Compute two orthonormal basis vectors perpendicular to the view vector
    perp_vector1 = basis_vector - np.dot(basis_vector, view_vector) * view_vector
    perp_vector1 /= np.linalg.norm(perp_vector1)
    perp_vector2 = np.cross(view_vector, perp_vector1)

    # Convert cone angle from degrees to radians
    cone_angle_rad = np.radians(cone_angle_deg)

    # Generate vectors around the cone
    cone_vectors = []
    for i in range(num_vectors):
        phi = 2 * np.pi * i / num_vectors  # Angle around the cone
        cone_vector = (
            np.cos(cone_angle_rad) * view_vector +
            np.sin(cone_angle_rad) * (np.cos(phi) * perp_vector1 + np.sin(phi) * perp_vector2)
        )
        cone_vectors.append(cone_vector / np.linalg.norm(cone_vector))  # Ensure unit length

    return np.array(cone_vectors)
###############################################################################
#%%
###############################################################################
def get_fov_footprint(view_vector, satellite_position, time, cone_angle, num_vectors):
    """
    Calculates the footprint polygon of the area on Earth's surface 
    seen by the satellite's field of view (FOV).

    Parameters:
        view_vector (np.array): The satellite's main viewing direction (unit vector in ECI coordinates).
        satellite_position (np.array): The satellite's position in ECI coordinates.
        time (datetime): The observation time (needed for coordinate conversion).
        cone_angle (float): The field of view (FOV) cone opening angle in degrees.
        num_vectors (int): Number of vectors to sample around the cone.

    Returns:
        shapely.geometry.Polygon: The footprint polygon of the observed area on Earth's surface.
    """
    # Generate cone boundary vectors
    cone_vectors = generate_fov_cone_vectors(view_vector, cone_angle, num_vectors)
    
    # Compute intersection points on Earth's surface
    footprint_longitudes = []
    footprint_latitudes = []
    for vector in cone_vectors:
        longitude, latitude = get_clara_fov_lonlat(vector, satellite_position, time)
        if latitude is not None:
            footprint_longitudes.append(longitude)
            footprint_latitudes.append(latitude)
    
    # Close the polygon by repeating the first point
    if footprint_longitudes:
        footprint_longitudes.append(footprint_longitudes[0])
        footprint_latitudes.append(footprint_latitudes[0])
    
    footprint_poly = Polygon(zip(footprint_longitudes, footprint_latitudes)).buffer(0)
    
    footprint_colongitudes = [mc.longitude_to_positive(lon) for lon in footprint_longitudes]
    footprint_colatitudes = [mc.latitude_to_colatitude(lat) for lat in footprint_latitudes]
    
    if (all(lon <0 for lon in footprint_longitudes) or all(lon >=0 for lon in footprint_longitudes)):
        match_poly = MultiPolygon([Polygon(zip(footprint_colongitudes,footprint_colatitudes))])
    else:
        delta_lon = (max(footprint_longitudes) - min(footprint_longitudes))
        delta_colon = (max(footprint_colongitudes) - min(footprint_colongitudes))
        if delta_lon > delta_colon:
            match_poly = MultiPolygon([Polygon(zip(footprint_colongitudes,footprint_colatitudes))])
        else:
            full_poly = Polygon([Point(0,0), Point(360,0), Point(360,180), Point(0,180)])
            left_poly = full_poly.intersection(Polygon(zip(footprint_longitudes,footprint_colatitudes)).buffer(0))
    
            high_longitudes = [lon+360 for lon in footprint_longitudes]
            right_poly = full_poly.intersection(Polygon(zip(high_longitudes,footprint_colatitudes)).buffer(0))
            
            all_polys = (mc.extract_polygons_from_multipolygon(left_poly) 
                         + mc.extract_polygons_from_multipolygon(right_poly))
            
            match_poly = MultiPolygon(all_polys)
    
    return footprint_poly, match_poly
###############################################################################
#%%
###############################################################################
""" Time & Coordinate Handling """
###############################################################################
#%%
###############################################################################
def get_coordinates(data, theta, n, lim_z_angle):
    """
    Main function to calculate the satellite's coordinates, the field of view,
    and other related metrics like the footprint, zenith angle, and subsatellite position.
    
    Parameters:
        data (pd.DataFrame): DataFrame containing satellite quaternion and position data
        theta (float): Angle for generating vectors in the field of view cone
        n (int): Number of vectors to generate for the field of view
        lim_z_angle (float): Limiting zenith angle beyond which FOV is filtered
        
    Returns:
        pd.DataFrame: DataFrame with calculated coordinates and satellite information

    Functions:
    ---------
    coordinates.py (this file)
    - get_clara_fov_lonlat           - get_fov_footprint
    math_and_conversion.py:
    - create_quaternion,             - get_satellite_eci_coordinates
    - get_satellite_pointing_vector  - compute_view_zenith_angle
    - longitude_to_positive          - latitude_to_colatitude
    
    """
    # Set initial z-axis direction (satellite is pointing downward in ECI coordinates)
    z_b = [0, 0, -1]

    # Create quaternions from quaternion parameters
    data['quaternion_bo'] = data.apply(mc.create_quaternion, axis=1)
    data = data.dropna().copy()
    
    if data.empty:
        return pd.DataFrame()  # If data is empty after dropping NaN, return an empty DataFrame

    # Get satellite data: ECI coordinates and pointing vectors
    data['CLARA_satellite_eci_coordinates'] = data.apply(mc.get_satellite_eci_coordinates, axis=1)
    data['CLARA_satellite_pointing_vector_eci'] = data.apply(lambda row: mc.get_satellite_pointing_vector(row, z_b), axis=1)
    
    # Compute the view zenith angle (angle between satellite pointing vector and ECI coordinates)
    data['CLARA_view_zenith_angle'] = data.apply(mc.compute_view_zenith_angle, axis=1)
    
    # Filter out data where the zenith angle is beyond the limiting angle
    filtered_data = data[['CLARA_time_utc', 
                          'CLARA_satellite_eci_coordinates', 
                          'CLARA_satellite_pointing_vector_eci', 
                          'CLARA_view_zenith_angle'
                         ]][data['CLARA_view_zenith_angle'] < lim_z_angle]
    
    if filtered_data.empty:
        return pd.DataFrame()  # If data is empty after dropping NaN, return an empty DataFrame
    
    # Get the field of view coordinates (latitude, longitude, altitude)
    filtered_data[['CLARA_fov_longitude', 
                   'CLARA_fov_latitude']] = filtered_data.apply(
        lambda row: [val.item() for val in get_clara_fov_lonlat(
            row['CLARA_satellite_pointing_vector_eci'], 
            row['CLARA_satellite_eci_coordinates'], 
            row['CLARA_time_utc']
        )], axis=1, result_type='expand')

    # Check if the limiting zenith angle is valid (i.e., satellite is actually seeing the Earth)
    if filtered_data['CLARA_fov_latitude'].isna().any():
        print("Limiting z-angle is too large.")
        print("FOV satellite does not see the Earth.")
        print(f"Date: {filtered_data['CLARA_time_utc'][0]}")
        return pd.DataFrame()

    # Get the subsatellite position (coordinates directly under the satellite)
    filtered_data[['CLARA_subsatellite_longitude', 
                   'CLARA_subsatellite_latitude']] = filtered_data.apply(
        lambda row: [val.item() for val in get_clara_fov_lonlat(
            row['CLARA_satellite_eci_coordinates'], 
            row['CLARA_satellite_eci_coordinates'], 
            row['CLARA_time_utc']
        )], axis=1, result_type='expand')

    # Convert longitude [-180, 180] and latitude [-90, 90] to positive longitude [0, 360) and colatitude [0, 180]
    filtered_data['CLARA_fov_colatitude'] = filtered_data['CLARA_fov_latitude'].apply(mc.latitude_to_colatitude)
    filtered_data['CLARA_fov_longitude_positive'] = filtered_data['CLARA_fov_longitude'].apply(mc.longitude_to_positive)
    filtered_data['CLARA_subsatellite_colatitude'] = filtered_data['CLARA_subsatellite_latitude'].apply(mc.latitude_to_colatitude)
    filtered_data['CLARA_subsatellite_longitude_positive'] = filtered_data['CLARA_subsatellite_longitude'].apply(mc.longitude_to_positive)
    
    # Generate the FOV footprint polygon based on the satellite's pointing vector and position
    filtered_data[['CLARA_fov_footprint','CLARA_fov_footprint_match']] = filtered_data.apply(
        lambda row: get_fov_footprint(
            row['CLARA_satellite_pointing_vector_eci'], 
            row['CLARA_satellite_eci_coordinates'], 
            row['CLARA_time_utc'], theta, n
        ), axis=1, result_type='expand')

    return filtered_data
###############################################################################
#%%
###############################################################################
tf = TimezoneFinder()
def calculate_local_time_after_midnight(lon, lat, utc_time):
    """
    Calculates the time (in hours) after midnight for a given location on Earth in local time.

    Parameters:
        lon (float): Longitude of the location in degrees, in the range [-180, 180].
        lat (float): Latitude of the location in degrees, in the range [-90, 90].
        utc_time (datetime): The UTC time (a timezone-aware datetime object).

    Returns:
        float: The number of hours after midnight in local time, or None if the time zone cannot be found.
    
    Notes:
        The function uses the TimezoneFinder library to find the time zone based on the given longitude and latitude.
        Then, it converts the provided UTC time to the local time zone and computes the hours after midnight.
    """
    # Get the time zone name based on longitude and latitude
    tz_name = tf.timezone_at(lng=lon, lat=lat)
    
    if tz_name:
        # Convert the UTC time to local time based on the time zone
        local_tz = timezone(tz_name)
        local_time = utc_time.tz_localize('UTC').tz_convert(tz_name)
        
        # Calculate the hours after midnight
        hours_after_midnight = local_time.hour + local_time.minute / 60 + local_time.second / 3600
        
        return hours_after_midnight
    else:
        # If time zone is not found, return None
        return None
###############################################################################
#%%
###############################################################################
