###############################################################################
""" Import the necessary Packages """
import numpy as np
import pandas as pd

from shapely.geometry import Point, Polygon, MultiPolygon

from scipy.spatial.transform import Rotation as R

###############################################################################
#%%
###############################################################################
""" Numerical & Unit Conversions """
###############################################################################
#%%
###############################################################################
def safe_float_conversion(value):
    """
    Safely converts a string to a float, returning 'No Float' if conversion fails.

    Parameters:
        value (str): The string to convert.

    Returns:
        float or str: Converted float value or "No Float" if conversion fails.
    """
    if value == 'nan':
        return "No Float"
    try:
        return float(value)
    except ValueError:
        return "No Float"
###############################################################################
#%%
###############################################################################
def latitude_to_colatitude(latitude):
    """Convert latitude (−90° to 90°) to colatitude (0° to 180°)."""
    return 90 - latitude
def colatitude_to_latitude(colatitude):
    """Convert colatitude (0° to 180°) back to latitude (−90° to 90°)."""
    return 90 - colatitude
###############################################################################
#%%
###############################################################################
def longitude_to_positive(longitude):
    """Convert longitude from the range (−180° to 180°) to (0° to 360°)."""
    return (longitude + 360) % 360
def positive_longitude_to_standard(longitude):
    """Convert longitude from the range (0° to 360°) to (−180° to 180°)."""
    return (longitude + 180) % 360 - 180
###############################################################################
#%%
###############################################################################
def clara_radiance_conversion(theta):
    """
    Conversion of the radiance from W/m2 to W/m2/sr depending on the opening
    angle of the CLARA fov.
    """
    theta_rad = theta/180*np.pi
    Omega = 2*np.pi*(1-np.cos(theta_rad))
    return 1/Omega
###############################################################################
#%%
###############################################################################
""" Quaternion & Vector Operations """
###############################################################################
#%%
###############################################################################
def create_quaternion(row):
    """
    Converts quaternion parameters to a normalized quaternion.
    Checks for unit and 0 quaternions.
    """
    qx, qy, qz, qw = row['qx'], row['qy'], row['qz'], row['qw']
    
    if ((qx == qy) & (qx == qz) & (qw == 1)):  # Invalid quaternion check
        return None
    
    norm = np.sqrt(qx**2 + qy**2 + qz**2 + qw**2)
    if norm == 0:  # Zero norm check
        return None
    
    return R.from_quat([qx, qy, qz, qw])
###############################################################################
#%%
###############################################################################
"""
Returns the satellite ECI coordinates based on the satellite position in ECI frame.
"""
def get_satellite_eci_coordinates(row):
    return np.array([-row['sx'], -row['sy'], -row['sz']])
###############################################################################
#%%
###############################################################################
"""
Returns the satellite pointing vector in the ECI frame based on quaternion.
"""
def get_satellite_pointing_vector(row, z_b):
    return row['quaternion_bo'].apply(z_b)
###############################################################################
#%%
###############################################################################
""" 
Computes the view zenith angle between the satellite's pointing vector 
and the satellite's position vector.
"""
def compute_view_zenith_angle(row):
    return np.arccos(np.dot(row['CLARA_satellite_eci_coordinates'] / np.linalg.norm(row['CLARA_satellite_eci_coordinates']), 
                           row['CLARA_satellite_pointing_vector_eci'])) * 180 / np.pi
###############################################################################
#%%
###############################################################################
""" Geospatiel Data Handling """
###############################################################################
#%%
###############################################################################
def extract_polygon_coords(polygons, max_vertices):
    """
    Converts a list of Shapely polygons and multipolygons into a NumPy array formatted for netCDF.
    
    Parameters:
        polygons (pd.Series): A Pandas Series containing Shapely Polygon objects.
        max_vertices (int): The maximum number of vertices in any polygon.
    
    Returns:
        np.ndarray: Array of shape (time, footprint_vertices, lon_lat)
    """
    num_obs = len(polygons)  # Number of time steps
    lon_lat = 2  # Latitude & Longitude
    
    # Initialize an empty array filled with NaN
    fov_array = np.full((num_obs, max_vertices, lon_lat), np.nan, dtype=np.float32)
    
    for i, poly in enumerate(polygons):
        if isinstance(poly, Polygon):
            coords = np.array(poly.exterior.coords[:])
        elif isinstance(poly, MultiPolygon):
            coords_list = []
            for geom in poly.geoms:
                coords_list.extend(geom.exterior.coords)
            unique_coords = []
            for c in coords_list:
                if c not in unique_coords:
                    unique_coords.append(c)
            coords = np.array(unique_coords)
        
        num_coords = len(coords)
    
        if num_coords > max_vertices:
            raise ValueError(f"Polygon at index {i} has more vertices ({num_coords}) than allowed ({max_vertices})!")
    
        fov_array[i, :num_coords, :] = coords[:max_vertices]
    
    return fov_array
###############################################################################
#%%
###############################################################################
def extract_polygons_from_multipolygon(geom):
    """ Extracts Polygon(s) from a Polygon or MultiPolygon. """
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    elif isinstance(geom, Polygon):
        return [geom]
    else:
        raise ValueError('Invalid geometry type: Must be Polygon or MultiPolygon')
