# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:36:37 2024

@author: lrussell
"""

import pandas as pd
import numpy as np
from shapely.geometry import Point,LineString,Polygon,MultiPolygon,shape,box
import geopandas as gpd
import rasterio as rio
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio import features
from rasterio.warp import Resampling
from rasterio.warp import reproject,calculate_default_transform
from rasterio.io import MemoryFile
from pysheds.grid import Grid
from scipy.integrate import trapezoid
from scipy.interpolate import RBFInterpolator,interp1d
import matplotlib.pyplot as plt
from matplotlib import colors
import copy
try:
    from segment_anything import SamAutomaticMaskGenerator, sam_model_registry
    from skimage.measure import regionprops_table
except:
    print('Missing dependencies to run pebble counts.')
    
class Terrain:
    def __init__(self, grid, dem, fdir, acc, crs):
        """
        Initialize a Terrain object with essential hydrological and geographical data.

        Parameters:
        grid : object
            The spatial grid object containing the DEM data.
        dem : ndarray
            Digital Elevation Model as a 2D array.
        fdir : ndarray
            Flow direction data derived from the DEM.
        acc : ndarray
            Flow accumulation data, showing how much water accumulates over the surface.
        crs : object
            Coordinate Reference System of the DEM data.
        """
        
        self.grid = grid
        self.dem = dem
        self.fdir = fdir
        self.acc = acc
        self.crs = crs

def terrain(dem_file):
    """
    Processes a DEM file to generate a terrain object that includes pit-filled DEM, flow directions,
    and flow accumulation.

    Parameters:
    dem_file : str
        Path to the DEM raster file.

    Returns:
    terrain : Terrain object
        A terrain object containing attributes for DEM, flow direction, accumulation, and CRS.
    """
    
    # Load the DEM grid
    grid = Grid.from_raster(dem_file)  # Ensure Grid.from_raster is a valid method in your library
    dem = grid.read_raster(dem_file)  # Reads the raster data into an array
    
    # Extract coordinate reference system from DEM
    with rio.open(dem_file) as src:
        crs = src.crs
        
    print('Grid loaded')
    
    # Hydrological preprocessing of the DEM
    # Fill all pits in the DEM to prevent artificial sinks
    pit_filled_dem = grid.fill_pits(dem)
    print('Pits filled')
    
    # Fill all depressions in the DEM to ensure proper flow routing
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    print('Depressions filled')
    
    # Resolve flat areas in the DEM, crucial for accurate flow direction determination
    inflated_dem = grid.resolve_flats(flooded_dem, eps=1e-12, max_iter=1e9)
    print('Flats resolved')
    
    # Calculate flow direction using the resolved DEM
    fdir = grid.flowdir(inflated_dem)
    print('Flow direction computed')
    
    # Calculate flow accumulation based on the flow direction
    acc = grid.accumulation(fdir)
    print('Flow accumulation computed')
    
    # Create a Terrain object which includes the processed DEM and hydrological attributes
    terrain_object = Terrain(grid, dem, fdir, acc, crs)
    
    return terrain_object

def streamlines(terrain, pour_point, threshold=None, snap_threshold=None, save=None):
    """
    Extract stream networks from a digital elevation model (DEM) based on flow accumulation thresholds and pour points.

    Parameters:
    terrain : object
        Terrain object containing grid, accumulation array, flow direction array, and crs.
    pour_point : str or GeoDataFrame
        Path to a shapefile or a GeoDataFrame containing the pour point geometry.
    threshold : float, optional
        Flow accumulation threshold for defining streams. Default is grid size/100.
    snap_threshold : float, optional
        Threshold for snapping the pour point to the nearest high accumulation cell. Default is grid size/10.
    save : str, optional
        File path where the resulting stream network GeoDataFrame is saved. If None, returns the GeoDataFrame.

    Returns:
    streams : GeoDataFrame or None
        GeoDataFrame of extracted stream networks, unless saved to file.
    """

    grid = copy.deepcopy(terrain.grid)
    
    # Set default thresholds if not provided
    if threshold is None:
        threshold = grid.size / 100
    if snap_threshold is None:
        snap_threshold = grid.size / 10

    # Load and snap pour points
    if isinstance(pour_point, str):
        pour_point = gpd.read_file(pour_point)
    pour_point = pour_point.to_crs(terrain.crs)
    pour_point = pour_point.geometry.values[0]
    
    xy = grid.snap_to_mask(terrain.acc > snap_threshold, [pour_point.x, pour_point.y], return_dist=False)
    x, y = xy[0, 0], xy[0, 1]

    # Delineate the catchment area
    catchment = grid.catchment(x=x, y=y, fdir=terrain.fdir, xytype='coordinate')
    grid.clip_to(catchment)

    # Extract stream network using flow direction and accumulation arrays
    streams = grid.extract_river_network(terrain.fdir, terrain.acc > threshold)
    
    # Convert the extracted streams to a GeoDataFrame
    lines = gpd.GeoDataFrame(streams, columns=['geometry'], crs=terrain.crs)
    lines['FID'] = range(len(lines))
    
    # Find upstream and downstream neighbors for each stream segment
    for idx, geom in lines.iterrows():
        # Define upstream and downstream points
        upstream_point = Point(geom.geometry.coords[0])
        downstream_point = Point(geom.geometry.coords[-1])
        
        # Identify neighboring stream segments
        lines['upstream'] = lines.geometry.apply(lambda x: x.distance(upstream_point) < 0.01 and not np.array_equal(x.coords[-1], geom.geometry.coords[0]))
        lines['downstream'] = lines.geometry.apply(lambda x: x.distance(downstream_point) < 0.01 and not np.array_equal(x.coords[0], geom.geometry.coords[-1]))

    # Calculate distances to pour point using a recursive function
    def calculate_distance_to_pour_point(line_fid, accumulated_distance=0):
        current_line = lines.loc[line_fid]
        lines.loc[line_fid, 'distance_to_pour_point'] = accumulated_distance
        downstream_fid = lines.loc[line_fid, 'downstream']
        
        if not downstream_fid.empty:
            next_distance = accumulated_distance + current_line.geometry.length
            calculate_distance_to_pour_point(downstream_fid.values[0], next_distance)

    # Start with the stream segment closest to the pour point
    pour_point_fid = lines[lines['downstream'].isnull()].index[0]
    calculate_distance_to_pour_point(pour_point_fid)

    if save:
        # Save the stream network to a file
        lines.to_file(save)
        print(f'Stream network saved to {save}')
    else:
        # Return the stream network GeoDataFrame
        return lines
    
def catchment(terrain, pour_point, snap_threshold=None, save=None):
    """
    Delineate a catchment area from a terrain model using a specified pour point.

    Parameters:
    terrain : Terrain object
        Terrain object must have attributes 'grid', 'crs', 'acc', and 'fdir'.
    pour_point : str or GeoDataFrame
        Path to a shapefile or a GeoDataFrame containing the pour point geometry.
    snap_threshold : float, optional
        Threshold for snapping the pour point to the nearest high accumulation cell.
        If None, defaults to one-tenth of the grid size.
    save : str, optional
        Path where the catchment shapefile is saved. If None, returns the catchment as a GeoDataFrame.

    Returns:
    cm : GeoDataFrame or None
        GeoDataFrame representing the catchment polygon, unless saved to file.
    """

    # Deep copy the grid to ensure no modification to the original data
    grid = copy.deepcopy(terrain.grid)

    # Set a default snapping threshold if none provided
    if snap_threshold is None:
        snap_threshold = grid.size / 10

    # Load the pour point if provided as a path, and transform to the terrain CRS
    if isinstance(pour_point, str):
        pour_point = gpd.read_file(pour_point)
    pour_point = pour_point.to_crs(terrain.crs)
    pour_point = pour_point.geometry.values[0]

    # Extract coordinates of the pour point
    x0, y0 = pour_point.x, pour_point.y

    # Find the nearest flow accumulated cell that exceeds the snap threshold
    xy = grid.snap_to_mask(terrain.acc > snap_threshold, np.column_stack([x0, y0]), return_dist=False)
    x, y = xy[0, 0], xy[0, 1]

    # Delineate the catchment using the snapped coordinates
    catch = grid.catchment(x=x, y=y, fdir=terrain.fdir, xytype='coordinate')

    # Clip the grid to the delineated catchment
    grid.clip_to(catch)

    # Polygonize the clipped grid to get the catchment area
    shapes = grid.polygonize()
    coords = np.asarray(next(shapes)[0]['coordinates'][0])
    catchment_poly = Polygon(coords)

    # Create a GeoDataFrame for the catchment polygon
    cm = gpd.GeoDataFrame([], geometry=[catchment_poly], crs=terrain.crs)

    # Save or return the GeoDataFrame
    if save:
        cm.to_file(save)
        print('Catchment delineated and saved to:', save)
    else:
        print('Catchment delineated')
        return cm
    
def get_xs(cl, xs_length, spacing, save=None):
    """
    Generate perpendicular cross-sectional lines along a centerline at specified intervals.

    Parameters:
    cl : str or GeoDataFrame
        Path to the centerline shapefile or a GeoDataFrame containing the centerline geometries.
    xs_length : float
        Half the total length of each cross-sectional line. The full length will be twice this value.
    spacing : float
        Distance between the centers of consecutive cross-sectional lines along the centerline.
    save : str, optional
        File path where the resulting cross-sectional lines GeoDataFrame is saved. If None, the function returns the GeoDataFrame.

    Returns:
    xs_lines : GeoDataFrame or None
        GeoDataFrame containing the generated cross-sectional lines, unless saved to file.
    """

    # Load centerline from file if a path is provided
    if isinstance(cl, str):
        cl = gpd.read_file(cl)

    # Prepare an empty GeoDataFrame to store cross-sectional lines
    xs_lines = gpd.GeoDataFrame([], columns=['CSID', 'Distance', 'geometry'], crs=cl.crs)

    # Iterate through each line in the centerline GeoDataFrame
    for idx, row in cl.iterrows():
        num_xsecs = int(row.geometry.length / spacing)
        
        # Generate cross-sections for each segment defined by the spacing
        for xs in range(1, num_xsecs):
            point = row.geometry.interpolate(xs * spacing)
            pointA = row.geometry.interpolate(xs * spacing - xs_length)
            pointB = row.geometry.interpolate(xs * spacing + xs_length)
            
            # Calculate the angle of the line to find the perpendicular
            deltaX = pointA.x - pointB.x
            deltaY = pointA.y - pointB.y
            if deltaX == 0:
                deltaX = 0.0000000000001  # Prevent division by zero
            slope = deltaY / deltaX
            theta = np.arctan(slope)
            new_theta = (theta + np.pi / 2) % np.pi  # Perpendicular angle
            
            # Calculate endpoints of the cross-sectional line
            line_end1 = Point(point.x + xs_length * np.cos(new_theta), point.y + xs_length * np.sin(new_theta))
            line_end2 = Point(point.x - xs_length * np.cos(new_theta), point.y - xs_length * np.sin(new_theta))
            
            line = LineString([line_end1, line_end2])
            
            # Store the cross-section in the GeoDataFrame
            xs_lines.loc[len(xs_lines)] = [idx, xs * spacing, line]
    
    # Optionally save to file or return the GeoDataFrame
    if save:
        xs_lines.to_file(save)
        print(f"Cross-sectional lines saved to {save}")
    else:
        return xs_lines

def get_points(lines, spacing, ep=True, save=None):
    """
    Generate points at specified intervals along line geometries.

    Parameters:
    lines : str or GeoDataFrame
        Path to the line shapefile or a GeoDataFrame containing line geometries.
    spacing : float
        Distance between points along the lines.
    ep : bool, optional
        If True, includes the endpoint of the line. Default is True.
    save : str, optional
        Path where the resulting points GeoDataFrame is saved. If None, the function returns the GeoDataFrame.

    Returns:
    points : GeoDataFrame or None
        GeoDataFrame of points if not saved to file.
    """

    # Load lines from file if a path is provided
    if isinstance(lines, str):
        lines = gpd.read_file(lines)

    # Define a function to create points including the endpoint
    def _get_points_ep(row):
        # Create points including the end point of the line
        intervals = np.arange(0, row.geometry.length, spacing).tolist() + [row.geometry.length]
        points = [row.geometry.interpolate(distance) for distance in intervals]
        return gpd.GeoDataFrame({
            'CSID': row.name,
            'Station': intervals,
            'geometry': points
        }, crs=lines.crs)

    # Define a function to create points excluding the endpoint
    def _get_points(row):
        # Create points excluding the end point of the line
        intervals = np.arange(0, row.geometry.length, spacing)
        points = [row.geometry.interpolate(distance) for distance in intervals]
        return gpd.GeoDataFrame({
            'CSID': row.name,
            'Station': intervals,
            'geometry': points
        }, crs=lines.crs)

    # Apply the appropriate function to each line
    if ep:
        lines['points'] = lines.apply(_get_points_ep, axis=1)
    else:
        lines['points'] = lines.apply(_get_points, axis=1)

    # Concatenate all points into a single GeoDataFrame
    points = pd.concat(lines['points'].tolist(), ignore_index=True)

    # Save or return the points GeoDataFrame
    if save:
        points.to_file(save)
        print(f"Points saved to {save}")
    else:
        return points

def inundate(raster, rel_wse_list, largest_only=False, invert=False, remove_holes=False, save=None):
    """
    Create inundation polygons from a raster based on relative water surface elevation (WSE) levels.

    Parameters:
    raster : str
        Path to the raster file.
    rel_wse_list : list
        List of relative water surface elevations to map.
    largest_only : bool, optional
        If True, only the largest polygon for each WSE level is retained.
    invert : bool, optional
        If True, inverts the inundation logic to find non-inundated areas.
    remove_holes : bool, optional
        If True, removes holes in the polygons.
    save : str, optional
        File path to save the resulting GeoDataFrame as a file.

    Returns:
    inundated : GeoDataFrame or None
        GeoDataFrame containing the inundation polygons, unless saved to file.
    """

    # Open the raster and read essential properties
    with rio.open(raster) as src:
        array = src.read(1)
        affine = src.transform
        nodata = src.nodata
        crs = src.crs

    # Mask the array where nodata values are present
    masked_array = np.ma.masked_where(array == nodata, array)
    inundated = gpd.GeoDataFrame([], crs=crs, columns=['WSE', 'geometry'])

    # Process each water surface elevation
    for wse in rel_wse_list:
        # Apply the inundation condition
        if invert:
            inundation_array = np.ma.where(masked_array <= wse, 0, 1)
        else:
            inundation_array = np.ma.where(masked_array <= wse, 1, 0)
        inundation_array = inundation_array.astype(np.uint8)

        # Extract polygons from the binary array
        if largest_only:
            largest_poly = Polygon([])
            for shapedict, value in features.shapes(inundation_array, mask=inundation_array == 1, transform=affine):
                polygon = shape(shapedict).buffer(0)
                if polygon.area > largest_poly.area:
                    largest_poly = polygon
            inundated.loc[len(inundated)] = [wse, largest_poly]
        else:
            polys = [shape(shapedict).buffer(0) for shapedict, value in features.shapes(inundation_array, mask=inundation_array == 1, transform=affine) if value == 1]
            geom = MultiPolygon(polys)
            inundated.loc[len(inundated)] = [wse, geom]

        print(f'Completed inundation mapping for {wse} feet above WSE')

    # Remove holes in polygons if requested
    if remove_holes:
        def _remove_holes(poly):
            if poly.geom_type == 'Polygon':
                return Polygon(poly.exterior.coords)
            elif poly.geom_type == 'MultiPolygon':
                return MultiPolygon([Polygon(geom.exterior.coords) for geom in poly.geoms])
            return poly

        inundated.geometry = inundated.geometry.apply(_remove_holes)

    # Save or return the GeoDataFrame
    if save:
        inundated.to_file(save)
        print(f'Inundation data saved to {save}')
    else:
        return inundated

def edit_raster(raster, output, crs=None, resample=None, resample_method='bilinear', clip=None, nodata=None):
    """
    Edit a raster for various transformations including clipping, nodata handling, resampling, and reprojection.

    Parameters:
    raster : str
        Path to the input raster file.
    output : str
        Path to save the modified raster file.
    crs : CRS or str, optional
        The target coordinate reference system to which to reproject the raster.
    resample : float or None, optional
        Scaling factor to resample the raster; less than 1 to downscale, more than 1 to upscale.
    resample_method : str, optional
        Method used for resampling ('bilinear', 'nearest', 'mean').
    clip : str or GeoDataFrame, optional
        Path to a shapefile or a GeoDataFrame used for clipping the raster.
    nodata : float or int, optional
        New nodata value to set in the raster.

    Returns:
    None
    """

    # Define resampling methods
    mdict = {
        'bilinear': Resampling.bilinear,
        'nearest': Resampling.nearest,
        'mean': Resampling.average
    }

    # Open the source raster and initialize the first MemoryFile
    with rio.open(raster) as src:
        profile = src.profile.copy()
        if nodata is not None:
            profile['nodata'] = nodata
        profile['driver'] = 'GTiff'
        
        with MemoryFile() as memfile1:
            with memfile1.open(**profile) as mem_dst:
                mem_dst.write(src.read())
        
            # Handling the clip operation
            if clip is not None:
                with MemoryFile() as memfile2:
                    with memfile1.open() as src:
                        if isinstance(clip, str):
                            clip = gpd.read_file(clip)
                        clip = clip.to_crs(src.crs)
                        masked_array, masked_transform = rio.mask.mask(src, clip.geometry, crop=True)
                        profile.update({
                            'transform': masked_transform,
                            'width': masked_array.shape[2],
                            'height': masked_array.shape[1]
                        })
                        with memfile2.open(**profile) as dst:
                            dst.write(masked_array)
                    memfile1 = memfile2
            
            # Handling the resample operation
            if resample:
                with MemoryFile() as memfile2:
                    with memfile1.open() as src:
                        new_height = int(src.height * resample)
                        new_width = int(src.width * resample)
                        data = src.read(
                            out_shape=(src.count, new_height, new_width),
                            resampling=mdict[resample_method]
                        )
                        transform = src.transform * src.transform.scale(
                            (src.width / data.shape[-1]),
                            (src.height / data.shape[-2])
                        )
                        profile.update({
                            'transform': transform,
                            'width': new_width,
                            'height': new_height
                        })
                        with memfile2.open(**profile) as dst:
                            dst.write(data)
                    memfile1 = memfile2
            
            # Handling the reprojection operation
            if crs:
                with MemoryFile() as memfile2:
                    with memfile1.open() as src:
                        transform, width, height = calculate_default_transform(
                            src.crs, crs, src.width, src.height, *src.bounds
                        )
                        profile.update({
                            'crs': crs,
                            'transform': transform,
                            'width': width,
                            'height': height
                        })
                        with memfile2.open(**profile) as dst:
                            for i in range(1, src.count + 1):
                                reproject(
                                    source=rio.band(src, i),
                                    destination=rio.band(dst, i),
                                    src_transform=src.transform,
                                    src_crs=src.crs,
                                    dst_transform=transform,
                                    dst_crs=crs,
                                    resampling=mdict[resample_method]
                                )
                    memfile1 = memfile2

            # Saving the final raster
            with memfile1.open() as final_src:
                with rio.open(output, 'w', **final_src.profile) as final_dst:
                    final_dst.write(final_src.read())

    print('Raster saved to {0}'.format(output))

def merge_rasters(rasters, output, method='first', resample_method='bilinear', compression=None, nodata=None):
    """
    Merge multiple raster files into a single output raster.

    Parameters:
    rasters: list
        List of paths to the raster files to be merged.
    output: str
        Path to save the merged output raster.
    method: str, optional
        Method for merging rasters ('first', 'last', 'min', 'max', 'mean').
    resample_method: str, optional
        Resampling method used during merging ('bilinear', 'nearest', 'mean').
    compression: str, optional
        Compression method for the output file (e.g., 'lzw', 'deflate').
    nodata: int or float, optional
        Value representing nodata in the output raster.

    Returns:
    None
    """
    # Mapping of resampling methods
    mdict = {
        'bilinear': Resampling.bilinear,
        'nearest': Resampling.nearest,
        'mean': Resampling.average
    }
    
    open_datasets = []
    for i, raster in enumerate(rasters):
        src = rio.open(raster)
        open_datasets.append(src)
        
        # Copy the profile from the first raster
        if i == 0:
            profile = src.profile.copy()
    
    # If nodata value is not provided, use the nodata value from the first raster
    if not nodata:
        nodata = profile['nodata']
    
    # Merge the rasters
    merged_array, merged_transform = merge(open_datasets, method=method, nodata=nodata, resampling=mdict[resample_method])
    
    # Update profile with new dimensions and transform
    profile['transform'] = merged_transform
    profile['width'] = merged_array.shape[2]
    profile['height'] = merged_array.shape[1]
    profile['driver'] = 'GTiff'
    profile['nodata'] = nodata
    
    # Apply compression if specified
    if compression:
        profile['compress'] = compression
    
    # Close all open datasets
    for src in open_datasets:
        src.close()
    
    # Write the merged raster to the output file
    with rio.open(output, 'w', **profile, BIGTIFF='YES') as dst:
        dst.write(merged_array)
        
    print(f'Merged rasters at {output}')

def difference_rasters(r1, r2, output, match_affine='first', method='nearest'):
    """
    Calculate the difference between two raster files and save the result as a new raster.

    Parameters:
    r1: str
        Path to the first raster file.
    r2: str
        Path to the second raster file.
    output: str
        Path to save the output difference raster.
    match_affine: str, optional
        Which raster's affine transform to match ('first' or 'last').
    method: str, optional
        Resampling method ('nearest' or 'bilinear').

    Returns:
    None
    """
    # Define resampling methods
    mdict = {'bilinear': Resampling.bilinear,
             'nearest': Resampling.nearest}
    
    # Open the first raster
    with rio.open(r1) as src1:
        profile1 = src1.profile.copy()
        
        # Open the second raster
        with rio.open(r2) as src2:
            profile2 = src2.profile.copy()
            
            # Determine the target profile based on match_affine argument
            if match_affine == 'first':
                target_profile = profile1
                src = src2
                a1 = src1.read(1)
            else:
                target_profile = profile2
                src = src1
                a2 = src2.read(1)
            
            # Use MemoryFile to reproject the raster to match the target profile
            with MemoryFile() as memfile:
                dst_match = memfile.open(**target_profile)
             
                # Reproject the raster to match the affine transform and CRS
                for i in range(1, src.count + 1):
                    reproject(
                        source=rio.band(src, i),
                        destination=rio.band(dst_match, i),
                        src_transform=src.profile['transform'],
                        src_crs=src.crs,
                        dst_transform=target_profile['transform'],
                        dst_crs=target_profile['crs'],
                        resampling=mdict[method])
                    
                # Read the reprojected raster and prepare for difference calculation
                if match_affine == 'first':
                    a2 = dst_match.read(1)
                    a1_mask = np.where(a1 == src1.nodata, 0, 1)
                    a2_mask = np.where(a2 == dst_match.nodata, 0, 1)
                elif match_affine == 'last':
                    a1 = dst_match.read(1)
                    a1_mask = np.where(a1 == dst_match.nodata, 0, 1)
                    a2_mask = np.where(a2 == src2.nodata, 0, 1)
                else:
                    return print('Incorrect affine target: first or last')
                 
            # Calculate the difference between the two rasters
            final = a2 - a1
            
            # Mask the difference raster to handle nodata values
            final_mask = np.where(a1_mask & a2_mask, 1, 0)
            final = np.where(final_mask, final, target_profile['nodata'])
            
            # Save the final difference raster
            with rio.open(output, 'w', **target_profile) as dst:
                dst.write(final, indexes=1)
                
def zonal_stats(df, raster, metric='mean', profile=None, band=1):
    """
    Calculate zonal statistics for vector geometries on a raster.

    Parameters:
    df: GeoDataFrame
        GeoDataFrame containing the vector geometries.
    raster: str or ndarray
        Path to the raster file or a NumPy array representing the raster.
    metric: str, optional
        The statistic to compute ('mean', 'min', 'max', 'sum', 'nonzero').
    profile: dict, optional
        Profile information required if `raster` is provided as an ndarray.
    band: int, optional
        The band of the raster to read (default is 1).

    Returns:
    stats: ndarray
        Array of calculated statistics for each geometry in the GeoDataFrame.
    """

    # Load the raster or use the provided ndarray
    if isinstance(raster, np.ndarray):
        a = raster
        if not profile:
            raise Exception("Need affine transformation with array")
        df = df.to_crs(profile['crs'])
        transform = profile['transform']
    else:
        with rio.open(raster) as src:
            a = src.read(band)
            transform = src.transform
            a = np.where(a == src.nodata, np.nan, a)
            df = df.to_crs(src.crs)

    # Dictionary to map metric names to NumPy functions
    met_dict = {
        'mean': np.nanmean,
        'min': np.nanmin,
        'max': np.nanmax,
        'sum': np.nansum,
        'nonzero': np.count_nonzero
    }
    
    # Function to apply the metric function to the raster within the geometry
    def run_func(func, geom, arr, trans):
        # Create a mask where the geometry overlaps the raster
        mask = rio.features.geometry_mask([geom], arr.shape, trans, invert=True)
        
        # If all values in the masked area are NaN, return NaN
        if np.isnan(arr[mask]).all():
            return np.nan
        
        # Apply the metric function to the masked array
        val = func(arr[mask])
        return val
    
    # Apply the run_func to each geometry in the GeoDataFrame
    stats = df.geometry.apply(lambda geom: run_func(met_dict[metric], geom, a, transform))
    stats = np.array(stats).astype(float)

    return stats

def sample_raster(raster, points, buff=None, metric='min', multiband=False, crs=None):
    """
    Sample values from a raster at specified point locations with optional buffering and metric calculation.

    Parameters:
    raster: str
        Path to the raster file.
    points: str or GeoDataFrame
        Path to the points shapefile or a GeoDataFrame of points.
    buff: float, optional
        Buffer distance (in units of CRS) around each point for sampling.
    metric: str, optional
        Metric to calculate within the buffered area ('min', 'max', 'mean', 'mode').
    multiband: bool, optional
        If True, samples values from all raster bands at each point.
    crs: CRS or str, optional
        Coordinate Reference System to reproject points.

    Returns:
    sampled_data: list
        List of sampled values (or list of lists if multiband is True).
    """
    
    def _get_metric(buff_geom, src, metric):
        """
        Compute the specified metric within the buffered geometry on the raster.
        """
        masked_array = mask(src, [buff_geom], crop=True)[0]
        
        # Combine conditions for nodata and NaN values
        condition1 = (masked_array == src.nodata)
        condition2 = np.isnan(masked_array)
        combined_condition = np.logical_or(condition1, condition2)
        
        masked_array = np.ma.masked_where(combined_condition, masked_array)
        
        # Calculate the desired metric
        if metric == 'min':
            value = masked_array.min()
        elif metric == 'max':
            value = masked_array.max()
        elif metric == 'mean':
            value = masked_array.mean()
        elif metric == 'mode':
            value = masked_array.mode()
        return value
    
    # Load the points if a file path is provided
    if isinstance(points, str):
        points = gpd.read_file(points)
    
    with rio.open(raster) as src:
        # Reproject points to match the raster CRS if necessary
        if crs:
            points = points.to_crs(crs)
        else:
            points = points.to_crs(src.crs)
        
        # Buffer the points if a buffer distance is provided
        if buff:
            points.geometry = points.geometry.buffer(buff)
            sampled_data = points.geometry.apply(lambda buffer_geom: _get_metric(buffer_geom, src, metric))
        else:
            # Sample values directly at points, with optional multiband sampling
            if multiband:
                sampled_data = [x.tolist() for x in src.sample(zip(points.geometry.x, points.geometry.y), masked=True)]
            else:
                sampled_data = [x.tolist()[0] for x in src.sample(zip(points.geometry.x, points.geometry.y), masked=True)]
    
    return sampled_data

def create_REM(dem, xs, output, Zcol=None, sample_dist=3, smooth_window=5, buffer=1000, limits=[-50, 50], method='min', vb=None, wse_path=None):
    """
    Create a Relative Elevation Model (REM) from a DEM and cross-sections.

    Parameters:
    dem: str
        Path to the DEM raster file.
    xs: str or GeoDataFrame
        Path to the cross-sections file or a GeoDataFrame of cross-sections.
    output: str
        Path where the output REM raster file will be saved.
    Zcol: str, optional
        Column name for the elevation values in the cross-sections. If None, elevation will be sampled from the DEM.
    sample_dist: float, optional
        Distance between sample points along the cross-sections (in units of CRS).
    smooth_window: int, optional
        Window size for smoothing elevation values along cross-sections.
    buffer: int, optional
        Buffer distance around cross-sections for clipping the DEM.
    limits: list, optional
        Minimum and maximum values for REM; values outside these limits will be set to -9999.
    method: str, optional
        Method to aggregate sampled elevations ('min' or 'mean').
    vb: str, optional
        Path to the valley bottom polygon shapefile.
    wse_path: str, optional
        Path to save the Water Surface Elevation (WSE) raster file.

    Returns:
    None
        Saves the REM raster and optionally the WSE raster to the specified paths.
    """

    # Load DEM and get CRS
    with rio.open(dem) as src:
        crs = src.crs
    
    # Load cross-sections and clip to the valley bottom polygon if provided
    if isinstance(xs, str):
        xs = gpd.read_file(xs).to_crs(crs)
    else:
        xs = xs.to_crs(crs)
    
    if vb:
        vb = gpd.read_file(vb).to_crs(crs)
        xs.geometry = xs.geometry.intersection(vb.geometry.unary_union)
    
    # Sample points along cross-sections to get elevation
    if Zcol:
        xs['Z'] = xs[Zcol]
    else:
        xsp = get_points(xs, sample_dist)
        xsp['Z'] = sample_raster(dem, xsp)
        
        if method == 'min':    
            xs['Z'] = xsp.groupby('CSID')['Z'].min()
        elif method == 'mean':
            xs['Z'] = xsp.groupby('CSID')['Z'].mean()
        
        for i, group in xsp.groupby('CSID'):
            xs.loc[xs.CSID == i, 'Z'] = xs.loc[xs.CSID == i, 'Z'].rolling(smooth_window, center=True, min_periods=1).mean()
    
    # Create points for interpolation along cross-sections
    def _get_int_points(row):
        intervals = np.linspace(0, row.geometry.length, num=4)
        temp = gpd.GeoDataFrame([], columns=['CSID', 'Station', 'geometry'], crs=crs)
        temp['geometry'] = [row.geometry.interpolate(x) for x in intervals]
        temp['Station'] = intervals
        temp['CSID'] = row.name
        return temp
    
    int_points = xs.apply(lambda row: _get_int_points(row), axis=1)
    int_points = gpd.GeoDataFrame(pd.concat(int_points.values), crs=crs)
    int_points['Z'] = int_points.CSID.map(dict(zip(xs.index, xs.Z)))
    int_points = int_points.dropna()
    int_points.index = range(len(int_points))
    int_points = int_points.iloc[int_points.round(1).drop_duplicates(subset=['CSID', 'Station']).index, :].reset_index(drop=True)

    # Clip the raster to the buffer of the cross-sections
    with rio.open(dem, 'r') as src:
        masked_array, masked_transform = mask(src, shapes=[box(*xs.total_bounds).buffer(buffer)], crop=True, nodata=-9999, indexes=1)
        nd_mask = (masked_array == -9999)
        
        # Update the masked raster profile
        masked_profile = src.meta.copy()
        masked_profile.update({
            "driver": "GTiff",
            "height": masked_array.shape[0],
            "width": masked_array.shape[1],
            "transform": masked_transform,
            "nodata": -9999
        })
           
    # Create linear RBF interpolator
    int_r, int_c = rio.transform.rowcol(masked_transform, int_points.geometry.x, int_points.geometry.y)
    rc = pd.DataFrame(np.array([int_c, int_r]).T, columns=['rows', 'cols']).drop_duplicates()
    rbfi = RBFInterpolator(rc.values, int_points.loc[rc.index, 'Z'], kernel='linear')
    
    # Interpolate WSE over cells within the buffered mask
    ci, ri = np.meshgrid(np.arange(masked_array.shape[1]), np.arange(masked_array.shape[0]))
    ci = ci[~nd_mask]
    ri = ri[~nd_mask]
    
    int_wse = rbfi(np.stack([ci, ri]).T)
    wse = masked_array.copy()
    wse[~nd_mask] = int_wse
    
    # Calculate REM by subtracting WSE from DEM
    rem = masked_array - wse
    
    # Apply limits to the REM
    if limits:
        rem = np.where((rem > limits[1]) | (rem < limits[0]), -9999, rem)
        
    rem[nd_mask] = -9999
    
    # Save the REM raster
    with rio.open(output, 'w', **masked_profile) as dst:
        dst.write(rem, indexes=1)
            
    # Optionally save the WSE raster
    if wse_path:
        with rio.open(wse_path, 'w', **masked_profile) as wse_dst:
            wse_dst.write(wse, indexes=1)
            
    print('\nREM created at {0}'.format(output))

def htab_1D(station, elevation, slope, D50, max_depth=10, breaks=None, save=None):
    """
    Calculate 1D Hydraulic Table for a cross-section.

    Parameters:
    station: array-like
        Station points along the cross-section.
    elevation: array-like
        Elevation points corresponding to the station points.
    slope: float
        Slope of the channel or cross-section.
    D50: float
        Median grain size for sediment.
    max_depth: float, optional
        Maximum depth to calculate the water surface elevation (WSE). Default is 10.
    breaks: array-like, optional
        Station points where there are changes in the cross-section, typically banks or breaks in the channel.
    save: str, optional
        Path to save the resulting hydraulic table (htab).

    Returns:
    DataFrame
        Hydraulic table with columns for depth, WSE, flow rate (Q), velocity (u), roughness (n),
        slope (S), cross-sectional area (A), wetted perimeter (P), hydraulic radius (R), and number of channels (n_chan).
    """

    # Convert station and elevation to numpy arrays
    station = np.array(station)
    elevation = np.array(elevation)
    
    # The lowest elevation in the cross-section (top of water ground level)
    twg = np.min(elevation)
    
    # Initialize an empty DataFrame for storing hydraulic data
    htab = pd.DataFrame([], columns=['d', 'wse', 'Q', 'u', 'n', 'S', 'A', 'P', 'R', 'n_chan'])
    
    # Iterate over depths from 0.1 to max_depth
    for d in np.linspace(0.1, max_depth, 25):
        wse = twg + d  # Water Surface Elevation for the current depth
        
        boo = elevation < wse  # Boolean array where elevation is below WSE
        
        if boo.all():
            # If all elevations are below WSE, the cross-section is unconfined, so skip this depth
            print('Error - Unconfined XS at {0}'.format(round(wse, 2)))
            continue
        
        # Find indices where the elevation crosses the WSE
        indices = np.nonzero(boo[1:] != boo[:-1])[0] + 1
        
        # Adjust indices to ensure proper splitting of station and elevation arrays
        indices = np.array([i-1 if np.where(indices==i)[0] % 2 == 0 else i for i in indices])
        
        if breaks:
            # Identify indices corresponding to station breaks (e.g., banks)
            chan_ind = np.array([np.argmin(abs(station - bound)) for bound in breaks])
            chan_ind = np.delete(chan_ind, elevation[chan_ind] > wse)
            indices = np.sort(np.unique(np.concatenate([indices, chan_ind]).flatten()))
        
        # Split station and elevation arrays at the calculated indices
        sta = np.split(station, indices)
        sta = [s for s in sta if (elevation[np.searchsorted(station, s)] < wse).any()]
        
        # Extend each array with the next station value after max(s), if applicable
        sta = [np.append(s, station[np.argwhere(station == max(s)) + 1][-1])
            if np.argwhere(station == max(s)) < len(station) - 1 else s
            for s in sta]
        
        # Store the index of the last element in each sub-array in end_sta
        end_sta = [np.argwhere(station == s[-1])[0][0] for s in sta]
        
        elev = np.split(elevation, indices)
        elev = [e for e in elev if (e < wse).any()]
        elev = [np.append(e, elevation[end_sta[i]]) for i, e in enumerate(elev)]
        
        # Create a list of channel segments as arrays of station and elevation
        channels = [np.array([sta[i], elev[i]]).astype(float) for i in range(len(sta))]
        channels = [chan for chan in channels if chan.size > 0]
        
        # Initialize lists to store calculated values
        aas, ps, rs, qs, us, ns = [], [], [], [], [], []
        
        # Calculate hydraulic properties for each channel segment
        for chan in channels:
            boo = chan[1] < wse  # Boolean array where channel elevation is below WSE
            
            # Interpolate the start and end points to match the WSE
            if wse < chan[1, 0]:
                chan[0, 0] = interp1d(chan[1, :2], chan[0, :2])(wse)
            if wse < chan[1, -1]:
                chan[0, -1] = interp1d(chan[1, -2:], chan[0, -2:])(wse)
                
            chan[1, [0, -1]] = wse  # Set the elevations at the start and end of the channel to WSE
            
            # Calculate cross-sectional area (A) using the trapezoidal rule
            a = trapezoid(wse - chan[1], chan[0])
            # Calculate wetted perimeter (P)
            p = np.sqrt(np.sum((chan.T[1:] - chan.T[:-1]) ** 2, -1)).sum()
            # Calculate hydraulic radius (R)
            r = a / p
            
            # Calculate friction factor (f) and Manning's roughness coefficient (n)
            f = 8 / (5.62 * np.log10(a / (max(chan[0]) - min(chan[0])) / D50) + 4) ** 2
            n = (1.49 / (8 * 32.17 * r * slope / f) ** 0.5) * r ** (2 / 3) * slope ** 0.5
            if n > 0.2:
                n = 0.2
            
            # Calculate flow velocity (u) and flow rate (q)
            u = (1.49 / n) * r ** (2 / 3) * slope ** 0.5
            q = a * u
            
            # Store calculated values for this channel segment
            aas.append(a)
            ps.append(p)
            rs.append(r * q)
            qs.append(q)
            us.append(u * q)
            ns.append(n * q)
            
        # Append calculated values for the current depth to the hydraulic table (htab)
        htab.loc[len(htab)] = [
            d, wse, sum(qs), sum(us) / sum(qs), sum(ns) / sum(qs),
            slope, sum(aas), sum(ps), sum(rs) / sum(qs), len(channels)
        ]
        
    # Save the hydraulic table to an Excel file or return it
    if save:
        htab.to_excel(save)
    else:
        return htab
    
def htab_2D(rem, extents, vcl, slope, roughness, channel_area=None, channel_roughness=0.035, save=None):
    """
    Calculate 2D Hydraulics based on terrain data.

    Parameters:
    rem: str
        Path to the raster elevation model (REM).
    extents: str or GeoDataFrame
        Extents within which calculations are performed.
    vcl: str or GeoDataFrame
        Vector contour lines for length calculations.
    slope: float or str
        The slope value or the name of the column in 'extents' containing slope values.
    roughness: float
        Manning's roughness coefficient for the floodplain.
    channel_area: str or GeoDataFrame, optional
        Area representing the active channel.
    channel_roughness: float, optional
        Manning's roughness coefficient for the channel. Default is 0.035.
    save: str, optional
        Path to save the resulting GeoDataFrame.

    Returns:
    GeoDataFrame
        Updated 'extents' with flow rate (Q) calculations if 'save' is not specified.
    """
    
    # Open the raster elevation model and read the data
    with rio.open(rem) as src:
        a = src.read(1)  # Elevation data array
        profile = src.profile.copy()  # Copy of the profile for CRS and transform information
        
    # Convert 'extents' to the same CRS as the raster
    if isinstance(extents, str):
        extents = gpd.read_file(extents).to_crs(profile['crs'])
    else:
        extents = extents.to_crs(profile['crs'])
        
    # Convert 'vcl' to the same CRS as the raster
    if isinstance(vcl, str):
        vcl = gpd.read_file(vcl).to_crs(profile['crs'])
    else:
        vcl = vcl.to_crs(profile['crs'])
        
    # Process the channel area if provided
    if channel_area is not None:
        if isinstance(channel_area, str):
            channel_area = gpd.read_file(channel_area).to_crs(profile['crs'])
        else:
            channel_area = channel_area.to_crs(profile['crs'])
        
        # Create a mask for the channel area
        channel_mask = rio.features.geometry_mask([channel_area.geometry.unary_union], a.shape, profile['transform'], invert=True)
    else:
        # If no channel area, create a mask filled with True (no channel)
        channel_mask = np.full(a.shape, True)
    
    # Calculate the slope of the surface from the elevation data
    px, py = np.gradient(a, profile['transform'][0])
    slope_surf = np.sqrt(px ** 2 + py ** 2)
    
    # Iterate through each row in 'extents' for flow rate calculations
    for i, row in extents.iterrows():
        # Create a mask for the current extent and remove the channel area from it
        mask = rio.features.geometry_mask([row.geometry], a.shape, profile['transform'], invert=True)
        mask = mask & ~channel_mask
        
        # Calculate the relative water surface elevation (WSE)
        rel_arr = row.WSE - a
        
        ### Floodplain calculations
        vol = rel_arr[mask].sum() * profile['transform'][0]**2  # Volume of water
        A = vol / vcl.geometry.length.sum()  # Cross-sectional area
        
        surface_area = np.sqrt(slope_surf[mask]**2 + profile['transform'][0]**2).sum()  # Wetted surface area
        P = surface_area / vcl.geometry.length.sum()  # Wetted perimeter
        
        # Determine the slope for the current extent
        if isinstance(slope, float):
            S = slope
        else:
            S = row[slope]
        
        # Calculate flow rate for the floodplain using Manning's equation
        FP_q = (1.49 / roughness) * A * (A / P)**(2/3) * S**0.5
        
        ### Active channel calculations (if channel area is provided)
        if channel_area is not None:
            vol = rel_arr[channel_mask].sum() * profile['transform'][0]**2  # Volume of water in the channel
            A = vol / vcl.geometry.length.sum()  # Cross-sectional area in the channel
            
            surface_area = np.sqrt(slope_surf[channel_mask]**2 + profile['transform'][0]**2).sum()  # Wetted surface area in the channel
            P = surface_area / vcl.geometry.length.sum()  # Wetted perimeter in the channel
            
            # Calculate flow rate for the channel using Manning's equation
            AC_q = (1.49 / roughness) * A * (A / P)**(2/3) * S**0.5
        else:
            AC_q = 0
        
        # Update the 'extents' GeoDataFrame with the total flow rate (Q)
        extents.loc[extents.index == i, 'Q'] = FP_q + AC_q
        
        # Print the result for the current water surface elevation
        print('Q for RWSE of {0} = {1}'.format(round(row.WSE, 2), round(FP_q + AC_q, 0)))
    
    # Save the updated 'extents' GeoDataFrame or return it
    if save:
        extents.to_file(save)
    else:
        return extents
    
def pebble_count(photo_file, obj_color='yellow', obj_size_mm=190.5, min_area=50):
    """
    Count pebbles in an image and calculate their sizes.

    Parameters:
    photo_file: str
        Path to the image file (JPEG).
    obj_color: str or tuple, optional
        The color of the reference object in the image. Default is 'yellow'.
    obj_size_mm: float, optional
        The size of the reference object in millimeters. Default is 190.5 mm.
    min_area: int, optional
        Minimum area for detecting a region in the image. Default is 50.

    Returns:
    None
        Saves the grain data and labeled image as files.
    """

    # Convert object color to RGBA format if it's a string
    if isinstance(obj_color, str):
        obj_color = colors.to_rgba(obj_color)
    
    # Load SAM model for automatic mask generation
    sam = sam_model_registry["vit_h"](checkpoint=r"Z:/Shared/W2r/Library/Python/sam_vit_h_4b8939.pth")
    mask_generator = SamAutomaticMaskGenerator(
        sam,
        points_per_side=64, 
        points_per_batch=32,
        pred_iou_thresh=0.7,
        stability_score_thresh=0.9,
        min_mask_region_area=min_area
    )
    
    # Load the image using Rasterio
    with rio.open(photo_file) as src:
        image = src.read()[:3, :, :]  # Read the first three bands (RGB)
        image = np.transpose(image, axes=[1, 2, 0])  # Transpose to get (height, width, channels)
        image = image.astype(np.uint8)
        
    # Generate masks for the image
    masks = mask_generator.generate(image)
    
    def merge_masks(anns):
        """
        Merge masks from the mask generator into a single mask image.
        
        Parameters:
        anns: list of dicts
            List of annotations with segmentation masks.

        Returns:
        img: ndarray
            Merged mask image.
        """
        sorted_anns = sorted(anns, key=lambda x: x['area'], reverse=True)
        img = np.zeros((sorted_anns[0]['segmentation'].shape[0], sorted_anns[0]['segmentation'].shape[1]))
        for i, ann in enumerate(sorted_anns):
            m = ann['segmentation']
            img[m] = i + 1
        return img
    
    # Merge masks and convert to integer format
    img = merge_masks(masks).astype(int)
    
    # Extract properties of the labeled regions
    props = regionprops_table(
        img, image,
        properties=('label', 'area', 'bbox', 'major_axis_length', 'minor_axis_length', 'intensity_mean')
    )
    grain_data = pd.DataFrame(props)
    grain_data[[col for col in grain_data.columns if 'intensity' in col]] /= 255
    
    # Calculate color closeness to the reference object color
    grain_data['closeness'] = np.abs(grain_data[[col for col in grain_data.columns if 'intensity' in col]] - obj_color[:3]).sum(axis=1)
    
    # Identify the reference object by finding the closest color match
    nb_label = grain_data.at[grain_data['closeness'].idxmin(), 'label']
    bbox = grain_data.loc[[grain_data['closeness'].idxmin()], ['bbox-1', 'bbox-0', 'bbox-3', 'bbox-2']].apply(lambda x: [i for i in x], axis=1).values[0]
    
    # Calculate pixel size in millimeters based on the reference object
    grain_data['pixel_size'] = obj_size_mm / grain_data.at[grain_data['closeness'].idxmin(), 'major_axis_length']
    grain_data['mm'] = grain_data['minor_axis_length'] * grain_data['pixel_size']
            
    # Annotate the reference object in the data
    grain_data['notes'] = np.nan
    grain_data.loc[grain_data.label == nb_label, 'notes'] = 'Reference Obj'
    
    # Remove unnecessary columns
    grain_data = grain_data.drop(columns=[col for col in grain_data.columns if 'intensity' in col] + ['closeness'])
    
    # Save grain data to a CSV file
    grain_data.to_csv(photo_file.replace('.jpg', '_grain_data.csv'), index=False)
    
    # Plot the labeled image and grain size distribution
    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(10, 8))
    ax[0].imshow(image)
    nb_img = np.where(img == nb_label, 1, np.nan)
    img = np.where(img == nb_label, 0, img)
    ax[0].imshow(img, alpha=0.5, cmap='tab20')
    ax[0].imshow(img == 0, alpha=0.5, cmap='bone_r')  
    ax[0].imshow(nb_img, cmap='autumn')
    gpd.GeoDataFrame([], geometry=[box(*bbox)]).plot(ax=ax[0], ec='r', fc='None', lw=1.2)
    
    ax[0].tick_params(
        axis='both',
        which='both',
        bottom=True,
        top=False,
        labelbottom=False,
        labelleft=False
    )
    
    # Plot grain size distribution
    bins = np.geomspace(0.1, 1000, 1000)
    counts, bin_edges = np.histogram(grain_data.mm, bins=bins, range=None, density=None, weights=grain_data.area)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax[1].plot(centers, np.cumsum(counts) / sum(counts) * 100, color='k')
    
    ax[1].set_xscale('log')
    ax[1].set_xlim(2, 1000)
        
    # Add shaded regions for grain size categories
    ax[1].axvspan(2, 4, color='k', alpha=0.25)
    ax[1].axvspan(8, 16, color='k', alpha=0.25)
    ax[1].axvspan(32, 64, color='k', alpha=0.25)
    ax[1].axvspan(128, 256, color='k', alpha=0.25)
    ax[1].axvspan(512, 1048, color='k', alpha=0.25)
    
    # Add labels for grain size categories
    ax[1].text(2.83, 99, 'Very Fine\nGravel', va='top', ha='center')
    ax[1].text(5.66, 99, 'Fine\nGravel', va='top', ha='center')
    ax[1].text(11.3, 99, 'Medium\nGravel', va='top', ha='center')
    ax[1].text(22.6, 99, 'Coarse\nGravel', va='top', ha='center')
    ax[1].text(45.3, 99, 'Very Coarse\nGravel', va='top', ha='center')
    ax[1].text(90.5, 99, 'Small\nCobbles', va='top', ha='center')
    ax[1].text(181, 99, 'Large\nCobbles', va='top', ha='center')
    ax[1].text(362, 99, 'Small\nBoulders', va='top', ha='center')
    ax[1].text(724, 99, 'Large\nBoulders', va='top', ha='center')
    
    # Interpolate grain size distribution for D84, D50, and D16
    gsd = interp1d(centers, np.cumsum(counts) / sum(counts) * 100)
    
    ax[1].text(0.99, 0.65, 'D84 = {0}'.format(gsd(84).round(1)), va='center', ha='right', transform=ax[1].transAxes)
    ax[1].text(0.99, 0.5, 'D50 = {0}'.format(gsd(50).round(1)), va='center', ha='right', transform=ax[1].transAxes)
    ax[1].text(0.99, 0.35, 'D16 = {0}'.format(gsd(16).round(1)), va='center', ha='right', transform=ax[1].transAxes)

    ax[1].set_xlabel('Grain Size (mm)')
    ax[1].set_ylabel('Percent Finer Than')
    ax[1].axhline(16, ls='--', c='k', lw=0.5)
    ax[1].axhline(50, ls='--', c='k', lw=0.5)
    ax[1].axhline(84, ls='--', c='k', lw=0.5)
    
    fig.tight_layout()
    
    # Save the labeled image
    fig.savefig(photo_file.replace('.jpg', '_grain_labels.png'), dpi=400)
    
    print('{0} pebbles counted for {1}'.format(len(grain_data), photo_file))

def delineate_trees(ch_file, output, canopy_floor=15, min_ht=50, max_ht=100, min_area=20, combine_dist=5):
    """
    Identify trees from a canopy height model (CHM) raster file and export the results to a shapefile.

    Parameters:
    ch_file : str
        Path to the canopy height model raster file.
    output : str
        Path to save the resulting tree locations shapefile.
    canopy_floor : int
        Height threshold to consider the ground level of the canopy (default 15).
    min_ht : int
        Minimum height to consider a detection as a tree (default 50).
    max_ht : int
        Maximum height to still consider the detection as a tree (default 100).
    min_area : int
        Minimum area for a tree detection to be considered valid (default 20).
    combine_dist : int
        Distance within which to combine tree points (default 5).

    Returns:
    None
    """

    with rio.open(ch_file) as src:
        ch = src.read(1)
        profile = src.profile.copy()
        nodata = profile['nodata']
        # Invert the CHM to facilitate pit detection as tree tops
        a = np.where(ch != nodata, ch * -1, np.nan)
        # Apply a Gaussian filter to smooth the CHM
        a = gaussian_filter(a, sigma=1)
        # Filter out areas below the canopy floor threshold
        a = np.where(a > -canopy_floor, np.nan, a)

    # Initialize the Grid and Raster objects
    chi_r = Raster(a)
    grid = Grid.from_raster(chi_r)
    fdir = grid.flowdir(chi_r)

    # Detect pits in the inverted CHM, which correspond to tree tops
    pits = grid.detect_pits(chi_r)
    idx = np.argwhere(pits)

    # Filter pits based on tree height thresholds
    pit_ht = ch[pits]
    tree_bool = (pit_ht > min_ht) & (pit_ht < max_ht)
    idx = idx[tree_bool]

    # Generate a GeoDataFrame of tree points
    points = gpd.GeoDataFrame([], geometry=[Point(rc[1], rc[0]) for rc in idx], crs=profile['crs'])
    points['Z'] = points.geometry.apply(lambda p: ch[int(p.y), int(p.x)])
    points['TID'] = points.index + 1

    # Function to filter and combine nearby points
    def filter_points(gdf, dist=combine_dist):
        for i, row in gdf.iterrows():
            buffer = row.geometry.buffer(dist)
            neighbors = gdf[gdf.geometry.intersects(buffer) & (gdf.index != i)]
            for j, neighbor in neighbors.iterrows():
                if neighbor['Z'] > row['Z']:
                    gdf.at[i, 'TID'] = neighbor['TID']
                    break
        return gdf

    points = filter_points(points)

    # Create the tree canopy areas
    trees = np.zeros(ch.shape, dtype=np.int32)
    for i, row in points.iterrows():
        catch = grid.catchment(x=row.geometry.x, y=row.geometry.y, fdir=fdir, xytype='index')
        if np.sum(catch) < min_area:
            continue
        trees[catch] = row['TID']

    mask = trees != 0
    results = ({'properties': {'TID': v}, 'geometry': s} for s, v in shapes(trees, mask=mask, transform=profile['transform']))

    # Convert results to a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features(list(results), crs=profile['crs'])

    # Save the GeoDataFrame to a file
    gdf.to_file(output)

def get_volume(terrain, df, target_elev='ELEVATION', method='cut', units='ft', save=None):
    """
    Calculate the volumes of cut or fill required for grading within polygons to a target elevation.

    Parameters:
    terrain : str
        Path to the terrain raster file.
    df : str or GeoDataFrame
        Path to the shapefile or a GeoDataFrame containing the polygons.
    target_elev : str or numeric
        The target elevation or the column in df specifying target elevations for each polygon.
    method : str
        Method to calculate volumes ('cut' or 'fill').
    units : str
        Units of measurement for volume output ('ft' for cubic yards, 'm' for cubic meters).
    save : str, optional
        File path where the resulting GeoDataFrame with volume calculations is saved, otherwise return a dataframe.

    Returns:
    df : GeoDataFrame
        GeoDataFrame with additional column 'Vol' for volumes, unless saved to file.
    """
    
    # Load the terrain raster
    with rio.open(terrain) as src:
        arr = src.read(1)
        transform = src.transform
        crs = src.crs

    # Load the polygon data
    if isinstance(df, str):
        df = gpd.read_file(df)
    df = df.to_crs(crs)
    
    # Assuming target_elev is a scalar numeric value...
    if np.isscalar(target_elev) and isinstance(target_elev, (int, float)):
        arr -= target_elev
        
        if method == 'cut':
            arr = np.where(arr < 0, 0, arr)
        elif method == 'fill':
            arr = np.where(arr > 0, 0, -arr)
        else:
            print("Method must be 'cut' or 'fill'.")
            return

        df['Vol'] = [np.sum(arr[rio.features.geometry_mask([geom], arr.shape, transform, invert=True)]) for geom in df.geometry]

    # Assuming target_elev is a column name in df...
    elif isinstance(target_elev, str): 
        df['Vol'] = np.nan

        for idx, row in df.iterrows():
            target = row[target_elev]
            arr_mod = arr - target
            
            if method == 'cut':
                arr_mod = np.where(arr_mod < 0, 0, arr_mod)
            elif method == 'fill':
                arr_mod = np.where(arr_mod > 0, 0, -arr_mod)
            else:
                print("Method must be 'cut' or 'fill'.")
                return

            mask = rio.features.geometry_mask([row.geometry], arr.shape, transform, invert=True)
            df.at[idx, 'Vol'] = np.sum(arr_mod[mask])
    else:
        print('Target elevation must be a scalar or column name.')
        return

    # Convert volumes to specified units
    if units == 'ft':
        df['Vol_KCY'] = df['Vol'] / transform[0] ** 3 / 1000
    elif units == 'm':
        df['Vol_KCY'] = df['Vol'] / transform[0] ** 3 * 1.308 / 1000
    else:
        print("Units must be 'ft' or 'm'.")
        return

    # Save or return the GeoDataFrame
    if save:
        df.to_file(save)
        print(f'Volumes calculated and saved to {save}')
    else:
        return df
