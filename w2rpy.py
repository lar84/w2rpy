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
from rasterio.features import shapes
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio import features
from rasterio.warp import Resampling
from rasterio.warp import reproject,calculate_default_transform
from rasterio.io import MemoryFile
from pysheds.grid import Grid
from pysheds.sview import Raster
from scipy.ndimage import gaussian_filter,zoom
from scipy.integrate import trapezoid
from scipy.interpolate import RBFInterpolator,interp1d
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path
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
    
    all_lines = gpd.GeoDataFrame([], columns=['WSID','FID','upstream','downstream','dist_to_pp','geometry'], crs=terrain.crs)
    for i,row in pour_point.interrows():
        pp = row.geometry
        
        xy = grid.snap_to_mask(terrain.acc > snap_threshold, [pp.x, pp.y], return_dist=False)
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
            lines.loc[line_fid, 'dist_to_pp'] = accumulated_distance
            downstream_fid = lines.loc[line_fid, 'downstream']
            
            if not downstream_fid.empty:
                next_distance = accumulated_distance + current_line.geometry.length
                calculate_distance_to_pour_point(downstream_fid.values[0], next_distance)
    
        # Start with the stream segment closest to the pour point
        pour_point_fid = lines[lines['downstream'].isnull()].index[0]
        calculate_distance_to_pour_point(pour_point_fid)
        
        lines['WSID'] = i
        all_lines = pd.concat([all_lines,lines])
    
    all_lines.crs = terrain.crs
    if save:
        # Save the stream network to a file
        all_lines.to_file(save)
        print(f'Stream network saved to {save}')
    else:
        # Return the stream network GeoDataFrame
        return all_lines
    
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
    
    
    all_cm = gpd.GeoDataFrame([], columns=['WSID','geometry'], crs=terrain.crs)
    for i,row in pour_point.interrows():
        pp = row.geometry

        # Extract coordinates of the pour point
        x0, y0 = pp.x, pp.y
    
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
        cm['WSID'] = i
        
        all_cm = pd.concat([all_cm,cm])
        
    all_cm.crs = terrain.crs
    
    # Save or return the GeoDataFrame
    if save:
        all_cm.to_file(save)
        print('Catchment delineated and saved to:', save)
    else:
        print('Catchment delineated')
        return all_cm
    
def xs(cl, xs_length, spacing, save=None):
    """
    Generate perpendicular cross-sectional lines along a centerline at specified intervals.

    Parameters:
    cl : str or GeoDataFrame
        Path to the centerline shapefile or a GeoDataFrame containing the centerline geometries.
    xs_length : float
        Total length of each cross-sectional line.
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
    xs_lines = gpd.GeoDataFrame(columns=['CSID', 'Distance', 'geometry'], crs=cl.crs)

    # Iterate through each line in the centerline GeoDataFrame
    for idx, row in cl.iterrows():
        geom = row.geometry
        num_xsecs = int(geom.length // spacing)
        
        for xs in range(1, num_xsecs + 1):
            distance = xs * spacing
            point = geom.interpolate(distance)  # Center of the cross-section
            
            # Tangent vector at the point (derivative of the curve)
            tangent = geom.interpolate(distance + 0.01).coords[0]
            center = point.coords[0]
            dx = tangent[0] - center[0]
            dy = tangent[1] - center[1]
            norm = np.hypot(dx, dy)
            unit_tangent = (dx / norm, dy / norm)

            # Perpendicular vector (rotated by 90 degrees)
            perp_vector = (-unit_tangent[1], unit_tangent[0])

            # Endpoints of the cross-sectional line
            end1 = (center[0] + xs_length/2 * perp_vector[0], center[1] + xs_length/2 * perp_vector[1])
            end2 = (center[0] - xs_length/2 * perp_vector[0], center[1] - xs_length/2 * perp_vector[1])
            line = LineString([end1, end2])

            # Append the line to the GeoDataFrame
            xs_lines.loc[len(xs_lines)] = [idx, distance, line]

    # Save or return the GeoDataFrame
    xs_lines.crs = cl.crs
    if save:
        xs_lines.to_file(save)
        print(f"Cross-sectional lines saved to {save}")
    else:
        return xs_lines

def points(lines, spacing, ep=True, save=None):
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
    def _points_ep(row):
        # Create points including the end point of the line
        intervals = np.arange(0, row.geometry.length, spacing).tolist() + [row.geometry.length]
        points = [row.geometry.interpolate(distance) for distance in intervals]
        return gpd.GeoDataFrame({
            'CSID': row.name,
            'Station': intervals,
            'geometry': points
        }, crs=lines.crs)

    # Define a function to create points excluding the endpoint
    def _points(row):
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
        lines['points'] = lines.apply(_points_ep, axis=1)
    else:
        lines['points'] = lines.apply(_points, axis=1)

    # Concatenate all points into a single GeoDataFrame
    points = pd.concat(lines['points'].tolist(), ignore_index=True)

    del lines['points']
    points.crs = lines.crs
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
    inundated.geometry = gpd.GeoSeries(inundated.geometry)
    inundated.crs = crs
    if save:
        inundated.to_file(save)
        print(f'Inundation data saved to {save}')
    else:
        return inundated

def edit_raster(raster, output, crs=None, resample=None, resample_method='bilinear', clip=None, nodata=None, match=None, match_transform=None, res=None):
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
        Pixel size to resample the raster to.
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
        og_profile = src.profile.copy()
        profile = og_profile.copy()
        
        if nodata is not None:
            profile['nodata'] = nodata
        profile['driver'] = 'GTiff'
        profile['BIGTIFF'] = 'YES'
        
        memfile1 = MemoryFile()
        memfile1.open(**profile).write(src.read())
        
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
                memfile1 = MemoryFile(memfile2.read())
                
        
        # Handling the resample operation
        if resample:
            with MemoryFile() as memfile2:
                with memfile1.open() as src:
                    transform = src.transform
                    
                    new_transform = rio.Affine(resample, transform.b, transform.c, transform.d, -resample, transform.f)
                    
                    new_width = int((src.bounds.right - src.bounds.left) / resample)
                    new_height = int((src.bounds.top - src.bounds.bottom) / resample)
                    
                    data = src.read(
                        out_shape=(src.count, new_height, new_width),
                        resampling=mdict[resample_method]
                    )

                    profile.update({
                        'transform': new_transform,
                        'width': new_width,
                        'height': new_height
                    })
                    
                    with memfile2.open(**profile) as dst:
                        dst.write(data)
                memfile1 = MemoryFile(memfile2.read())
        
        # Handling the reprojection operation
        if crs:
            with MemoryFile() as memfile2:
                with memfile1.open() as src:
                    if res is None:
                        res = src.res
                    transform, width, height = calculate_default_transform(
                        src.crs, crs, src.width, src.height, *src.bounds, resolution=res
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
                                                    
                memfile1 = MemoryFile(memfile2.read())
                
        # Handling the reprojection operation to match a raster
        if match:
            with rio.open(match) as src:
                target_profile = src.profile.copy()
                
            with MemoryFile() as memfile2:
                with memfile1.open() as src:
                    target_profile['dtype'] = src.profile['dtype']
                    target_profile['count'] = src.profile['count']
                    target_profile['nodata'] = src.profile['nodata']
                    
                    with memfile2.open(**target_profile) as dst_match:
                        # Reproject the raster to match the affine transform and CRS
                        for i in range(1, src.count + 1):
                            reproject(
                                source=rio.band(src, i),
                                destination=rio.band(dst_match, i),
                                src_transform=src.profile['transform'],
                                src_crs=src.crs,
                                dst_transform=target_profile['transform'],
                                dst_crs=target_profile['crs'],
                                resampling=mdict[resample_method])
                                                    
                memfile1 = MemoryFile(memfile2.read())

    # Saving the final raster
    with memfile1.open() as final_src:
        with rio.open(output, 'w', **final_src.profile) as final_dst:
            final_dst.write(final_src.read())

    del memfile1
    del memfile2
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
    profile['BIGTIFF'] = 'YES'
    
    # Apply compression if specified
    if compression:
        profile['compress'] = compression
    
    # Close all open datasets
    for src in open_datasets:
        src.close()
    
    # Write the merged raster to the output file
    with rio.open(output, 'w', **profile) as dst:
        dst.write(merged_array)
        
    print(f'Merged rasters at {output}')

def difference_rasters(r1, r2, output, match_affine='first', method='nearest'):
    """
    Subtract r1 from r2 and save the result as a new raster.

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
                
def zonal_stats(df, raster, metric='mean', profile=None, band=1, quant=None):
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
        'nonzero': np.count_nonzero,
        'quantile':np.nanquantile
    }
    
    # Function to apply the metric function to the raster within the geometry
    def run_func(func, geom, arr, trans, quant):
        # Create a mask where the geometry overlaps the raster
        mask = rio.features.geometry_mask([geom], arr.shape, trans, invert=True)
        
        # If all values in the masked area are NaN, return NaN
        if np.isnan(arr[mask]).all():
            return np.nan
        
        # Apply the metric function to the masked array
        if quant:
            val = func(arr[mask],quant)
        else:
            val = func(arr[mask])
        return val
    
    # Apply the run_func to each geometry in the GeoDataFrame
    stats = df.geometry.apply(lambda geom: run_func(met_dict[metric], geom, a, transform, quant))
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

def create_REM(dem, xs, output, Zcol=None, CSID=None, sample_dist=3, smooth_window=5, buffer=1000, limits=[-50, 50], method='min', vb=None, wse_path=None, save_xs=None):
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
    CSID: str, optional
        Column name for grouping cross-sections to smooth profile, use this to denote different streamlines/centerlines.
        If None, all cross-sections are assumed to be on the same centerline.
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
        profile = src.profile.copy()
    
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
        xsp = points(xs, sample_dist)
        xsp = xsp[~np.isnan(xsp.geometry.x)]
        xsp = xsp[~np.isnan(xsp.geometry.y)]
        
        xsp['Z'] = sample_raster(dem, xsp)
        xsp = xsp.dropna(subset='Z')
        
        if method == 'min':    
            xs['Z'] = xsp.groupby('CSID')['Z'].min()
        elif method == 'mean':
            xs['Z'] = xsp.groupby('CSID')['Z'].mean()
        
        if not CSID:
            xs['Z'] = xs['Z'].rolling(smooth_window, center=True, min_periods=1).mean()
        else:
            for i, group in xsp.groupby('CSID'):
                xs.loc[xs[CSID] == i, 'Z'] = xs.loc[xs[CSID] == i, 'Z'].rolling(smooth_window, center=True, min_periods=1).mean()

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
    int_points = int_points[~np.isnan(int_points.geometry.x)]
    int_points = int_points[~np.isnan(int_points.geometry.y)]
    int_points.index = range(len(int_points))

    # Clip the raster to the buffer of the cross-sections
    with rio.open(dem, 'r') as src:
        masked_array, masked_transform = mask(src, shapes=[box(*xs.total_bounds).buffer(buffer)], crop=True, nodata=-9999, indexes=1)
        masked_array = np.where(np.isnan(masked_array),-9999,masked_array)
        nd_mask = (masked_array == -9999)
        
        # Update the masked raster profile
        masked_profile = src.meta.copy()
        masked_profile.update({
            "driver": "GTiff",
            "height": masked_array.shape[0],
            "width": masked_array.shape[1],
            "transform": masked_transform,
            "nodata": -9999,
            'BIGTIFF': 'YES',
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
            
    # Optionally save the XS shapefile
    if save_xs:
        xs.to_file(save_xs)
            
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
    
def pebble_count(photos, obj_size_mm, method='rapid',model_path=None):
    """
    Count pebbles in an image and calculate their sizes.

    Parameters:
    photos: list
        List of image file paths (JPG or PNG).
    obj_size_mm: float, optional
        The major axis length of the reference object in millimeters. Default is 190.5 mm.
    method: str, optional
        Sample all pebbles in photo ('detailed'), or sample roughly 100 grains ('rapid')
    model_path: str, optional
        File path to SegmentAnything model checkpoint,
        if None then function uses a file path only available for w2r employees.
        https://github.com/facebookresearch/segment-anything?tab=readme-ov-file#model-checkpoints

    Returns:
    None
        Saves the grain data and labeled image as files.
    """
    def draw_polygon(image):

        # Store the polygon points
        polygon_points = []
    
        # Function to capture the polygon vertices
        def onselect(verts):
            nonlocal polygon_points
            polygon_points = verts
            print(f"Polygon coordinates: {polygon_points}")
    
        # Create the figure and axis for displaying the image
        fig, ax = plt.subplots()
        ax.imshow(image)
    
        # Allow the user to draw a polygon on the image
        poly_selector = PolygonSelector(ax, onselect)
    
        # Display the plot and let the user interact with it
        plt.show(block=True)
            
        path = Path(polygon_points)

        # Create a grid of coordinates corresponding to the image pixels
        height, width = image.shape[:2]
        xx, yy = np.meshgrid(np.arange(width), np.arange(height))
        coords = np.vstack((xx.flatten(), yy.flatten())).T
    
        # Check if each point is inside the polygon
        inside = path.contains_points(coords)
    
        # Reshape the result to match the image shape (height, width)
        mask_arr = inside.reshape((height, width)).astype('uint8')
        
        return mask_arr,polygon_points

    
    if method == 'rapid':
        params = [12,256,0.9,0.7,10,0.3]
    elif method == 'detailed':
        params = [64,256,0.9,0.9,10,0.3]
    else:
        print("Please choose the 'rapid' or 'detailed' maethod")
        return
    
    # Load SAM model for automatic mask generation
    model = None
    if not model_path:
        sam = sam_model_registry["vit_h"](checkpoint=r"Z:/Shared/W2r/Library/Python/sam_vit_h_4b8939.pth")
        model = 'vit_h'
    else:
        models = ['vit_b','vit_l','vit_h']
        for m in models:
            if m in model_path:
                model = m
                sam = sam_model_registry[m](checkpoint=model_path)
                break
    
    if not model:
        print('Model checkpoint needs to conatin version tag in file path: vit_b, vit_l, or vit_h.')
        return

    mask_generator = SamAutomaticMaskGenerator(
        sam,
        points_per_side=params[0],  # Reduced points
        points_per_batch=params[1],  # Increased batch size for speed
        pred_iou_thresh=params[2],  # Slightly lower threshold
        stability_score_thresh=params[3],  # Relaxed stability
        min_mask_region_area=params[4],  # Adjusted for smallest pebble size
        crop_overlap_ratio=params[5],  # Minimal overlap
    )
    
    for photo_file in photos:
    
        # Load the image using Rasterio
        with rio.open(photo_file) as src:
            image = src.read()[:3, :, :]  # Read the first three bands (RGB)
            image = np.transpose(image, axes=[1, 2, 0])  # Transpose to get (height, width, channels)
            scaling_factors = (0.5, 0.5, 1.0)
            image = zoom(image, zoom=scaling_factors, order=1)  # Scale down by 50%
            image = image.astype(np.uint8)
        
        ref_mask,ref_coords = draw_polygon(image)
        rect = Polygon(ref_coords).minimum_rotated_rectangle
        
        # calculate the length of each side of the minimum bounding rectangle
        ref_pixels = max([LineString((ref_coords[i], ref_coords[i+1])).length for i in range(len(ref_coords) - 1)])
        
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
        img[ref_mask.astype(bool)] = 0
        
        # Extract properties of the labeled regions
        props = regionprops_table(
            img, image,
            properties=('label', 'area', 'perimeter', 'centroid', 'orientation', 'major_axis_length', 'minor_axis_length')
        )
        grain_data = pd.DataFrame(props)
        grain_data = grain_data[grain_data.perimeter/grain_data.area<0.6]
        
        # Calculate pixel size in millimeters based on the reference object
        grain_data['pixel_size'] = obj_size_mm / ref_pixels
        grain_data['mm'] = grain_data['minor_axis_length'] * grain_data['pixel_size']
                
        # Save grain data to a CSV file
        grain_data.to_csv(photo_file.split('.')[0]+'_grain_labels.csv', index=False)
        
        # Plot the labeled image and grain size distribution
        fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], figsize=(8, 8))
        ax[0].imshow(image)
        ax[0].imshow(img, alpha=0.5, cmap='tab20')
        ax[0].imshow(img == 0, alpha=0.3, cmap='bone_r')  
        
        ref_mask = ref_mask.astype('float')
        ref_mask[ref_mask==0] = np.nan
        ax[0].imshow(ref_mask, cmap='autumn',alpha=0.5)
        gpd.GeoSeries([rect]).plot(ax=ax[0],ec='r',fc='None',zorder=10)
        
        for _, row in grain_data.iterrows():
            # Get properties for the fitted ellipse
            y0, x0 = row['centroid-0'], row['centroid-1']
            orientation = row['orientation'] + np.pi/2

            x1 = x0 + row['major_axis_length'] / 2 * np.cos(orientation)
            y1 = y0 - row['major_axis_length'] / 2 * np.sin(orientation)
            ax[0].plot([x0, x1], [y0, y1], 'k')

            x2 = x0 - row['minor_axis_length'] / 2 * np.sin(orientation)
            y2 = y0 - row['minor_axis_length'] / 2 * np.cos(orientation)
            ax[0].plot([x0, x2], [y0, y2], 'k')
            
            ax[0].plot(x0, y0, '.k', markersize=5)
        
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
        if method == 'rapid':
            counts, bin_edges = np.histogram(grain_data.mm, bins=bins, range=None, density=None)
        else:
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
        gsd = interp1d(np.cumsum(counts) / sum(counts) * 100, centers)
        
        ax[1].text(0.99, 0.65, 'D84 = {0}'.format(gsd(84).round(1)), va='center', ha='right', transform=ax[1].transAxes)
        ax[1].text(0.99, 0.5, 'D50 = {0}'.format(gsd(50).round(1)), va='center', ha='right', transform=ax[1].transAxes)
        ax[1].text(0.99, 0.35, 'D16 = {0}'.format(gsd(16).round(1)), va='center', ha='right', transform=ax[1].transAxes)
    
        ax[1].set_xlabel('Grain Size (mm)')
        ax[1].set_ylabel('Percent Finer Than')
        ax[1].axhline(16, ls='--', c='k', lw=0.5)
        ax[1].axhline(50, ls='--', c='k', lw=0.5)
        ax[1].axhline(84, ls='--', c='k', lw=0.5)
        
        ax[1].text(0.98, 0.02, 'n = {0}'.format(len(grain_data)), va='bottom', ha='right', transform=ax[1].transAxes)
        
        fig.tight_layout()
        
        # Save the labeled image
        fig.savefig(photo_file.split('.')[0]+'_grain_labels.png', dpi=400)
        
        print('{0} pebbles counted for {1}'.format(len(grain_data), photo_file))
        
def delineate_trees(ch_file, output, canopy_floor=15, min_ht=50, max_ht=100, min_area=20, combine_dist=5, veg_mask=None):
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
        
    if veg_mask is not None:            
        a = np.where(veg_mask,a,np.nan)

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

def HSI(dep_raster,vel_raster,curve,output):
    #time to define the preference curves for each salmon species and lifestage
    #all lists are indexed the same 0=depth, 1=depth preference value, 2=velocity, 3=velocity preference value
    ChnSpwnLR = [[0,0.55,1.05,1.55,5.05,10,30,35,99], #depth
                    [0,0,0.75,1,1,0,0,0,0], #depth pref
                    [0,0.55,0.75,1.55,3.55,4.95,6.55,7,99], #velocity
                    [0,0,0.79,1,1,0,0,0,0]] #velocity pref
    
    ChnSpwnSR = [[0,0.35,0.95,1.25,1.75,2.75,99],
                    [0,0,0.8,0.94,1,0.4,0.4],
                    [0,0.55,0.65,1.15,2.25,2.35,3.75,3.85,5,99],
                    [0,0,0.1,0.2,1,1,0.5,0.2,0,0]]
    
    ChnJuv = [[0,0.45,1.05,1.65,2.05,2.45,99],
                 [0,0,0.3,0.85,0.95,1,1],
                 [0,0.15,0.55,0.95,1.05,1.85,3.65,99],
                 [0.24,0.3,0.85,1,1,0.45,0,0]]
    
    CohoSpwn = [[0,0.15,0.55,0.85,1.15,1.55,1.95,2.75,99],
                   [0,0,0.65,1,1,0.9,0.53,0.35,0.35],
                   [0,0.45,1.25,1.45,4.25,5,99],
                   [0,0.53,1,1,0.62,0,0]]
    
    CohoJuv = [[0,0.1,0.25,1.55,2.5,3.25,3.9,4,99],
               [0,0,0.25,0.9,1,1,0.9,0.27,0.27],
               [0,0.15,0.3,0.45,0.6,1.2,2,99],
               [0.78,1,0.96,0.31,0.2,0.16,0,0]]
    
    SockSpwn = [[0,0.15,0.55,1.15,1.25,1.55,99],
                    [0,0,0.6,1,1,0.45,0.45],
                    [0,0.05,0.25,0.85,1.25,2.35,3.95,99],
                    [0,0,0.5,1,1,0.26,0,0]]
    
    RainbowTrout = [[0,0.55,1.55,2.25,2.6,2.75,3.4,4.75,99],
               [0,0,0.45,0.5,0.65,1,1,0.66,0.66],
               [0,0.85,1.75,2.65,3.7,5.25,99],
               [0.25,1,0.45,0.4,0.1,0,0]]
    
    SpringChnHold = [[0,0.8,2,6.5,99],
                    [0,0,0.1,1,1],
                    [0,2.4,3.8,4.8,6,99],
                    [1,1,0.8,0.2,0,0]]
    
    OmykissJuv = [[0,0.15,0.65,1.35,2.65,99],
                  [0,0,0.1,0.63,1,1],
                  [0,0.75,0.95,1.15,1.55,1.85,3.15,3.85,5,99],
                  [0.55,1,1,0.87,0.78,0.54,0.3,0.07,0,0]]
    
    labels = ['Adult Chinook Spawning Large River',
              'Adult Chinook Spawning Small River',
              'Juvenile Chinook Rearing',
              'Adult Coho Spawning',
              'Juvenile Coho Rearing',
              'Adult Sockeye Spawning',
              'Juvenile/Adult Rainbow Trout Rearing',
              'Spring Chinook Holding',
              'O. mykiss Juvenile']
    
    curves = [ChnSpwnLR,ChnSpwnSR,ChnJuv,CohoSpwn,CohoJuv,SockSpwn,RainbowTrout,SpringChnHold,OmykissJuv]
    curve_dict = dict(zip(labels,curves))
    
    dep_int = interp1d(curve_dict[curve][0],curve_dict[curve][1])
    vel_int = interp1d(curve_dict[curve][2],curve_dict[curve][3])
    
    with rio.open(dep_raster) as src:
        profile = src.profile.copy()

        dep = src.read(1)
        extent = np.where(dep==src.nodata,1,0)
        dep[extent] = 0
        dep = dep_int(dep)
        
    with rio.open(vel_raster) as src:
        vel = src.read(1)
        vel[extent] = 0
        vel = vel_int(vel)
    
    comp = np.sqrt(dep*vel)
    comp[extent] = -1
    
    profile['driver'] = 'GTIFF'
    profile['nodata'] = -1
    with rio.open(output,'w',**profile) as dst:
        dst.write(comp, indexes=1)
    
def Dcrit(shear_raster,output,tc=0.06,psed=162.3128,pwater=62.428):
    """
    Compute the critical diameter (Dcrit) from a shear stress raster and save the result as a new raster.
 
    Parameters:
    ----------
    shear_raster : str
        Path to the input raster file containing shear stress values.
    output : str
        Path to the output raster file where the Dcrit values will be saved.
    tc : float, optional
        Critical shear stress threshold. Default is 0.06.
    psed : float, optional
        Sediment density in the same units as pwater (e.g., lb/ft³ or kg/m³). Default is 162.3128.
    pwater : float, optional
        Water density in the same units as psed. Default is 62.428.
 
    Returns:
    -------
    None
        Writes the processed raster to the specified output file.
 
    Notes:
    ------
    - The function assumes the input raster has a valid nodata value.
    """
   
    with rio.open(shear_raster) as src:
        array = src.read()
        profile = src.profile.copy()
                
        mask = np.where(array==src.nodata,1,0)
        new_array = (array)/((psed - pwater)*tc)
        new_array[mask] = src.nodata
            
        profile['driver'] = 'GTIFF'
        with rio.open(output,'w',**profile) as dst:
            dst.write(new_array)
    
    
