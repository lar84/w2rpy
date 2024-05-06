# w2rpy
Water Resource Functions For Fluvial Systems

Install to any python env using:  
	pip install w2rpy  
To upgrade the package to stay up to date:  
	pip install --upgrade w2rpy  


terrain(dem_file)  
  •	Returns Terrain object. Terrain object contains the following rasters:  
    o	DEM, flow accumulation, flow direction, pysheds Grid, and coordinate system  
  •	Input options:   
    o	dem_file: file path to a DEM raster (.tif, .vrt, etc…)  
    
streamlines(terrain, pour_point, threshold=None, snap_threshold=None, save=None)  
  •	Returns a geodataframe of delineated stream lines.  
  •	Input options:  
    	o	terrain: Terrain object  
    o	pour_point: file path to or geodataframe of points at the bottom of the desired watershed. Will only use the first point.  
    o	threshold: threshold for stream delineation in number of flow accumulated pixels. Defaults to number of pixels/100.  
    o	snap_threshold: threshold for snapping pour point to the flow accumulation raster. Defaults to number of pixels/10.  
    o	save: optional file path to save streamlines as shapefile rather than returning a geodataframe object.  

catchment(terrain, pour_point, snap_threshold=None, save=None)  
  •	Returns a geodataframe of delineated catchment.  
  •	Input options:  
    o	terrain: Terrain object  
    o	pour_point: file path to or geodataframe of points at the bottom of the desired watershed. Will only use the first point.  
    o	snap_threshold: threshold for snapping pour point to the flow accumulation raster. Defaults to number of pixels/10.  
    o	save: optional file path to save catchment as shapefile rather than returning a geodataframe object.  

get_xs(cl, xs_length, spacing, save=None)  
  •	Returns a geodataframe of cross-sections.  
  •	Input options:  
    o	cl: centerline, or any lines, to generate cross-sections along.  
    o	xs_length: length of cross-section from cl. Full XS length is xs_length*2  
    o	spacing: distance between cross-sections on cl.   
    o	save: optional file path to save cross-sections as shapefile rather than returning a geodataframe object.  

get_points(lines, spacing, ep=True, save=None)  
  •	Returns a geodataframe of cross-sections.  
  •	Input options:  
    o	lines: centerline, or any lines, to generate cross-sections along.  
    o	spacing: distance between points on lines.   
    o	ep: option to include end point on line. Either True or False.  
    o	save: optional file path to save points as shapefile rather than returning a geodataframe object.  

inundate(raster, rel_wse_list, largest_only=False, invert=False, remove_holes=False, save=None)  
  •	Returns a geodataframe of polygons showing inundation extent for given relative WSEs.   
  •	Input options:  
    o	raster: file path to DEM or REM to perform inundation on.   
    o	rel_wse_list: list for WSEs to create inundation polygons for. (Example: [1,3.5,7,10])  
    o	largest_only: True/False to return only the largest inundated polygon. Helpful for channel delineation.  
    o	invert: True/False to select the are below or above a given WSE.  
    o	remove_holes: True/False option to remove all holes within inundation extent polygons.  
    o	save: optional file path to save the inundation extents as shapefile rather than returning a geodataframe object.  

sample_raster(raster, points, buff=None, metric='min', multiband=False, crs=None)  
  •	Returns list of values sampled from raster at location of points.  
  •	Input options:  
    o	raster: file path to DEM or imagery to sample data.   
    o	points: geodatabase or file path of points to sample.  
    o	buff: optional buffer distance to around each point to sample.   
    o	metric: suite of options to sample for including ‘min’, ’mean’, ’max’, and ‘mode’  
    o	multiband: True/False option to include samples from all bands, otherwise it only returns the first band.   
    o	crs: option to convert points to specific coordinate system, only use if raster has no set coordinate system.   

zonal_stats(df,raster,metric='mean',band=1)  
  •	Returns list of values sampled from raster within extent of given polygons.  
  •	Input options:  
    o	raster: file path to DEM or imagery to sample data.   
    o	df: geodatabase or file path of polygons to sample.  
    o	metric: suite of options including ‘min’, ’mean’, ’max’, ‘sum’, and ‘nonzero’  
    o	band: integer value of which raster band to sample data from.    

edit_raster(raster, output, crs=None, resample=None, resample_method='bilinear', clip=None, nodata=None)  
  •	All in one function to edit or warp rasters in various ways. Saves output to a new tif file.   
  •	Input options:  
    o	raster: file path to DEM or imagery to warp/edit.  
    o	output: file path to save new raster.  
    o	crs: coordinate system to reproject raster to.  
    o	Resample: pixel size to resample raster to.  
    o	resample_method: options for resampling (‘bilinear’,’nearest’,’mean’)  
    o	clip: geodataframe or path to polygons to clip raster to.   
    o	nodata: will set the nodata value for the output raster.   

merge_rasters(rasters, output, method='first', compression=None, nodata=None)  
  •	Merges rasters into a single tif file. This does not build a mosaic, be wary of resulting file sizes.   
  •	Input options:  
    o	rasters: list of file paths to DEM or imagery to merge.  
    o	output: file path to save new raster.  
    o	method: options for merging rasters including ‘first’, ‘last’, ‘max’ and ‘min’. See rasterio documentation for more info.   
    o	compression: option to compress merged raster. Best to use are ‘JPEG’ or ‘LZW’.  
    o	nodata: will set the nodata value for the output raster.   

difference_rasters(r1, r2, output, match_affine='first', method='nearest')  
  •	Will match the extent and resolution of both rasters then difference. Output = r2 – r1.  
  •	Input options:  
    o	r1: file path to DEM or imagery.  
    o	r2: file path to DEM or imagery.  
    o	output: file path to save new raster.  
    o	match_affine: If ‘first’ then will warp r2 to match r1 extent and resolution. If ‘last’ then will warp r1 to match r2 extent and resolution.  
    o	method: method for raster transformation includes ‘nearest’ and ‘bilinear’  

create_REM(dem, xs, output, sample_dist=3, smooth_window=5, buffer=1000, limits=[-50,50], method='min', vb=None, ret_xs=False, wse_path=None)  
  •	Will create a REM from a given DEM and set of cross-sections.  
  •	Input options:  
    o	dem: file path to DEM.  
    o	xs: file path or geodataframe of cross-sections to use.  
    o	output: file path to save new REM raster.  
    o	sample_dist: Length along cross-sections to sample elevation data.  
    o	smooth_window: number of cross-sections to average using rolling window.  
    o	buffer: clips new raster to include buffer around total bounding box of cross-sections.  
    o	limits: sets all data above or below limits to nodata.  
    o	method: Option to use ‘min’ or ‘mean’ sampled elevation per cross-section  
    o	vb: geodataframe or path to valley bottom polygon to clip cross-sections.  
    o	ret_xs: will return XS if set to true, no need to use this.    
    o	wse_path: file path to save WSE/GGL surface raster.  

htab(station, elevation, slope, D50, max_depth=10, breaks=None, save=None)  
  •	Returns dataframe of hydraulic table (htab) calcs for a given cross-section. Htab table returns either the sum (for discharge, area, and perimeter) or discharge-weighted average (for roughness, hydraulic radius, and velocity) across all flowing channels.   
  •	Input options:  
    o	station: list , array, or series of station data.  
    o	elevation: list , array, or series of elevation data.  
    o	output: file path to save new raster.  
    o	slope: dimensionless slope for hydraulic calcs. Constant across all flows/depths.   
    o	D50: representative grain size used for bathhurst calcs to vary roughness with depth.   
    o	max_depth: depth from thalweg to perform calcs from. Splits this range in 25 WSEs to build rating curve.   
    o	breaks: station values to split hydraulic calcs by. Helpful for separating channel from floodplain.   
    o	save: optional file path to save htab as excel rather than returning a dataframe object.  

