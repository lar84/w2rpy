# w2rpy

**[GitHub](https://github.com/lar84/w2rpy)**

**[PyPI](https://pypi.org/project/w2rpy)**

To install or update:

```python
pip install w2rpy
pip install --upgrade w2rpy
```

PyTorch and Meta's SegmentAnything model must be installed separately to run the pebble_count function.
Download a pytorch wheel specific to your GPU [here](https://pytorch.org/get-started/locally/).
Download SegmentAnything model checkpoints [here](https://github.com/facebookresearch/segment-anything?tab=readme-ov-file#model-checkpoints).
And install the SegmentAnything python package:
```python
pip install https://github.com/facebookresearch/segment-anything/archive/master.zip
```

To create a REM and estimate channel extent:

```python
import w2rpy as w2r
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.plot import show

dem = '/dem.tif'
valley_centerline = '/VCL.shp'

xs = w2r.xs(valley_centerline,1500,100)

haws = r'C:/users/lrussell/downloads/haws.tif'
w2r.create_REM(dem,xs,haws)

channel = w2r.inundate(haws,[3],remove_holes=True,largest_only=True)


fig,ax = plt.subplots(1,1)

with rio.open(haws) as src:
    show(src,ax=ax,cmap='viridis',vmin=-1,vmax=10,alpha=0.75)
    
xs.plot(ax=ax,color='k',lw=0.5)
channel.plot(ax=ax,ec='b',fc='gray',alpha=0.5)
```
![image](https://github.com/lar84/w2rpy/blob/main/images/w2rpy_fig1.png)


To estimate a ratinge curve:

```python
import w2rpy as w2r
import geopandas as gpd
import matplotlib.pyplot as plt

import geopandas as gpd

# Get inundation polygons over a range of relative elevations
inundation_polygons = w2r.inundate(haws,range(3,11),largest_only=True)

# clip polygons to area of hydraulic calcs
aoi = gpd.read_file(r'C:/Users/lrussell/Downloads/aoi_clip.shp')
inundation_polygons.geometry = inundation_polygons.intersection(aoi.union_all())
channel.geometry = channel.intersection(aoi.union_all())

# run htab calcs
htab = w2r.htab_2D(haws,
            extents=inundation_polygons,
            vcl=valley_centerline,
            slope=0.01,
            roughness=0.1,
            channel_area=channel,
            channel_roughness=0.032)


fig,ax = plt.subplots(1,2)

inundation_polygons.plot(ax=ax[0],cmap='Blues_r',ec='k',alpha=0.1)

ax[1].plot(htab.Q,htab.WSE,color='k')
ax[1].scatter(htab.Q,htab.WSE,color='k')

ax[1].set_xlabel('Discharge (cfs)')
ax[1].set_ylabel('Relative WSE (ft)')

ax[1].grid()
```
![image](https://github.com/lar84/w2rpy/blob/main/images/w2rpy_fig2.png)

To run a photo-based pebble count:
```python
w2r.pebble_count(['/pebble_photo.png'],obj_size_mm=180)
```
![image](https://github.com/lar84/w2rpy/blob/main/images/w2rpy_fig3.png)

To run a habitat suitability (HSI) analysis:
```python
dep = r'/depth_raster.tif'
vel = r'/velocity_raster.tif'
w2r.HSI(dep,vel,'Adult Coho Spawning',r'/HSI_output.tif')
```
![image](https://github.com/lar84/w2rpy/blob/main/images/w2rpy_fig4.png)

To run a tree delineation analysis for trees between 50 and 100 feet tall:
```python
ch = '/canopy_height.tif'
output = r'/trees.shp'
w2r.delineate_trees(ch, output)
```
![image](https://github.com/lar84/w2rpy/blob/main/images/w2rpy_fig5.png)

## Documentation

**terrain(dem_file)**
- Returns Terrain object. Terrain object contains the following objects:
 	- DEM, flow accumulation, flow direction, pysheds Grid, and coordinate system
- Input options: 
 	- dem_file: file path to a DEM raster (.tif, .vrt, etc…)

**streamlines(terrain, pour_point, threshold=None, snap_threshold=None, save=None)**
- Returns a geodataframe of delineated stream lines.
- Input options:
 	- terrain: Terrain object
 	- pour_point: file path to or geodataframe of points at the bottom of the desired watershed. Will only use the first point.
 	- threshold: threshold for stream delineation in number of flow accumulated pixels. Defaults to number of pixels/100.
 	- snap_threshold: threshold for snapping pour point to the flow accumulation raster. Defaults to number of pixels/10.
 	- save: optional file path to save streamlines as shapefile rather than returning a geodataframe object.

**catchment(terrain, pour_point, snap_threshold=None, save=None)**
- Returns a geodataframe of delineated catchment.
- Input options:
 	- terrain: Terrain object
 	- pour_point: file path to or geodataframe of points at the bottom of the desired watershed. Will only use the first point.
 	- snap_threshold: threshold for snapping pour point to the flow accumulation raster. Defaults to number of pixels/10.
 	- save: optional file path to save catchment as shapefile rather than returning a geodataframe object.

**xs(cl, xs_length, spacing, save=None)**
- Returns a geodataframe of cross-sections.
- Input options:
 	- cl: centerline, or any lines, to generate cross-sections along.
 	- xs_length: length of cross-section from cl. Full XS length is xs_length*2
 	- spacing: distance between cross-sections on cl. 
 	- save: optional file path to save cross-sections as shapefile rather than returning a geodataframe object.

**points(lines, spacing, ep=True, save=None)**
- Returns a geodataframe of cross-sections.
- Input options:
 	- lines: centerline, or any lines, to generate cross-sections along.
 	- spacing: distance between points on lines. 
 	- ep: option to include end point on line. Either True or False.
 	- save: optional file path to save points as shapefile rather than returning a geodataframe object.

**inundate(raster, rel_wse_list, largest_only=False, invert=False, remove_holes=False, save=None)**
- Returns a geodataframe of polygons showing inundation extent for given relative WSEs. 
- Input options:
 	- raster: file path to DEM or REM to perform inundation on. 
 	- rel_wse_list: list for WSEs to create inundation polygons for. (Example: [1,3.5,7,10])
 	- largest_only: True/False to return only the largest inundated polygon. Helpful for channel delineation.
 	- invert: True/False to select the are below or above a given WSE.
 	- remove_holes: True/False option to remove all holes within inundation extent polygons.
 	- save: optional file path to save the inundation extents as shapefile rather than returning a geodataframe object.

**sample_raster(raster, points, buff=None, metric="min", multiband=False, crs=None)**
- Returns list of values sampled from raster at location of points.
- Input options:
 	- raster: file path to DEM or imagery to sample data. 
 	- points: geodatabase or file path of points to sample.
 	- buff: optional buffer distance to around each point to sample. 
 	- metric: suite of options to sample for including "min", "mean", "max", and "mode"
 	- multiband: True/False option to include samples from all bands, otherwise it only returns the first band. 
 	- crs: option to convert points to specific coordinate system, only use if raster has no set coordinate system. 

**zonal_stats(df,raster,metric="mean",band=1)**
- Returns list of values sampled from raster within extent of given polygons.
- Input options:
 	- raster: file path to DEM or imagery to sample data. 
 	- df: geodatabase or file path of polygons to sample.
 	- metric: suite of options including "min", "mean", "max", "sum", and "nonzero"
 	- band: integer value of which raster band to sample data from.  

**edit_raster(raster, output, crs=None, resample=None, resample_method="bilinear", clip=None, nodata=None)**
- All in one function to edit or warp rasters in various ways. Saves output to a new tif file. 
- Input options:
 	- raster: file path to DEM or imagery to warp/edit.
 	- output: file path to save new raster.
 	- crs: coordinate system to reproject raster to.
 	- resample: pixel size to resample raster to.
 	- resample_method: options for resampling ("bilinear","nearest","mean")
 	- clip: geodataframe or path to polygons to clip raster to. 
 	- nodata: will set the nodata value for the output raster. 

**merge_rasters(rasters, output, method="first", compression=None, nodata=None)**
- Merges rasters into a single tif file. This does not build a mosaic, be wary of resulting file sizes. 
- Input options:
 	- rasters: list of file paths to DEM or imagery to merge.
 	- output: file path to save new raster.
 	- method: options for merging rasters including "first", "last", "max" and "min". See rasterio documentation for more info. 
 	- compression: option to compress merged raster. Best to use are "JPEG" or "LZW".
 	- nodata: will set the nodata value for the output raster. 

**difference_rasters(r1, r2, output, match_affine="first", method="nearest")**
- Will match the extent and resolution of both rasters then difference. Output = r2 – r1.
- Input options:
 	- r1: file path to DEM or imagery.
 	- r2: file path to DEM or imagery.
 	- output: file path to save new raster.
 	- match_affine: If "first" then will warp r2 to match r1 extent and resolution. If "last" then will warp r1 to match r2 extent and resolution.
 	- method: method for raster transformation includes "nearest" and "bilinear"

**create_REM(dem, xs, output, sample_dist=3, smooth_window=5, buffer=1000, limits=[-50,50], method="min", vb=None, ret_xs=False, wse_path=None)**
- Will create a REM from a given DEM and set of cross-sections.
- Input options:
 	- dem: file path to DEM.
 	- xs: file path or geodataframe of cross-sections to use.
 	- output: file path to save new REM raster.
 	- sample_dist: Length along cross-sections to sample elevation data.
 	- smooth_window: number of cross-sections to average using rolling window.
 	- buffer: clips new raster to include buffer around total bounding box of cross-sections.
 	- limits: sets all data above or below limits to nodata.
 	- method: Option to use "min" or "mean" sampled elevation per cross-section
 	- vb: geodataframe or path to valley bottom polygon to clip cross-sections.
 	- ret_xs: will return XS if set to true, no need to use this.  
 	- wse_path: file path to save WSE/GGL surface raster.

**htab_1D(station, elevation, slope, D50, max_depth=10, breaks=None, save=None)**
- Returns dataframe of hydraulic table (htab) calcs for a given cross-section. Htab table returns either the sum (for discharge, area, and perimeter) or discharge-weighted average (for roughness, hydraulic radius, and velocity) across all flowing channels. 
- Input options:
 	- station: list , array, or series of station data.
 	- elevation: list , array, or series of elevation data.
 	- slope: dimensionless slope for hydraulic calcs. Constant across all flows/depths. 
 	- D50: representative grain size used for bathhurst calcs to vary roughness with depth. 
 	- max_depth: depth from thalweg to perform calcs from. Splits this range in 25 WSEs to build rating curve. 
 	- breaks: station values to split hydraulic calcs by. Helpful for separating channel from floodplain. 
 	- save: optional file path to save htab as excel rather than returning a dataframe object.

**htab_2D(rem,extents,vcl,slope,roughness,channel_area=None,channel_roughness=0.035,save=None)**
- Returns dataframe of reach-averaged hydraulic table (htab) calcs for a given 2D area. Htab table returns either the sum (for discharge, area, and perimeter) or discharge-weighted average (for roughness, hydraulic radius, and velocity) across all flowing channels. 
- Input options:
 	- rem: Relative elevation model raster file. Stages in the rating curve will be referencing this REM.
 	- extents: Shapefile of area to run hydraulic calcs. 
	- vcl: Shapefile of the valley center line, used to calculate reach length. 
 	- slope: dimensionless slope for hydraulic calcs. Constant across all flows/depths. 
 	- roughness: Averge or assumed floodplain roughness. Will apply to channel as well without channel_area provided. 
 	- channel_area: Optional shapefile of active channel area for separate hydraulic calcs from floodplain.
	- channel_roughness: Manning's roughness for channel, only used if channel_area is provided
 	- save: optional file path to save htab as excel rather than returning a dataframe object.

**pebble_count(photos, obj_size_mm, method='rapid')**
- Creates classified image with grain size distribution and CSV file of raw data in the same folder as your image.
- Input options:
	- photos: list of file paths to PNG/JPG of gravel bar with clear identifying object in the photo.
	- obj_size_mm: length of object along major axis in mm.
	- method: "rapid" for a faster Wolman pebble count or "detailed" for an area-weighted count of most visible grains 

**delineate_trees(ch_file, output, canopy_floor=10, min_ht=60, max_ht=120, min_area=20, combine_dist=5)**
- Saves a set of polygons representing delineated trees from a canopy height model.
- Input options:
	-  ch_file: Path to the canopy height model raster file.
    	-  output: Path to save the resulting tree locations shapefile.
    	-  canopy_floor: Height threshold to consider the ground level of the canopy (default 15).
    	-  min_ht: Minimum height to consider a detection as a tree (default 50).
    	-  max_ht: Maximum height to still consider the detection as a tree (default 100).
    	-  min_area: Minimum area for a tree detection to be considered valid (default 20).
    	-  combine_dist: Number of pixels within which to combine tree points (default 5).

**def HSI(dep_raster,vel_raster,curve,output)**
- Determines habitat suitability index given depth and velocity rasters. Values range from 0-1. Curves based on WDFW Instream Flow Guidelines.
- Input options:
	-  dep_raster: Path to the depth raster file.
	-  vel_raster: Path to the velocity raster file.
	-  curve: Species and lifestage curve, options are: 
		-"Adult Chinook Spawning Large River"
		-"Adult Chinook Spawning Small River"
		-"Juvenile Chinook Rearing"
		-"Adult Coho Spawning"
		-"Juvenile Coho Rearing"
		-"Adult Sockeye Spawning"
		-"Juvenile/Adult Rainbow Trout Rearing"
		-"Spring Chinook Holding"
		-"O. mykiss Juvenile"
    	-  output: Path to save the resulting HSI raster.
 
**Dcrit(shear_raster,output,tc=0.06,psed=162.3128,pwater=62.428)**
- Determines critical grain size (Dcrit) using Shields' equation. .
- Input options:
	-  shear_raster: Path to the shear stress (lbs/ft^2) raster file.
   	-  tc: dimensionless critical shear, use values between 0.045 and 0.06.
   	-  psed: density of water (lbs/ft^3)
   	-  pwater: density of water (lbs/ft^3)
    	-  output: Path to save the resulting Dcrit raster.

