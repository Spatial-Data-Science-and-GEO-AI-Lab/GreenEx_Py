# Geographic Data Handling
import geopandas as gpd
from geopandas.tools import sjoin
from vt2geojson.tools import vt_bytes_to_geojson  
import shapely.geometry as sg
from shapely.geometry import Point, box
from scipy.spatial import cKDTree 
import osmnx as ox  
import mercantile  

# Network Analysis and Routing
import osmnx as ox

# Data Manipulation and Analysis
import numpy as np
from numpy import zeros, column_stack
import pandas as pd

# Time and Date Handling
from time import time
from datetime import timedelta

# File and directory operations
import os

# Mathematical Functions
from math import exp, hypot

# Raster and Image Processing
from rasterio.transform import rowcol, xy
from skimage.draw import line, disk, circle_perimeter
import rasterio
from osgeo import gdal

# Image Processing and Analysis
from transformers import AutoImageProcessor, Mask2FormerForUniversalSegmentation  
from scipy.signal import find_peaks 
import torch  
from PIL import Image  
import requests  

# Libraries for Concurrency and File Manipulation
from concurrent.futures import ThreadPoolExecutor, as_completed 
import multiprocessing as mp 

# Progress Tracking
from tqdm import tqdm

##### MAIN FUNCTIONS
def get_viewshed_GVI(point_of_interest_file, greendata_raster_file, dtm_raster_file, dsm_raster_file, network_file=None, crs_epsg=None, 
                     polygon_type="neighbourhood", buffer_dist=None, viewing_dist=250, sample_dist=50, observer_height=1.7, write_to_file=True, 
                     output_dir=os.getcwd()):
    ### Step 1: Read and process user inputs, check conditions
    poi = gpd.read_file(point_of_interest_file)
    if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon'):
        geom_type = poi.iloc[0]['geometry'].geom_type
    else:
        raise ValueError("Please make sure all geometries are of 'Point' type or all geometries are of 'Polygon' type and re-run the function")

    # Make sure CRS of poi file is projected rather than geographic
    if not poi.crs.is_projected:
        if crs_epsg is None:
            print("Warning: The CRS of the PoI dataset is currently geographic, therefore it will now be projected to CRS with EPSG:3395")
            epsg = 3395
            poi.to_crs(f"EPSG:{epsg}", inplace=True)
        else:
            print(f"Warning: The CRS of the PoI dataset is currently geographic, therefore it will now be projected to EPSG:{crs_epsg} as specified")
            epsg = crs_epsg
            poi.to_crs(f"EPSG:{epsg}", inplace=True)
    else:
        epsg = poi.crs.to_epsg()

    # In case of house polygons, transform to centroids
    if geom_type == "Polygon":
        if polygon_type not in ["neighbourhood", "house"]:
            raise TypeError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
        if polygon_type == "house":
            print("Changing geometry type to Point by computing polygon centroids...")
            poi['geometry'] = poi['geometry'].centroid
            geom_type = poi.iloc[0]['geometry'].geom_type
            print("Done \n")

    # Make sure poi file contains ID columns to identify unique locations
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    # Validate user inputs
    if geom_type == "Point":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise TypeError("Please make sure that the buffer_dist argument is set to a positive integer")
    
    if not isinstance(viewing_dist, int) or (not viewing_dist > 0):
        raise TypeError("Please make sure that the viewing_dist argument is set to a positive integer")

    if not isinstance(sample_dist, (float, int)) or (not sample_dist > 0):
        raise TypeError("Please make sure that the sample_dist argument is set to a positive integer")
    
    if not isinstance(observer_height, (float, int)) or (not observer_height > 0):
        raise TypeError("Please make sure that the observer_height argument is set to a positive integer")

    # Read DSM, DTM and greenspace rasters
    with rasterio.open(dsm_raster_file) as src:
        dsm = src.read(1)
        dsm_crs = src.crs.to_epsg()
        dsm_bounds = src.bounds

    # Reproject if EPSG is not equal to CRS of poi file
    if not dsm_crs == epsg:
        print("Reprojecting the DSM file so that the CRS matches the CRS of the poi file...")
        gdal.Warp('/vsimem/reprojected.tif', dsm_raster_file, srcSRS=f"EPSG:{dsm_crs}", dstSRS=f"EPSG:{epsg}")
        with rasterio.open('/vsimem/reprojected.tif') as src:
            dsm = src.read(1)
            dsm_bounds = src.bounds
        gdal.Unlink('/vsimem/reprojected.tif')
        print("Done \n")

    # Make sure all points of interest are within or do at least intersect (in case of polygons) the DSM raster provided
    if not all(geom.within(sg.box(*dsm_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the DSM file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*dsm_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the DSM file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the DSM file provided, results will be based on intersecting part of polygons involved \n")
    
    with rasterio.open(dtm_raster_file) as src:
        dtm = src.read(1)
        dtm_crs = src.crs.to_epsg()
        dtm_bounds = src.bounds
    
    if not dtm_crs == epsg:
        print("Reprojecting the DTM file so that the CRS matches the CRS of the poi file...")
        gdal.Warp('/vsimem/reprojected.tif', dtm_raster_file, srcSRS=f"EPSG:{dtm_crs}", dstSRS=f"EPSG:{epsg}")
        with rasterio.open('/vsimem/reprojected.tif') as src:
            dtm = src.read(1)
            dtm_bounds = src.bounds
        gdal.Unlink('/vsimem/reprojected.tif')
        print("Done \n")
    
    # Make sure all points of interest are within or do at least intersect (in case of polygons) the DTM raster provided
    if not all(geom.within(sg.box(*dtm_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the DTM file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*dtm_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the DTM file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the DTM file provided, results will be based on intersecting part of polygons involved \n")
    
    meta = {
        'height': dtm.shape[0],
        'width': dtm.shape[1],
        'transform': src.transform,
        "driver": "GTiff",
        'count': 1,
        "crs": src.crs
    }

    with rasterio.open(greendata_raster_file) as src:
        green = src.read(1)
        green_crs = src.crs.to_epsg()
        green_bounds = src.bounds
    
    if not green_crs == epsg:
        print("Reprojecting the greenspace file so that the CRS matches the CRS of the poi file...")
        gdal.Warp('/vsimem/reprojected.tif', greendata_raster_file, srcSRS=f"EPSG:{green_crs}", dstSRS=f"EPSG:{epsg}")
        with rasterio.open('/vsimem/reprojected.tif') as src:
            green = src.read(1)
            green_bounds = src.bounds
        gdal.Unlink('/vsimem/reprojected.tif')
        print("Done \n")
        
    # Make sure all points of interest are within or do at least intersect (in case of polygons) the Greenspace raster provided
    if not all(geom.within(sg.box(*green_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the Greenspace file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*green_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the Greenspace file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the Greenspace file provided, results will be based on intersecting part of polygons involved \n")

    options = {
        'radius': viewing_dist, 
        'o_height': observer_height,  # 1.7 meters (typical height of a person)
        't_height': 0  # 0 meters (the target is the ground)
    }

    mask = {
        'meta': meta,
        'options': options,
        'dsm': dsm,
        'dtm': dtm,
        'green': green
    }

    ### Step 2: Retrieve network, use OSM if not provided
    if network_file is not None:
        if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
            raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
        elif network_file is not None and (os.path.splitext(network_file)[1] == ".gpkg"):
            graph_projected_edges = gpd.read_file(network_file, layer='edges')
        else: 
            graph_projected_edges = gpd.read_file(network_file)

        if not graph_projected_edges.crs.to_epsg() == epsg:
            print("Adjusting CRS of Network file to match with Point of Interest CRS...")
            graph_projected_edges.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Check if house locations are within network file provided
        bbox_network = graph_projected_edges.unary_union.envelope
        if not all(geom.within(bbox_network) for geom in poi['geometry']):
            raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")
    else:
        if buffer_dist is None:
            poi_polygon = sg.box(*poi.total_bounds)
        else:
            poi_polygon = sg.box(*poi.total_bounds).buffer(buffer_dist)
        polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") # Transform to 4326 for OSM
        wgs_polygon = polygon_gdf_wgs['geometry'].values[0] # Extract polygon in EPSG 4326

        print(f"Retrieving network within total bounds of {geom_type}(s) of interest, extended by the buffer_dist in case provided...")
        start_network_retrieval = time()
        network_graph = ox.graph_from_polygon(wgs_polygon, network_type='all')
        graph_projected = ox.project_graph(network_graph, to_crs=f"EPSG:{epsg}")
        graph_projected_edges = ox.graph_to_gdfs(graph_projected, nodes=False, edges=True)
        end_network_retrieval = time()
        elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")

    ### Step 3: Define sample points on road network for calculating GVI scores
    print("Computing sample points for roads within area of interest's network...")
    start_sample_points = time()
    poi['sampled_points'] = poi.apply(lambda row: get_network_sample_points(df_row=row, network_edges=graph_projected_edges, buffer_dist=buffer_dist, sample_dist=sample_dist), axis=1)
    # Explode the 'sampled_points' column
    sampled_points_exploded = poi.explode('sampled_points')[['id','sampled_points']].reset_index(drop=True).rename(columns={'sampled_points': 'geometry'})
    sampled_points_gdf = gpd.GeoDataFrame(sampled_points_exploded, crs=f'EPSG:{epsg}').reset_index(drop=True)
    end_sample_points = time()
    elapsed_sample_points = end_sample_points - start_sample_points
    print("Note: creation of sample points based on code by Ondrej Mlynarcik \nsource: https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/2.5D-GreenViewIndex-Netherlands/blob/main/sample_points_linestrings.ipynb")
    print(f"Done, running time: {str(timedelta(seconds=elapsed_sample_points))} \n")

    ### Step 4: Perform the Viewshed GVI calculation
    poi[['GVI', 'nr_of_points', 'GVI_list']] = poi.apply(lambda row: pd.Series(f(mask, row, geom_type)), axis=1)
    # Assign GVI scores to each sample road location
    sampled_points_gdf['GVI'] = poi['GVI_list'].explode().reset_index(drop=True)
    print("Note: calculation of Viewshed GVI based on code by Johnny Huck and Labib Labib \nsource: https://github.com/jonnyhuck/green-visibility-index/blob/master/gvi.py \n")

    # Drop irrelevant columns
    poi.drop(['sampled_points', 'GVI_list'], axis=1, inplace=True)

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
        poi.to_file(os.path.join(output_dir, f"{input_filename}_ViewshedGVI_added.gpkg"), driver="GPKG")
        sampled_points_gdf.to_file(os.path.join(output_dir, f"{input_filename}_ViewshedGVI_sampled_points.gpkg"), driver="GPKG")
        print("Done")

    return poi, sampled_points_gdf


def get_streetview_GVI(point_of_interest_file, access_token=None, crs_epsg=None, polygon_type="neighbourhood", buffer_dist=None,
                       sample_dist=50, network_file=None, write_to_file=True, output_dir=os.getcwd()):
    
    ### Step 1: Read and process user inputs, check conditions
    poi = gpd.read_file(point_of_interest_file)
    if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon'):
        geom_type = poi.iloc[0]['geometry'].geom_type
    else:
        raise ValueError("Please make sure all geometries are of 'Point' type or all geometries are of 'Polygon' type and re-run the function")

    if not poi.crs.is_projected:
        if crs_epsg is None:
            print("Warning: The CRS of the PoI dataset is currently geographic, therefore it will now be projected to CRS with EPSG:3395")
            epsg = 3395
            poi.to_crs(f"EPSG:{epsg}", inplace=True)
        else:
            print(f"Warning: The CRS of the PoI dataset is currently geographic, therefore it will now be projected to EPSG:{crs_epsg} as specified")
            epsg = crs_epsg
            poi.to_crs(f"EPSG:{epsg}", inplace=True)
    else:
        epsg = poi.crs.to_epsg()

    # In case of house polygons, transform to centroids
    if geom_type == "Polygon":
        if polygon_type not in ["neighbourhood", "house"]:
            raise TypeError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
        if polygon_type == "house":
            print("Changing geometry type to Point by computing polygon centroids...")
            poi['geometry'] = poi['geometry'].centroid
            geom_type = poi.iloc[0]['geometry'].geom_type
            print("Done \n")

    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    if geom_type == "Point":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise TypeError("Please make sure that the buffer_dist argument is set to a positive integer")

    if access_token is None:
        raise TypeError("Please make sure that an access token for Mapillary is provided")
    
    # Determine area of interest by taking bounding box of poi file, incl. buffer if specified
    if buffer_dist is None:
        poi_polygon = box(*poi.total_bounds)
        poi['buffer'] = poi['geometry']
    else:
        poi_polygon = box(*poi.total_bounds).buffer(buffer_dist)
        poi['buffer'] = poi['geometry'].buffer(buffer_dist)
    polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") # Transform to 4326 for OSM
    wgs_polygon = polygon_gdf_wgs['geometry'].values[0] # Extract polygon in EPSG 4326

    if network_file is not None:
        if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
            raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
        elif network_file is not None and (os.path.splitext(network_file)[1] == ".gpkg"):
            network_edges = gpd.read_file(network_file, layer='edges')
        else: 
            network_edges = gpd.read_file(network_file)

        if not network_edges.crs.to_epsg() == epsg:
            print("Adjusting CRS of Network file to match with Point of Interest CRS...")
            network_edges.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Check if house locations are within network file provided
        bbox_network = network_edges.unary_union.envelope
        if not all(geom.within(bbox_network) for geom in poi['geometry']):
            raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")
    else:
        print(f"Retrieving network within total bounds of {geom_type}(s) of interest, extended by the buffer_dist in case provided...")
        start_network_retrieval = time()
        network_edges = get_road_network_with_points(wgs_polygon, epsg=epsg)
        end_network_retrieval = time()
        elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")

    print("Computing sample points for roads within area of interest's network...")
    start_sample_points = time()
    road_points = select_points_on_road_network(network_edges)
    buffer_points = select_points_within_buffers(poi, road_points)
    end_sample_points = time()
    elapsed_sample_points = end_sample_points - start_sample_points
    print(f"Done, running time: {str(timedelta(seconds=elapsed_sample_points))} \n")
    
    print("Downloading StreetView images for road sample points...")
    start_images = time()
    features = get_features_on_points(buffer_points, access_token)
    gvi_per_point = download_images_for_points(features, access_token, epsg)
    end_images = time()
    elapsed_images = end_images - start_images
    print(f"Done, running time: {str(timedelta(seconds=elapsed_images))} \n")
    
    print("Calculating StreetView GVI score...")
    start_calc = time()
    poi, sampled_points_gdf = get_gvi_per_buffer(poi, gvi_per_point)
    end_calc = time()
    elapsed_calc = end_calc - start_calc
    print(f"Done, running time: {str(timedelta(seconds=elapsed_calc))} \n")

    print("Note: workflow for calculating Streetview GVI based on code by Ilse A. Vázquez Sánchez \nsource: https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/StreetView-NatureVisibility \n")
    
    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
        poi.to_file(os.path.join(output_dir, f"{input_filename}_StreetviewGVI_added.gpkg"), driver="GPKG")
        sampled_points_gdf.to_file(os.path.join(output_dir, f"{input_filename}_StreetviewGVI_sampled_points.gpkg"), driver="GPKG")
        print("Done")

    return poi, sampled_points_gdf


##### SUPPORTING FUNCTIONS
def coords2Array(a, x, y):
    """
    * convert between coords and array position
    *  returns row,col (y,x) as expected by rasterio
    """
    r, c = rowcol(a, x, y)
    return int(r), int(c)


def array2Coords(a, row, col):
    """
    * convert between array position and coords
    *  params are row,col (y,x) as expected by rasterio
    *  returns coords at the CENTRE of the cell
    """
    x, y = xy(a, row, col)
    return int(x), int(y)


def viewshed(r0, c0, radius_px, resolution, observerHeight, targetHeight, dsm_data, dtm_data, a):
    """
    * Use Bresenham's Circle / Midpoint algorithm to determine endpoints for viewshed
    """

    # create output array at the same dimensions as data for viewshed
    output = zeros(dtm_data.shape)

    # set the start location as visible automatically
    output[(r0, c0)] = 1

    # get pixels in the circle
    for r, c in column_stack(circle_perimeter(r0, c0, radius_px)):

        # calculate line of sight to each pixel
        output = lineOfSight(r0, c0, r, c, resolution, observerHeight, targetHeight, dsm_data, dtm_data, output)

    # return the resulting viewshed
    return output


def lineOfSight(r0, c0, r1, c1, observer_height, resolution, target_height, dsm_data, dtm_data, output):
    """
     * Runs a single ray-trace from one point to another point, returning a list of visible cells
    """

    # init variables for loop
    cur_dydx = 0 		  	# current dydx (base of object)
    max_dydx = 0 	  		# biggest dydx so far
    # top_dydx = 0 		    # current dydx (top of object)
    distance_travelled = 0  # how far we have travelled along the ray

    # get the viewer height
    height0 = dtm_data[(r0, c0)] + observer_height

    # get the pixels in the line (excluding the first one	)
    pixels = column_stack(line(r0, c0, r1, c1))[1:]

    # loop along the pixels in the line
    for r, c in pixels:

        # distance travelled so far
        distance_travelled = hypot(c0 - c, r0 - r)

        ''' comment this out as long as we use 0 as target offset '''
        ## set cell as visible if the height of the top of the object from the DTM > previous max
        # top_dydx = (dsm_data[(r, c)] - height0 + target_height) / distance_travelled
        # if (top_dydx >= max_dydx):
        # 	output[(r, c)] = 1
        #
        ## update max dydx the height of the base of the object on the DSM > previous max
        # cur_dydx = (dsm_data[(r, c)] - height0) / distance_travelled
        # if (cur_dydx > max_dydx):
        # 	max_dydx = cur_dydx

        # update max dydx the height of the base of the object on the DSM > previous max
        cur_dydx = (dsm_data[(r, c)] - height0) / (distance_travelled * resolution)
        if (cur_dydx > max_dydx):
            max_dydx = cur_dydx
            output[(r, c)] = 1

    # return updated output surface
    return output


def f(mask, df_row, geom_type):
    # create an output array at the same dimensions as data for output
    gvi = zeros((mask["meta"]["height"], mask["meta"]["width"]))

    # radius in pixels
    radius_px = int(mask["options"]["radius"] // mask['meta']['transform'][0])

    # build weighting mask
    weighting_mask = zeros((radius_px*2, radius_px*2))
    for r, c in column_stack(disk((radius_px, radius_px), radius_px, shape=weighting_mask.shape)):
        weighting_mask[(r, c)] = exp(-0.0003 * (hypot(radius_px - c, radius_px - r) * mask['meta']['transform'][0]))

    # get pixel references for aoi extents
    nr_of_points = len(df_row['sampled_points'])
    point_list = df_row['sampled_points']
    
    gvi_values = []
    for point in tqdm(point_list, desc=f'Calculating GVI for {geom_type} {df_row.id}'):
        r,c  = coords2Array(mask["meta"]["transform"], point.x, point.y)

        # call (weighted) viewshed
        output = viewshed(r, c, radius_px, 		# coords and radius in pixels
            mask['meta']['transform'][0],		# resolution of datasets
            mask["options"]["o_height"], 		# observer height
            mask["options"]["t_height"],		# target height
            mask["dsm"], 						# dsm dataset
            mask["dtm"],						# dtm dataset
            mask["meta"]["transform"])			# affine transform

        # extract the viewshed data from the output surface and apply weighting mask
        visible = output[r-radius_px:r+radius_px, c-radius_px:c+radius_px] * weighting_mask

        # multiply extract of (weighted) viewshed with extract of (weighted) green dataset
        visible_green = visible * (mask["green"][r-radius_px:r+radius_px, c-radius_px:c+radius_px] * weighting_mask)

        # get the ratio for greenness in the view
        gvi = visible_green.sum() / visible.sum()
        gvi_values.append(gvi)
    
    return np.mean(gvi_values).round(3), nr_of_points, gvi_values


# Get sample points within sub-network to calculate GVI
def get_network_sample_points(df_row, network_edges, buffer_dist, sample_dist):
    if buffer_dist is None:
        buffer_edges = network_edges[network_edges.intersects(df_row['geometry'])].reset_index(drop=True)
    else:
        buffer_edges = network_edges[network_edges.intersects(df_row['geometry'].buffer(buffer_dist))].reset_index(drop=True)
    
    edges_length_meters = buffer_edges.length
    edges_length_numeric = edges_length_meters.astype(float)
    
    sampled_points = []
    for i in range(len(buffer_edges)):
        if edges_length_numeric[i] < sample_dist:
            point = buffer_edges.geometry[i].centroid
            sampled_points.append(point)
        else:
            line = buffer_edges.geometry[i]
            num_points = int(line.length / sample_dist) + 1  # use this instead of distances to ensure more uniform sampling
            distances = np.linspace(0, line.length, num=num_points)
            points = [line.interpolate(distance) for distance in distances]
            sampled_points.extend(points)

    return sampled_points


def get_road_network_with_points(poi_polygon, epsg):
    # Get the road network within the bounding box
    G = ox.graph_from_polygon(poi_polygon, network_type='drive', simplify=True)

    #Project the graph from latitude-longitude coordinates to a local projection (in meters)
    G_proj = ox.project_graph(G, to_crs=f"EPSG:{epsg}")

    # Convert the projected graph to a GeoDataFrame
    edges = ox.graph_to_gdfs(G_proj, nodes=False, edges=True)

    return edges


# Get a list of points over the road map with a N distance between them
def select_points_on_road_network(network_edges, distance=50):
    # Initialize a list to store the points
    points = []

    # Loop through each road in the road network graph
    for road in network_edges.geometry:
        # Calculate the total length of the road
        road_length = road.length

        # Start at the beginning of the road
        current_position = 0

        # Loop through the road, adding points every 50 meters
        while current_position < road_length:
            # Get the point on the road at the current position
            current_point = road.interpolate(current_position)

            # Add the curent point to the list of points
            points.append(current_point)

            # Increment the position by the desired distance
            current_position += distance
    
    # Convert the list of points to a GeoDataFrame
    gdf_points = gpd.GeoDataFrame(geometry=points, crs=network_edges.crs)

    return gdf_points


def select_points_within_buffers(poi, road_points):
    points_within_buffers = sjoin(road_points, poi.set_geometry('buffer'), how='inner', predicate='within')

    # Get the unique points that fall within any buffer
    unique_points = points_within_buffers['geometry_left'].unique()

    # Create a new GeoDataFrame with the points that fall within any buffer
    return gpd.GeoDataFrame(geometry=[Point(p.x, p.y) for p in unique_points], crs=poi.crs)


def get_features_for_tile(tile, access_token):
    tile_url = f"https://tiles.mapillary.com/maps/vtp/mly1_public/2/{tile.z}/{tile.x}/{tile.y}?access_token={access_token}"
    response = requests.get(tile_url)
    result = vt_bytes_to_geojson(response.content, tile.x, tile.y, tile.z, layer="image")
    return [tile, result]


def get_features_on_points(buffer_points, access_token, zoom=14):
    # Transform the coordinate reference system to EPSG 4326
    buffer_points_wgs = buffer_points.copy(deep=True).to_crs(epsg=4326)

    # Add a new column to gdf_points that contains the tile coordinates for each point
    buffer_points['tile'] = [mercantile.tile(x, y, zoom) for x, y in zip(buffer_points_wgs.geometry.x, buffer_points_wgs.geometry.y)]

    # Group the points by their corresponding tiles
    groups = buffer_points.groupby('tile')

    # Download the tiles and extract the features for each group
    features = []
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []

        for tile, _ in groups:
            futures.append(executor.submit(get_features_for_tile, tile, access_token))
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading tiles"):
            result = future.result()
            features.append(result)

    pd_features = pd.DataFrame(features, columns=["tile", "features"])

    # Compute distances between each feature and all the points in gdf_points
    feature_points = pd.DataFrame(
        [(Point(f["geometry"]["coordinates"]), f) for row in pd_features["features"] for f in row["features"]],
        columns=["geometry", "feature"]
    )
    feature_tree = cKDTree(feature_points["geometry"].apply(lambda p: [p.x, p.y]).tolist())
    _, indices = feature_tree.query(buffer_points_wgs["geometry"].apply(lambda p: [p.x, p.y]).tolist())

    # Select the closest feature for each point
    buffer_points["feature"] = feature_points.loc[indices, "feature"].tolist()

    # Convert results to geodataframe
    buffer_points['tile'] = buffer_points['tile'].astype(str)
    
    return buffer_points


def get_models():
    processor = AutoImageProcessor.from_pretrained("facebook/mask2former-swin-large-cityscapes-semantic")
    model = Mask2FormerForUniversalSegmentation.from_pretrained("facebook/mask2former-swin-large-cityscapes-semantic")
    return processor, model


def run_length_encoding(in_array):
    image_array = np.asarray(in_array)
    length = len(image_array)
    if length == 0: 
        return (None, None, None)
    else:
        pairwise_unequal = image_array[1:] != image_array[:-1]
        change_points = np.append(np.where(pairwise_unequal), length - 1)   # must include last element posi
        run_lengths = np.diff(np.append(-1, change_points))       # run lengths
        return(run_lengths, image_array[change_points])


def get_road_pixels_per_column(prediction):
    road_pixels = prediction == 0.0 # The label for the roads is 0
    road_pixels_per_col = np.zeros(road_pixels.shape[1])
    
    for i in range(road_pixels.shape[1]):
        run_lengths, values = run_length_encoding(road_pixels[:,i])
        road_pixels_per_col[i] = run_lengths[values.nonzero()].max(initial=0)
    return road_pixels_per_col


def get_road_centres(prediction, distance=2000, prominence=100):
    road_pixels_per_col = get_road_pixels_per_column(prediction)
    peaks, _ = find_peaks(road_pixels_per_col, distance=distance, prominence=prominence)
    
    return peaks


def find_road_centre(segmentation):
	distance = int(2000 * segmentation.shape[1] // 5760)
	prominence = int(100 * segmentation.shape[0] // 2880)
	
	centres = get_road_centres(segmentation, distance=distance, prominence=prominence)
	
	return centres


def crop_panoramic_images(original_width, image, segmentation, road_centre):
    width, height = image.size

    # Find duplicated centres
    duplicated_centres = [centre - original_width for centre in road_centre if centre >= original_width]
            
    # Drop the duplicated centres
    road_centre = [centre for centre in road_centre if centre not in duplicated_centres]

    # Calculate dimensions and offsets
    w4 = int(width / 4) # 
    h4 = int(height / 4)
    hFor43 = int(w4 * 3 / 4)
    w98 = width + (w4 / 2)
    xrapneeded = int(width * 7 / 8)

    images = []
    pickles = []
    # Crop the panoramic image
    for centre in road_centre:
        # Wrapped all the way around
        if centre >= w98:
            xlo = int(centre - w4/2)
            cropped_image = image.crop((xlo, h4,  xlo+w4, h4 + hFor43))
            cropped_segmentation = segmentation[h4:h4+hFor43, xlo:xlo+w4]
        
        # Image requires assembly of two sides
        elif centre > xrapneeded:
            xlo = int(centre - (w4/2)) # horizontal_offset
            w4_p1 = width - xlo
            w4_p2 = w4 - w4_p1
            cropped_image_1 = image.crop((xlo, h4, xlo + w4_p1, h4 + hFor43))
            cropped_image_2 = image.crop((0, h4, w4_p2, h4 + hFor43))

            cropped_image = Image.new(image.mode, (w4, hFor43))
            cropped_image.paste(cropped_image_1, (0, 0))
            cropped_image.paste(cropped_image_2, (w4_p1, 0))

            cropped_segmentation_1 = segmentation[h4:h4+hFor43, xlo:xlo+w4_p1]
            cropped_segmentation_2 = segmentation[h4:h4+hFor43, 0:w4_p2]
            cropped_segmentation = torch.cat((cropped_segmentation_1, cropped_segmentation_2), dim=1)
        
        # Must paste together the two sides of the image
        elif centre < (w4 / 2):
            w4_p1 = int((w4 / 2) - centre)
            xhi = width - w4_p1
            w4_p2 = w4 - w4_p1

            cropped_image_1 = image.crop((xhi, h4, xhi + w4_p1, h4 + hFor43))
            cropped_image_2 = image.crop((0, h4, w4_p2, h4 + hFor43))

            cropped_image = Image.new(image.mode, (w4, hFor43))
            cropped_image.paste(cropped_image_1, (0, 0))
            cropped_image.paste(cropped_image_2, (w4_p1, 0))

            cropped_segmentation_1 = segmentation[h4:h4+hFor43, xhi:xhi+w4_p1]
            cropped_segmentation_2 = segmentation[h4:h4+hFor43, 0:w4_p2]
            cropped_segmentation = torch.cat((cropped_segmentation_1, cropped_segmentation_2), dim=1)
            
        # Straightforward crop
        else:
            xlo = int(centre - w4/2)
            cropped_image = image.crop((xlo, h4, xlo + w4, h4 + hFor43))
            cropped_segmentation = segmentation[h4:h4+hFor43, xlo:xlo+w4]
        
        images.append(cropped_image)
        pickles.append(cropped_segmentation)

    return images, pickles


def segment_images(image, processor, model):
    inputs = processor(images=image, return_tensors="pt")
    
    # Forward pass
    with torch.no_grad():
        outputs = model(**inputs)
    
    # You can pass them to processor for postprocessing
    segmentation = processor.post_process_semantic_segmentation(outputs, target_sizes=[image.size[::-1]])[0]

    return segmentation


def get_GVI(segmentations):
    green_percentage = 0
    for segment in segmentations:
        total_pixels = segment.numel()
        vegetation_pixels = (segment == 8).sum().item()
        green_percentage += vegetation_pixels / total_pixels
    
    return green_percentage / len(segmentations)


def process_images(image_url, is_panoramic, processor, model):
    try:
        image = Image.open(requests.get(image_url, stream=True).raw)

        if is_panoramic:
            # Get the size of the image
            width, height = image.size

            # Crop the bottom 20% of the image to cut the band on the bottom of the panoramic image
            bottom_crop = int(height * 0.2)
            image = image.crop((0, 0, width, height - bottom_crop))

        # Image segmentation
        segmentation = segment_images(image, processor, model)

        if is_panoramic:
            # Create a widened panorama by wrapping the first 25% of the image onto the right edge
            width, height = image.size
            w4 = int(0.25 * width)

            segmentation_25 = segmentation[:, :w4]
            # Concatenate the tensors along the first dimension (rows)
            segmentation_road = torch.cat((segmentation, segmentation_25), dim=1)
        else:
            segmentation_road = segmentation
        
        # Find roads to determine if the image is suitable for the analysis or not AND crop the panoramic images
        road_centre = find_road_centre(segmentation_road)

        if len(road_centre) > 0:
            if is_panoramic:
                images, pickles = crop_panoramic_images(width, image, segmentation_road, road_centre)
            else:
                images = [image]
                pickles = [segmentation]
        
            # Now we can get the Green View Index
            GVI = get_GVI(pickles)
            return [GVI, is_panoramic, False, False]
        else:
            # There are not road centres, so the image is unusable
            return [None, None, True, False]
    except:
        return [None, None, True, True]


def download_image(geometry, image_metadata, access_token, processor, model):
    header = {'Authorization': 'OAuth {}'.format(access_token)}

    image_id = image_metadata["properties"]["id"]
    is_panoramic = image_metadata["properties"]["is_pano"]
    
    url = 'https://graph.mapillary.com/{}?fields=thumb_original_url'.format(image_id)
    response = requests.get(url, headers=header)
    data = response.json()
    image_url = data["thumb_original_url"]

    result = process_images(image_url, is_panoramic, processor, model)
    result.insert(0, geometry)

    return result


def process_data(index, data_part, processor, model, access_token, max_workers):
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for _, row in data_part.iterrows():
            feature = row["feature"]
            geometry = row["geometry"]
            futures.append(executor.submit(download_image, geometry, feature, access_token, processor, model))
    
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Downloading images (Process {index})"):
            image_result = future.result()
            results.append(image_result)
        return results
    
'''
def download_images_for_points(gdf, access_token, epsg, max_workers=1):
    processor, model = get_models()
    
    gdf_wgs = gdf.copy(deep=True).to_crs("EPSG:4326")

    images_results = []

    # Split the dataset into parts
    num_processes = mp.cpu_count() # Get the number of CPU cores
    data_parts = np.array_split(gdf_wgs, num_processes) # Split the dataset
    
    with mp.get_context("spawn").Pool(processes=num_processes) as pool:
        # Apply the function to each part of the dataset using multiprocessing
        results = pool.starmap(process_data, [(index, data_part, processor, model, access_token, max_workers) for index, data_part in enumerate(data_parts)])

        # Combine the results from all parts
        images_results = [result for part_result in results for result in part_result]

        # Close the pool to release resources
        pool.close()
        pool.join()

    images_results_gdf = gpd.GeoDataFrame(images_results, columns=["geometry", "GVI", "is_panoramic", "missing", "error"], crs="EPSG:4326").to_crs(f"EPSG:{epsg}")
    return images_results_gdf
'''

def download_images_for_points(gdf, access_token, epsg, max_workers=1):
    processor, model = get_models()

    gdf_wgs = gdf.copy(deep=True).to_crs("EPSG:4326")

    images_results = []

    # Split the dataset into parts
    num_processes = 1
    data_parts = np.array_split(gdf_wgs, num_processes)

    # Sequential execution without multiprocessing
    for index, data_part in enumerate(data_parts):
        part_results = process_data(index, data_part, processor, model, access_token, max_workers)
        images_results.extend(part_results)

    images_results_gdf = gpd.GeoDataFrame(images_results, columns=["geometry", "GVI", "is_panoramic", "missing", "error"], crs="EPSG:4326").to_crs(f"EPSG:{epsg}")
    return images_results_gdf


def get_gvi_per_buffer(buffered_points, gvi_per_point):
    joined = gpd.sjoin(gvi_per_point, buffered_points.set_geometry('buffer'), how='inner', predicate='within').drop('index_right', axis=1)
    
    # Group the points by buffer
    grouped = joined.groupby('id', group_keys=True)
    # Convert 'grouped' to a DataFrame
    grouped_df = grouped.apply(lambda x: x.reset_index(drop=True))
    grouped_df = grouped_df[["geometry_left", "GVI", "is_panoramic", "missing"]].reset_index()
    # Convert grouped_df to a GeoDataFrame
    grouped_gdf = gpd.GeoDataFrame(grouped_df, geometry='geometry_left').rename(columns={'geometry_left':'geometry'}).drop('level_1', axis=1)
    grouped_gdf = grouped_gdf.set_geometry('geometry')
    # Calculate the average 'gvi' for each group
    avg_gvi = np.round(grouped['GVI'].mean(),3).reset_index()
    nr_of_points = grouped['GVI'].count().reset_index(name='nr_of_points')
    # Merge with the buffered_points dataframe to get the buffer geometries
    result = avg_gvi.merge(buffered_points, left_on='id', right_on='id')
    result = result.merge(nr_of_points, on='id')
    # Convert the result to a GeoDataFrame
    result = gpd.GeoDataFrame(result[['id', 'geometry', 'GVI', 'nr_of_points']])

    return result, grouped_gdf