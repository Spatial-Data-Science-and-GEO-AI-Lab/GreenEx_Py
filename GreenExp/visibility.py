# Geographic Data Handling
import geopandas as gpd
import shapely.geometry as sg
import pyproj

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
    
    if not isinstance(viewing_dist, int) or (not viewing_dist > 0):
        raise TypeError("Please make sure that the viewing_dist argument is set to a positive integer")

    if not isinstance(sample_dist, (float, int)) or (not sample_dist > 0):
        raise TypeError("Please make sure that the sample_dist argument is set to a positive integer")
    
    if not isinstance(observer_height, (float, int)) or (not observer_height > 0):
        raise TypeError("Please make sure that the observer_height argument is set to a positive integer")

    with rasterio.open(dsm_raster_file) as src:
        dsm = src.read(1)
        dsm_crs = src.crs.to_epsg()
        dsm_bounds = src.bounds
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

    ### Step 3: 
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

    ### Step 4:
    poi[['GVI', 'nr_of_points']] = poi.apply(lambda row: pd.Series(f(mask, row, geom_type)), axis=1)
    print("Note: calculation of Viewshed GVI based on code by Johnny Huck and Labib Labib \nsource: https://github.com/jonnyhuck/green-visibility-index/blob/master/gvi.py \n")

    poi.drop('sampled_points', axis=1, inplace=True)

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
        poi.to_file(os.path.join(output_dir, f"{input_filename}_ViewshedGVI_added.gpkg"), driver="GPKG")
        sampled_points_gdf.to_file(os.path.join(output_dir, f"{input_filename}_ViewshedGVI_sampled_points.gpkg"), driver="GPKG")
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
    
    return np.mean(gvi_values).round(3), nr_of_points

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