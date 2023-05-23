# Data manipulation and analysis
import numpy as np
import pandas as pd

# File and directory operations
import os

# Geospatial data processing and analysis
import geopandas as gpd
import osmnx as ox
import networkx as nx
import momepy
import rioxarray
import xrspatial
from rasterio.enums import Resampling
import pyproj
import shapely.geometry as sg
from shapely.ops import transform

# Geospatial data access and catalogs
import pystac_client
import planetary_computer
import odc.stac

# Date and time manipulation
from datetime import datetime

# Progress tracking
from tqdm import tqdm

##### MAIN FUNCTIONS
def get_mean_NDVI(point_of_interest_file, ndvi_raster_file=None, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None, network_type=None,
                  trip_time=None, travel_speed=None, output_dir=os.getcwd(), year=datetime.now().year):
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

    epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326") # Transformer to use planetary computer and OSM

    # Make sure poi dataframe contains ID column
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    # Retrieve NDVI raster, use planetary computer if not provided by user 
    if ndvi_raster_file is None:
        print("Retrieving NDVI raster through planetary computer...")
        aoi_source = sg.box(*poi.total_bounds).buffer(buffer_dist*1.1)
        bounding_box_pc = transform(epsg_transformer.transform, aoi_source).bounds  # transform CRS to comply with planetary computer requirements
        bounding_box_pc = [bounding_box_pc[1], bounding_box_pc[0], bounding_box_pc[3], bounding_box_pc[2]] # Swap coords order to match with planetary computer format

        # Query planetary computer
        catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
        # Obtain Area of Interest
        time_of_interest = f"{year}-01-01/{year}-12-30" 
        # Search Data
        search = catalog.search(collections=["sentinel-2-l2a"],
                                bbox=bounding_box_pc,
                                datetime=time_of_interest,
                                query={"eo:cloud_cover": {"lt": 20}})
        # Obtain Data
        items = search.item_collection()
        # Find the item with the lowest cloud cover
        selected_item = min(items, key=lambda item: item.properties["eo:cloud_cover"])
        # Obtain Bands of Interest
        selected_item_data = odc.stac.stac_load([selected_item], bands = ['red', 'green', 'blue', 'nir'], bbox = bounding_box_pc).isel(time=0)
        # Calculate NDVI values
        ndvi = xrspatial.multispectral.ndvi(selected_item_data['nir'], selected_item_data['red'])
        # Reproject to original poi CRS
        ndvi_src = ndvi.rio.reproject(f"EPSG:{epsg}", resampling= Resampling.nearest, nodata=np.nan)
        print("Done \n")
    else:
        ndvi_src = rioxarray.open_rasterio(ndvi_raster_file)
        if not ndvi_src.rio.crs.to_epsg() == epsg:
            print("Adjusting CRS of NDVI file to match with Point of Interest CRS...")
            ndvi_src.rio.write_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Make sure all points of interest are within or do at least intersect (in case of polygons) the NDVI raster provided
        if not all(geom.within(sg.box(*ndvi_src.rio.bounds())) for geom in poi['geometry']):
            if geom_type == "Point":
                raise ValueError("Not all points of interest are within the NDVI file provided, please make sure they are and re-run the function")
            else:
                if not all(geom.intersects(sg.box(*ndvi_src.rio.bounds())) for geom in poi['geometry']):
                    raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the NDVI file provided, please make sure they are/do and re-run the function")
                else:
                    print("Warning: Not all polygons of interest are completely within the area covered by the NDVI file provided, results will be based on intersecting part of polygons involved \n") 

    ### Step 2: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        if geom_type == "Polygon":
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
        else:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")
    else:
        if buffer_type not in ["euclidian", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")

        if buffer_type == "euclidian":
            if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                raise TypeError("Please make sure that the buffer distance is set as a positive integer")             

            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:            
            if network_file is None:
                if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                    raise TypeError("Please make sure that the buffer distance is set as a positive integer")

                if network_type not in ["walk", "bike", "drive", "all"]:
                    raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                
                if not isinstance(travel_speed, int) or (not travel_speed > 0):
                    raise TypeError("Please make sure that the travel speed is set as a positive integer")

                if not isinstance(trip_time, int) or (not trip_time > 0):
                    raise TypeError("Please make sure that the trip time is set as a positive integer")

                if geom_type == "Polygon":
                    print("Changing geometry type to Point by computing polygon centroids so that network can be retrieved...")
                    poi['geometry'] = poi['geometry'].centroid
                    print("Done \n")             
        
                meters_per_minute = travel_speed * 1000 / 60  # km per hour to m per minute

                aoi_geometry = []
                for geom in tqdm(poi['geometry'], desc = 'Retrieving network for point(s) of interest'):
                    latlon = epsg_transformer.transform(geom.x, geom.y) # Transform point geometry into latlon for OSMnx
                    graph = ox.graph_from_point(latlon, network_type=network_type, dist=buffer_dist) # Retrieve street network for desired network type and buffer distance surrounding poi
                    graph = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") # Project street network graph back to original poi CRS
                    center_node = ox.distance.nearest_nodes(graph, geom.x, geom.y) # Find node which is closest to point location as base for next steps
                    for _, _, _, data in graph.edges(data=True, keys=True): # Calculate the time it takes to cover each edge's distance
                        data["time"] = data["length"] / meters_per_minute
                    isochrone_poly = make_iso_poly(graph, center_node=center_node, trip_time=trip_time) # See separate function for line by line explanation
                    aoi_geometry.append(isochrone_poly)

                aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
                print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")   
            else:
                if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
                    raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
                elif os.path.splitext(network_file)[1] == ".gpkg":
                    network = gpd.read_file(network_file, layer='edges')
                else:
                    network = gpd.read_file(network_file)

                if not network.crs.to_epsg() == epsg:
                    print("Adjusting CRS of Network file to match with Point of Interest CRS...")
                    network.to_crs(f'EPSG:{epsg}', inplace=True)
                    print("Done \n")

                # Create bounding box for network file
                bbox_network = network.unary_union.convex_hull

                if not all(geom.within(bbox_network) for geom in poi['geometry']):
                    raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")

                aoi_gdf = gpd.GeoDataFrame(geometry=[bbox_network], crs=f'EPSG:{epsg}')            
    
    ### Step 3: Calculate mean NDVI values and write results to file
    print("Calculating mean NDVI values...")
    if not all(geom.within(sg.box(*ndvi_src.rio.bounds())) for geom in aoi_gdf['geometry']):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the NDVI raster, note that results will be based on the intersecting part of the buffer zone")
    poi['mean_NDVI'] = aoi_gdf.apply(lambda row: ndvi_src.rio.clip([row.geometry]).clip(min=0).mean().values.round(3), axis=1)
    print("Done \n")

    print("Writing results to new geopackage file in specified directory...")
    input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
    poi.to_file(os.path.join(output_dir, f"{input_filename}_ndvi_added.gpkg"), driver="GPKG")
    print("Done")
    
    return poi

def get_landcover_percentages(point_of_interest_file, landcover_raster_file, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None, 
                              network_type=None, trip_time=None, travel_speed=None, output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
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

    # Make sure poi dataframe contains ID column
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    landcover_src = rioxarray.open_rasterio(landcover_raster_file)
    if not landcover_src.rio.crs.to_epsg() == epsg:
        print("Adjusting CRS of Land Cover file to match with Point of Interest CRS...")
        landcover_src.rio.write_crs(f'EPSG:{epsg}', inplace=True)
        print("Done \n")

    # Make sure all points of interest are within or do at least intersect (in case of polygons) the NDVI raster provided
    if not all(geom.within(sg.box(*landcover_src.rio.bounds())) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the Land Cover file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*landcover_src.rio.bounds())) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the Land Cover file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the Land Cover file provided, results will be based on intersecting part of polygons involved \n")

    ### Step 2: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        if geom_type == "Polygon":
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
        else:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")
    else:
        if buffer_type not in ["euclidian", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")

        if buffer_type == "euclidian":
            if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                raise TypeError("Please make sure that the buffer distance is set as a positive integer")             

            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:         
            if network_file is None:
                if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                    raise TypeError("Please make sure that the buffer distance is set as a positive integer")

                if network_type not in ["walk", "bike", "drive", "all"]:
                    raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                
                if not isinstance(travel_speed, int) or (not travel_speed > 0):
                    raise TypeError("Please make sure that the travel speed is set as a positive integer")

                if not isinstance(trip_time, int) or (not trip_time > 0):
                    raise TypeError("Please make sure that the trip time is set as a positive integer") 

                if geom_type == "Polygon":
                    print("Changing geometry type to Point by computing polygon centroids so that network can be retrieved...")
                    poi['geometry'] = poi['geometry'].centroid
                    print("Done \n")            

                meters_per_minute = travel_speed * 1000 / 60  # km per hour to m per minute
                epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326") # EPSG transformer to use OSM

                aoi_geometry = []
                for geom in tqdm(poi['geometry'], desc='Retrieving network for point(s) of interest'):
                    latlon = epsg_transformer.transform(geom.x, geom.y) # Transform point geometry into latlon for OSMnx
                    graph = ox.graph_from_point(latlon, network_type=network_type, dist=buffer_dist) # Retrieve street network for desired network type and buffer distance surrounding poi
                    graph = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") # Project street network graph back to original poi CRS
                    center_node = ox.distance.nearest_nodes(graph, geom.x, geom.y) # Find node which is closest to point location as base for next steps
                    for _, _, _, data in graph.edges(data=True, keys=True): # Calculate the time it takes to cover each edge's distance
                        data["time"] = data["length"] / meters_per_minute
                    isochrone_poly = make_iso_poly(graph, center_node=center_node, trip_time=trip_time) # See separate function for line by line explanation
                    aoi_geometry.append(isochrone_poly)

                aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
                print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")   
            else:
                if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
                    raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
                elif os.path.splitext(network_file)[1] == ".gpkg":
                    network = gpd.read_file(network_file, layer='edges')
                else:
                    network = gpd.read_file(network_file)
                
                if not network.crs.to_epsg() == epsg:
                    print("Adjusting CRS of Network file to match with Point of Interest CRS...")
                    network.to_crs(f'EPSG:{epsg}', inplace=True)
                    print("Done \n")

                # Create bounding box for network file
                bbox_network = network.unary_union.convex_hull

                if not all(geom.within(bbox_network) for geom in poi['geometry']):
                    raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")

                aoi_gdf = gpd.GeoDataFrame(geometry=[bbox_network], crs=f'EPSG:{epsg}')

    ### Step 3: Perform calculations and write results to file
    print("Calculating land cover class percentages...")
    if not all(geom.within(sg.box(*landcover_src.rio.bounds())) for geom in aoi_gdf['geometry']):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the Land Cover file, note that results will be based on the intersecting part of the buffer zone")
       
    # apply the function to each geometry in the GeoDataFrame and create a new Pandas Series
    landcover_percentages_series = aoi_gdf.geometry.apply(lambda x: pd.Series(calculate_landcover_percentages(landcover_src=landcover_src, geometry=x)))
    # rename the columns with the land cover class values
    landcover_percentages_series.columns = ["class_" + str(col) for col in landcover_percentages_series.columns]
    # concatenate the new series to the original dataframe
    poi = pd.concat([poi, landcover_percentages_series], axis=1)
    print("Done \n")

    print("Writing results to new geopackage file in specified directory...")
    input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
    poi.to_file(os.path.join(output_dir, f"{input_filename}_LCperc_added.gpkg"), driver="GPKG")
    print("Done")

    return poi

def get_canopy_percentage(point_of_interest_file, canopy_vector_file, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None,
                          network_type=None, trip_time=None, travel_speed=None, output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
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

    # Make sure poi dataframe contains ID column
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    canopy_src = gpd.read_file(canopy_vector_file)
    if not (canopy_src['geometry'].geom_type.isin(['Polygon', 'MultiPolygon']).all()):
        raise ValueError("Please make sure all geometries of the tree canopy file are of 'Polygon' or 'MultiPolygon' type and re-run the function")

    if not canopy_src.crs.to_epsg() == epsg:
        print("Adjusting CRS of Greenspace file to match with Point of Interest CRS...")
        canopy_src.to_crs(f'EPSG:{epsg}', inplace=True)
        print("Done \n")

    # Make sure all points of interest are within or do at least intersect (in case of polygons) the NDVI raster provided
    if not all(geom.within(sg.box(*canopy_src.total_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the tree canopy file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*canopy_src.total_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the tree canopy file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the tree canopy file provided, results will be based on intersecting part of polygons involved \n")

    ### Step 2: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        if geom_type == "Polygon":
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
        else:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")
    else:
        if buffer_type not in ["euclidian", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")

        if buffer_type == "euclidian":
            if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                raise TypeError("Please make sure that the buffer distance is set as a positive integer")             

            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:            
            if network_file is None:
                if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                    raise TypeError("Please make sure that the buffer distance is set as a positive integer")

                if network_type not in ["walk", "bike", "drive", "all"]:
                    raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                
                if not isinstance(travel_speed, int) or (not travel_speed > 0):
                    raise TypeError("Please make sure that the travel speed is set as a positive integer")

                if not isinstance(trip_time, int) or (not trip_time > 0):
                    raise TypeError("Please make sure that the trip time is set as a positive integer")   

                if geom_type == "Polygon":
                    print("Changing geometry type to Point by computing polygon centroids so that network can be retrieved...")
                    poi['geometry'] = poi['geometry'].centroid
                    print("Done \n")          
        
                meters_per_minute = travel_speed * 1000 / 60  # km per hour to m per minute
                epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326")

                aoi_geometry = []
                for geom in tqdm(poi['geometry'], desc='Retrieving network for point(s) of interest'):
                    latlon = epsg_transformer.transform(geom.x, geom.y) # Transform point geometry into latlon for OSMnx
                    graph = ox.graph_from_point(latlon, network_type=network_type, dist=buffer_dist) # Retrieve street network for desired network type and buffer distance surrounding poi
                    graph = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") # Project street network graph back to original poi CRS
                    center_node = ox.distance.nearest_nodes(graph, geom.x, geom.y) # Find node which is closest to point location as base for next steps
                    for _, _, _, data in graph.edges(data=True, keys=True): # Calculate the time it takes to cover each edge's distance
                        data["time"] = data["length"] / meters_per_minute
                    isochrone_poly = make_iso_poly(graph, center_node=center_node, trip_time=trip_time) # See separate function for line by line explanation
                    aoi_geometry.append(isochrone_poly)

                aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
                print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")    
            else:
                if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
                    raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
                elif os.path.splitext(network_file)[1] == ".gpkg":
                    network = gpd.read_file(network_file, layer='edges')
                else:
                    network = gpd.read_file(network_file)

                if not network.crs.to_epsg() == epsg:
                    print("Adjusting CRS of Network file to match with Point of Interest CRS...")
                    network.to_crs(f'EPSG:{epsg}', inplace=True)
                    print("Done \n")

                # Create bounding box for network file
                bbox_network = network.unary_union.convex_hull

                if not all(geom.within(bbox_network) for geom in poi['geometry']):
                    raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")

                aoi_gdf = gpd.GeoDataFrame(geometry=[bbox_network], crs=f'EPSG:{epsg}')

    ### Step 3: Perform calculations and write results to file
    print("Calculating percentage of tree canopy coverage...")
    if not all(geom.within(sg.box(*canopy_src.total_bounds)) for geom in aoi_gdf['geometry']):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the tree canopy file, note that results will be based on the intersecting part of the buffer zone")

    # Calculate percentage of tree canopy cover   
    poi['canopy_cover'] = aoi_gdf.apply(lambda row: str(((canopy_src.clip(row.geometry).area.sum()/row.geometry.area)*100).round(2))+'%', axis=1)
    print("Done \n")

    print("Writing results to new geopackage file in specified directory...")
    input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
    poi.to_file(os.path.join(output_dir, f"{input_filename}_CanopyPerc_added.gpkg"), driver="GPKG")
    print("Done")

    return poi

def get_park_percentage(point_of_interest_file, park_vector_file=None, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None,
                        network_type=None, trip_time=None, travel_speed=None, output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
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

    # Make sure poi dataframe contains ID column
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    ### Step 2: Read park data, retrieve from OSM if not provided by user
    epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326")
    if park_vector_file is None:
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise TypeError("Please make sure that the buffer distance is set as a positive integer") 

        print(f"Retrieving parks within buffer distance for {geom_type}(s) of interest...")
        if geom_type == "Polygon":
            park_retrieval = gpd.GeoDataFrame(geometry=poi['geometry'].centroid)
        else:
            park_retrieval = gpd.GeoDataFrame(geometry=poi['geometry'])

        # Tags seen as Urban Greenspace (UGS) require the following:
        # 1. Tag represent an area
        # 2. The area is outdoor
        # 3. The area is (semi-)publically available
        # 4. The area is likely to contain trees, grass and/or greenery
        # 5. The area can reasonable be used for walking or recreational activities
        park_tags = {'landuse':['allotments','forest','greenfield','village_green'], 'leisure':['garden','fitness_station','nature_reserve','park','playground'],'natural':'grassland'}
        park_src = gpd.GeoDataFrame()
        for geom in park_retrieval['geometry']:
            latlon = epsg_transformer.transform(geom.x, geom.y)
            park_geom = ox.geometries_from_point(latlon, tags=park_tags, dist=buffer_dist)
            park_src = gpd.GeoDataFrame(pd.concat([park_src, park_geom], ignore_index=True), crs=park_geom.crs)
        park_src.to_crs(f"EPSG:{epsg}", inplace=True)
        print("Done \n")
    else:
        park_src = gpd.read_file(park_vector_file)
        if not (park_src['geometry'].geom_type.isin(['Polygon', 'MultiPolygon']).all()):
            raise ValueError("Please make sure all geometries of the park file are of 'Polygon' or 'MultiPolygon' type and re-run the function")
        
        if not park_src.crs.to_epsg() == epsg:
            print("Adjusting CRS of Greenspace file to match with Point of Interest CRS...")
            park_src.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

    # Make sure all points of interest are within or do at least intersect (in case of polygons) the NDVI raster provided
    if not all(geom.within(sg.box(*park_src.total_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the park file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*park_src.total_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the park file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the park file provided, results will be based on intersecting part of polygons involved \n")

    ### Step 3: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        if geom_type == "Polygon":
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
        else:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")
    else:
        if buffer_type not in ["euclidian", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidian' or 'network' and re-run the function")

        if buffer_type == "euclidian":
            if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                raise TypeError("Please make sure that the buffer distance is set as a positive integer")             

            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:
            if network_file is None:
                if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                    raise TypeError("Please make sure that the buffer distance is set as a positive integer")

                if network_type not in ["walk", "bike", "drive", "all"]:
                    raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                
                if not isinstance(travel_speed, int) or (not travel_speed > 0):
                    raise TypeError("Please make sure that the travel speed is set as a positive integer")

                if not isinstance(trip_time, int) or (not trip_time > 0):
                    raise TypeError("Please make sure that the trip time is set as a positive integer") 

                if geom_type == "Polygon":
                    print("Changing geometry type to Point by computing polygon centroids so that network can be retrieved...")
                    poi['geometry'] = poi['geometry'].centroid
                    print("Done \n")            
        
                meters_per_minute = travel_speed * 1000 / 60  # km per hour to m per minute

                aoi_geometry = []
                for geom in tqdm(poi['geometry'], desc='Retrieving network for point(s) of interest'):
                    latlon = epsg_transformer.transform(geom.x, geom.y) # Transform point geometry into latlon for OSMnx
                    graph = ox.graph_from_point(latlon, network_type=network_type, dist=buffer_dist) # Retrieve street network for desired network type and buffer distance surrounding poi
                    graph = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") # Project street network graph back to original poi CRS
                    center_node = ox.distance.nearest_nodes(graph, geom.x, geom.y) # Find node which is closest to point location as base for next steps
                    for _, _, _, data in graph.edges(data=True, keys=True): # Calculate the time it takes to cover each edge's distance
                        data["time"] = data["length"] / meters_per_minute
                    isochrone_poly = make_iso_poly(graph, center_node=center_node, trip_time=trip_time) # See separate function for line by line explanation
                    aoi_geometry.append(isochrone_poly)

                aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
                print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")  
            else:
                if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
                    raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
                elif os.path.splitext(network_file)[1] == ".gpkg":
                    network = gpd.read_file(network_file, layer='edges')
                else:
                    network = gpd.read_file(network_file)

                if not network.crs.to_epsg() == epsg:
                    print("Adjusting CRS of Network file to match with Point of Interest CRS...")
                    network.to_crs(f'EPSG:{epsg}', inplace=True)
                    print("Done \n")

                # Create bounding box for network file
                bbox_network = network.unary_union.convex_hull

                if not all(geom.within(bbox_network) for geom in poi['geometry']):
                    raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")

                aoi_gdf = gpd.GeoDataFrame(geometry=[bbox_network], crs=f'EPSG:{epsg}')

    ### Step 4: Perform calculations and write results to file
    print("Calculating percentage of park area coverage...")
    if not all(geom.within(sg.box(*park_src.total_bounds)) for geom in aoi_gdf['geometry']):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the park file, note that results will be based on the intersecting part of the buffer zone")

    # Calculate percentage of park area cover   
    poi['park_cover'] = aoi_gdf.apply(lambda row: str(((park_src.clip(row.geometry).area.sum()/row.geometry.area)*100).round(2))+'%', axis=1)
    print("Done \n")

    print("Writing results to new geopackage file in specified directory...")
    input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
    poi.to_file(os.path.join(output_dir, f"{input_filename}_ParkPerc_added.gpkg"), driver="GPKG")
    print("Done")

    return poi

##### SUPPORTING FUNCTIONS
# Function to create isochrone polygon of network
def make_iso_poly(G, edge_buff=25, node_buff=0, center_node=None, trip_time=None):
    #Note: based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb
    subgraph = nx.ego_graph(G, center_node, radius=trip_time, distance="time") # Create sub graph of the street network which contains only parts which can be reached within specified travel parameters

    node_points = [sg.Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)] # Create list of point geometries existing of x and y coordinates for each node in subgraph retrieved from previous step
    nodes_gdf = gpd.GeoDataFrame({"id": list(subgraph.nodes)}, geometry=node_points) # Create geodataframe containing data from previous step
    nodes_gdf = nodes_gdf.set_index("id") # Set index to node ID

    edge_lines = []
    for n_fr, n_to in subgraph.edges(): # Iterate over edges in subgraph
        f = nodes_gdf.loc[n_fr].geometry # Retrieve geometry of the 'from' node of the edge
        t = nodes_gdf.loc[n_to].geometry # Retrieve geometry of the 'to' node of the edge
        edge_lookup = G.get_edge_data(n_fr, n_to)[0].get("geometry", sg.LineString([f, t])) # Retrieve edge geometry between from and to nodes
        edge_lines.append(edge_lookup) # Append edge geometry to list of edge lines

    n = nodes_gdf.buffer(node_buff).geometry # Create buffer around the nodes
    e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry # Create buffer around the edges
    all_gs = list(n) + list(e) # Concatenate nodes and edges
    isochrone_poly = gpd.GeoSeries(all_gs).unary_union # Create polygon of the concatenated nodes and edges

    isochrone_poly = sg.Polygon(isochrone_poly.exterior) # try to fill in surrounded areas so shapes will appear solid and blocks without white space inside them
    
    return isochrone_poly

# Function to calculate land cover percentages for a single geometry
def calculate_landcover_percentages(landcover_src, geometry):
    clipped = landcover_src.rio.clip([geometry]).clip(min=0) # Clip landcover raster to area of interest
    unique, counts = np.unique(clipped.values, return_counts=True) # Count the occurrences of all unique raster values
    total = counts.sum() # Calculate total nr. of occurrences 
    percentages = {value: str((count / total * 100).round(3)) + "%" for value, count in zip(unique, counts)} # Calculate percentages for each class
    return percentages