# Data manipulation and analysis
import numpy as np
import pandas as pd

# File and directory operations
import os

# Geospatial data processing and analysis
import geopandas as gpd
import osmnx as ox
import rioxarray
import pyproj
import shapely.geometry as sg
import networkx as nx
import momepy

##### MAIN FUNCTIONS
def get_shortest_distance_park(point_of_interest_file, crs_epsg=None, target_dist=300, park_vector_file=None, destination="centroids", 
                               network_file=None, network_type=None, output_dir=os.getcwd()):
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
        print("Changing geometry type to Point by computing polygon centroids so that network distance can be retrieved...")
        poi['geometry'] = poi['geometry'].centroid
        print("Done \n")

    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    if not isinstance(target_dist, int) or (not target_dist > 0):
        raise TypeError("Please make sure that the target distance is set as a positive integer")
    
    if destination not in ["centroids", "entrance"]:
        raise TypeError("Please make sure that the destination argument is set to either 'centroids' or 'entrance'")

    ### Step 2: Obtain bounding box in which all points of interest are located, including 1000m buffer to account for edge effects
    poi_polygon = sg.box(*poi.total_bounds).buffer(target_dist*1.5)
    polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") # Transform to 4326 for OSM
    wgs_polygon = polygon_gdf_wgs['geometry'].values[0] # Extract polygon in EPSG 4326

    ### Step 3: Read park polygons, retrieve from OSM if not provided by user 
    if park_vector_file is not None:
        park_src = gpd.read_file(park_vector_file)
        if not park_src.crs.to_epsg() == epsg:
            print("Adjusting CRS of Park file to match with Point of Interest CRS...")
            park_src.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")
    else:
        print(f"Retrieving parks within total bounds of point(s) of interest, extended by a {target_dist*1.5}m buffer to account for edge effects...")
        # Tags seen as Urban Greenspace (UGS) require the following:
        # 1. Tag represent an area
        # 2. The area is outdoor
        # 3. The area is (semi-)publically available
        # 4. The area is likely to contain trees, grass and/or greenery
        # 5. The area can reasonable be used for walking or recreational activities
        park_tags = {'landuse':['allotments','forest','greenfield','village_green'], 'leisure':['garden','fitness_station','nature_reserve','park','playground'],'natural':'grassland'}
        park_src = ox.geometries_from_polygon(wgs_polygon, tags=park_tags)
        park_src.to_crs(f"EPSG:{epsg}", inplace=True)
        print("Done \n")
    
    if destination == "centroids":
        park_src['centroid'] = park_src['geometry'].centroid
    
    park_src['park_id'] = list(range(len(park_src)))

    ### Step 3: Read network, retrieve from OSM if not provided by user 
    if network_file is not None:
        if os.path.splitext(network_file)[1] not in [".gpkg", ".shp"]:
            raise ValueError("Please provide the network file in '.gpkg' or '.shp' format")
        elif network_file is not None and (os.path.splitext(network_file)[1] == ".gpkg"):
            network = gpd.read_file(network_file, layer='edges')
        else: 
            network = gpd.read_file(network_file)

        if not network.crs.to_epsg() == epsg:
            print("Adjusting CRS of Network file to match with Point of Interest CRS...")
            network.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Check if house locations are within network file provided
        bbox_network = network.unary_union.envelope
        if not all(geom.within(bbox_network) for geom in poi['geometry']):
            raise ValueError("Not all points of interest are within the network file provided, please make sure they are and re-run the function")

        # Convert network to graph object using momempy
        network_graph = momepy.gdf_to_nx(network)
    else:
        if network_type not in ["walk", "bike", "drive", "all"]:
            raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
            
        print(f"Retrieving infrastructure network within total bounds of point(s) of interest, extended by a {target_dist*1.5}m buffer to account for edge effects...")
        network_graph = ox.graph_from_polygon(wgs_polygon, network_type=network_type)
        network_graph = ox.project_graph(network_graph, to_crs=f"EPSG:{epsg}")
        print("Done \n")
    
    ### Step 4: Perform calculations and write results to file
    print("Calculating shortest distances...")
    poi[[f'park_within_{target_dist}m', 'distance_to_park']] = poi.apply(lambda row: pd.Series(calculate_shortest_distance(df_row=row, target_dist=target_dist, network_graph=network_graph, park_src=park_src, destination=destination)), axis=1)
    print("Done \n")
    
    print("Writing results to new geopackage file in specified directory...")
    input_filename, _ = os.path.splitext(os.path.basename(point_of_interest_file))
    poi.to_file(os.path.join(output_dir, f"{input_filename}_ShortDistPark_added.gpkg"), driver="GPKG")
    print("Done")

    return poi

##### SUPPORTING FUNCTIONS
def calculate_shortest_distance(df_row=None, target_dist=None, network_graph=None, park_src=None, destination=None):   
    ### Step 1: Clip park boundaries to poi incl. buffer to minimize possible destination points
    park_src_buffer = park_src.clip(df_row['geometry'].buffer(target_dist))
    ### Step 2: Retrieve nearest network node for house location and calculate euclidian distance between these points
    # Euclidian distance will be added to network distance to minimize distance error
    nearest_node = ox.distance.nearest_nodes(network_graph, df_row['geometry'].x, df_row['geometry'].y)
    penalty_home = df_row['geometry'].distance(sg.Point(network_graph.nodes[nearest_node]['x'], network_graph.nodes[nearest_node]['y']))

    ### Step 3: Create subgraph to only consider network in point's proximity -- save time 
    subgraph = nx.ego_graph(network_graph, nearest_node, radius=target_dist*1.5, distance="length")

    ### Step 4: Create fake entry points for parks by getting network nodes that are within 20m of park boundaries
    pos = {n: (subgraph.nodes[n]['x'], subgraph.nodes[n]['y']) for n in subgraph.nodes} # Create dictionary to extract geometries for nodes of interest
    
    # For each park, retrieve the network nodes which are within 20m of the park boundary and store in dictionary
    park_boundary_nodes = {}
    for park_id, geom in zip(park_src_buffer['park_id'], park_src_buffer['geometry']):
        boundary_nodes = [node for node in subgraph.nodes() if sg.Point(pos[node]).distance(geom.boundary) < 20]
        park_boundary_nodes[park_id] = boundary_nodes
    
    ### Step 5: Calculate the network distances between the house location's nearest node and the fake park entry points
    # Add penalty_home as defined before to network distance, as well as penalty_centroid in case user defined destination argument as "centroids"
    distances = {}
    for park_id, boundary_nodes in park_boundary_nodes.items():
        for node in boundary_nodes:
            try:
                path = nx.shortest_path(subgraph, nearest_node, node, weight='length')
                if destination == "centroids": 
                    # Calculate euclidian distance between the fake park entry points and the corresponding park's centroid to minimize distance error
                    penalty_centroid = park_src_buffer[park_src_buffer['park_id'] == park_id]['centroid'].iloc[0].distance(sg.Point(subgraph.nodes[node]['x'], subgraph.nodes[node]['y']))
                    distance = sum([subgraph.edges[path[i], path[i+1], 0]['length'] for i in range(len(path)-1)]) + penalty_home + penalty_centroid
                else:
                    distance = sum([subgraph.edges[path[i], path[i+1], 0]['length'] for i in range(len(path)-1)]) + penalty_home
                distances[node] = distance
            except:
                continue

    # Get the minimum distance (house location to park)
    if distances:
        min_distance = round(min(distances.values()),0)
    else:
        min_distance = np.nan
    
    ### Step 6: Define result, if minimum distance smaller than/equal to target distance threshold --> Good
    if min_distance <= target_dist:
        outcome = "True"
    else:
        outcome = "False"
    
    return outcome, min_distance