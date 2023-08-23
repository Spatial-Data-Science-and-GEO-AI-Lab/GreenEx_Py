# Data manipulation and analysis
import numpy as np
import pandas as pd

# File and directory operations
import os

# Geospatial data processing and analysis
import geopandas as gpd
import osmnx as ox
import networkx as nx
import rioxarray
import xrspatial
from rasterio.enums import Resampling
import pyproj
import shapely.geometry as sg
from shapely.ops import transform
import folium
import branca.colormap as cmp

# Geospatial data access and catalogs
import pystac_client
import planetary_computer
import odc.stac
import ee

# Date and time manipulation
from datetime import datetime, timedelta
from time import time

# Progress tracking
from tqdm.auto import tqdm

# Images
from PIL import Image
import requests
from io import BytesIO
from IPython.display import display

##### MAIN FUNCTIONS
def get_mean_NDVI(point_of_interest, ndvi_raster_file=None, crs_epsg=None, polygon_type="neighbourhood", buffer_type=None, 
                  buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, year=datetime.now().year, plot_aoi=True,
                  write_to_file=True, save_ndvi=True, output_dir=os.getcwd()):
    ### Step 1: Read and process user inputs, check conditions
    if isinstance(point_of_interest, gpd.GeoDataFrame):
        poi = point_of_interest
    elif os.path.isfile(point_of_interest):
        poi = gpd.read_file(point_of_interest)
    else:
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise ValueError("Please make sure that the buffer_dist argument is set to a positive integer")
        print("Retrieving street network through OpenStreetMap...")
        start_streetnetwork_retrieval = time()
        # Retrieve graph in WGS84
        poi_graph_wgs = ox.graph_from_place(point_of_interest).to_undirected()
        # Extract edges and store in geodataframe
        poi = ox.graph_to_gdfs(poi_graph_wgs, nodes=False, edges=True).reset_index()
        poi = poi.drop(columns=['u', 'v', 'key', 'osmid'])
        end_streetnetwork_retrieval = time()
        elapsed_streetnetwork_retrieval = end_streetnetwork_retrieval - start_streetnetwork_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_streetnetwork_retrieval))} \n")

    # Make sure CRS is projected rather than geographic
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

    # Create epsg transformer to use planetary computer and OSM
    epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326")

    if isinstance(point_of_interest, gpd.GeoDataFrame) or os.path.isfile(point_of_interest):
        # Verify that locations are either all provided using point geometries or all provided using polygon geometries
        if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon') or all(poi['geometry'].geom_type == 'MultiPolygon'):
            geom_type = poi.iloc[0]['geometry'].geom_type
        else:
            raise TypeError("Please make sure all geometries are of 'Point' type, all geometries are of 'Polygon' type or all geometries are of 'MultiPolygon' type and re-run the function")
        
        # In case of house polygons, transform to centroids
        if geom_type == "Polygon" or geom_type == "MultiPolygon":
            if polygon_type not in ["neighbourhood", "house"]:
                raise ValueError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
            if polygon_type == "house":
                print("Changing geometry type to Point by computing polygon centroids...")
                poi['geometry'] = poi['geometry'].centroid
                geom_type = poi.iloc[0]['geometry'].geom_type
                print("Done \n")

        # Make sure buffer distance and type are set in case of point geometries
        if geom_type == "Point":
            if buffer_type not in ["euclidean", "network"]:
                raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function") 

        # Make sure the buffer_type argument has a valid value if not None
        if buffer_type is not None and buffer_type not in ["euclidean", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

        # If buffer type is set to euclidean, make sure that the buffer distance is set
        if buffer_type == "euclidean":
                if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                    raise ValueError("Please make sure that the buffer_dist argument is set to a positive integer")  
        
        # If buffer type is set to network, make sure that either the buffer distance is set or both trip_time and travel_speed are set
        if buffer_type == "network":
            if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
                if not isinstance(travel_speed, (float, int)) or (not travel_speed > 0) or (not isinstance(trip_time, int) or (not trip_time > 0)):
                    raise TypeError("Please make sure that either the buffer_dist argument is set to a positive integer or both the travel_speed and trip_time are set to positive numbers")
                else:
                    # Convert km per hour to m per minute
                    meters_per_minute = travel_speed * 1000 / 60 
                    # Calculate max distance that can be travelled based on argument specified by user
                    buffer_dist = trip_time * meters_per_minute 
            else:
                # Buffer_dist and combination of travel_speed and trip_time cannot be set at same time
                if isinstance(travel_speed, (float, int)) and travel_speed > 0 and isinstance(trip_time, int) and trip_time > 0:
                    raise ValueError("Please make sure that one of the following requirements is met:\
                                    \n1. If buffer_dist is set, travel_speed and trip_time should not be set\
                                    \n2. If travel_speed and trip_time are set, buffer_dist shoud not be set")

        # Create polygon in which all pois are located to extract data from PC/OSM, incl. buffer if specified
        if buffer_dist is None:
            poi_polygon = sg.box(*poi.total_bounds)
        else:
            poi_polygon = sg.box(*poi.total_bounds).buffer(buffer_dist*1.10) # Add 10% to account for edge effects
    else:
        # Set buffer of 10m around the OSM network edges
        aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        # Determine bounding box of area of interest
        poi_polygon = sg.box(*aoi_gdf.total_bounds)

    # Make sure poi dataframe contains ID column
    if "id" in poi.columns:
        if poi['id'].isnull().values.any():
            poi['id'] = poi['id'].fillna(pd.Series(range(1, len(poi) + 1))).astype(int)
    else:
        poi['id'] = pd.Series(range(1, len(poi) + 1)).astype(int)

    # Retrieve NDVI raster, use planetary computer if not provided by user 
    if ndvi_raster_file is None:
        print("Retrieving NDVI raster through planetary computer...")
        start_ndvi_retrieval = time()

        # Transform CRS to comply with planetary computer requirements
        bounding_box_pc = transform(epsg_transformer.transform, poi_polygon).bounds 
        # Swap coords order to match with planetary computer format
        bounding_box_pc = [bounding_box_pc[1], bounding_box_pc[0], bounding_box_pc[3], bounding_box_pc[2]] 

        # Query planetary computer
        catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
        # Obtain Area of Interest
        time_of_interest = f"{year}-01-01/{year}-12-31" 
        # Search Data
        search = catalog.search(collections=["sentinel-2-l2a"],
                                bbox=bounding_box_pc,
                                datetime=time_of_interest,
                                query={"eo:cloud_cover": {"lt": 20}})
        # Obtain Data
        items = search.item_collection()
        # Create dataframe from planetary computer's item collection dictionary
        items_df = gpd.GeoDataFrame.from_features(items.to_dict(), crs="epsg:4326")
        # Make sure only images are maintained that contain all points/polygons of interest
        items_df_poi = items_df[items_df.geometry.contains(sg.box(*bounding_box_pc))]
        # Determine lowest percentage of cloud cover among filtered items
        lowest_cloud_cover = items_df_poi['eo:cloud_cover'].min()
        # Filter the satellite image which has the lowest cloud cover percentage
        item_to_select = items_df_poi[items_df_poi['eo:cloud_cover'] == lowest_cloud_cover]
        # Select item that matches the filters above and will be used to compose ndvi raster
        selected_item = next(item for item in items if item.properties["s2:granule_id"] == item_to_select.iloc[0]['s2:granule_id'])
        # Obtain Bands of Interest
        selected_item_data = odc.stac.stac_load([selected_item], bands = ['red', 'green', 'blue', 'nir'], bbox = bounding_box_pc).isel(time=0)
        # Calculate NDVI values
        ndvi = xrspatial.multispectral.ndvi(selected_item_data['nir'], selected_item_data['red'])
        # Reproject to original poi CRS
        ndvi_src = ndvi.rio.reproject(f"EPSG:{epsg}", resampling= Resampling.nearest, nodata=np.nan)

        # Provide information on satellite image that was used to user
        print(f"Information on the satellite image retrieved from planetary computer, used to calculate NDVI values:\
              \n   Date on which image was generated: {selected_item.properties['s2:generation_time']}\
              \n   Percentage of cloud cover: {selected_item.properties['eo:cloud_cover']}\
              \n   Percentage of pixels with missing data: {selected_item.properties['s2:nodata_pixel_percentage']}")

        # Save satellite image that was used in case user specifies so
        if save_ndvi:
            # Retrieve the image URL
            image_url = selected_item.assets["rendered_preview"].href
            # Download the image data
            response = requests.get(image_url)
            # Create a PIL Image object from the downloaded image data
            image = Image.open(BytesIO(response.content))
            # Create directory if the one specified by the user does not yet exist
            os.makedirs(output_dir, exist_ok=True)
            # Ask for filename to use
            input_filename_satellite = input("Enter filename for satellite image: ")
            input_filename_ndvi = input("Enter filename for ndvi raster: ")
            # Save the image to a file
            image.save(os.path.join(output_dir, f"{input_filename_satellite}.tif"))
            ndvi_src.rio.to_raster(os.path.join(output_dir, f"{input_filename_ndvi}.tif"))
            print("Satellite image and created NDVI raster successfully saved to file")
        end_ndvi_retrieval = time()
        elapsed_ndvi_retrieval = end_ndvi_retrieval - start_ndvi_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_ndvi_retrieval))} \n")
    else:
        # Read ndvi raster provided by user
        ndvi_src = rioxarray.open_rasterio(ndvi_raster_file)
        # Make sure that ndvi raster has same CRS as poi file
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

    if isinstance(point_of_interest, gpd.GeoDataFrame) or os.path.isfile(point_of_interest):
        ### Step 2: Construct the Area of Interest based on the arguments as defined by user
        if buffer_type is None:
            # Buffer type == None implies that provided polygons serve as areas of interest
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
        else:
            if buffer_type == "euclidean":
                # Create area of interest based on euclidean distance
                aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
            else:
                # Make sure network type argument has valid value
                if network_type not in ["walk", "bike", "drive", "all"]:
                    raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                
                # If poi file still contains polygon geometries, compute centroids so that isochrone maps can be created
                if geom_type == "Polygon" or geom_type == "MultiPolygon":
                    print("Changing geometry type to Point by computing polygon centroids so that isochrones can be retrieved...")
                    poi['geometry'] = poi['geometry'].centroid
                    print("Done \n")
                
                print("Retrieving network within total bounds of point(s) of interest, extended by buffer distance as specified...")
                start_network_retrieval = time()
                # Transform total bounds polygon of poi file to 4326 for OSM
                polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") 
                # Extract polygon in EPSG 4326 
                wgs_polygon = polygon_gdf_wgs['geometry'].values[0]        
                # Retrieve street network for desired network type
                graph = ox.graph_from_polygon(wgs_polygon, network_type=network_type) 
                # Project street network graph back to original poi CRS
                graph_projected = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") 
                end_network_retrieval = time()
                elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
                print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")             

                # Compute isochrone areas for points of interest
                aoi_geometry = []
                for geom in tqdm(poi['geometry'], desc = 'Retrieving isochrone for point(s) of interest'):
                    # Find node which is closest to point location as base for next steps
                    center_node = ox.distance.nearest_nodes(graph_projected, geom.x, geom.y) 
                    # Create sub graph of the street network which contains only parts which can be reached within specified travel parameters
                    subgraph = nx.ego_graph(graph_projected, center_node, radius=buffer_dist, distance="length") 
                    # Compute isochrones, see separate function for line by line explanation
                    isochrone_poly = make_iso_poly(graph_projected=graph_projected, subgraph=subgraph) 
                    aoi_geometry.append(isochrone_poly)

                # Create geodataframe of isochrone geometries
                aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
                print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")             
    
    ### Step 3: Calculate mean NDVI values and write results to file
    print("Calculating mean NDVI values...")
    start_calc = time()

    # Check whether areas of interest were succesfully created for all PoIs
    if not all(geom is not None for geom in aoi_gdf['geometry']):
        print(f"Warning: Buffer zones could not be created for all {geom_type}s of Interest based on the current argument values, resulting in missing values for the mean NDVI")
    # Check whether areas of interest, created in previous steps, are fully covered by the ndvi raster, provide warning if not
    if not all(geom.within(sg.box(*ndvi_src.rio.bounds())) for geom in aoi_gdf['geometry'] if geom is not None):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the NDVI raster, note that results will be based on the intersecting part of the buffer zone")

    # Calculate mean ndvi for geometries in poi file
    poi[['mean_NDVI', 'std_NDVI']] = aoi_gdf.apply(lambda row: pd.Series([np.nan, np.nan]) if row.geometry is None else pd.Series([ndvi_src.rio.clip([row.geometry]).clip(min=0).mean().values.round(3), ndvi_src.rio.clip([row.geometry]).clip(min=0).std().values.round(3)]), axis=1)
    end_calc = time()
    elapsed_calc = end_calc - start_calc
    print(f"Done, running time: {str(timedelta(seconds=elapsed_calc))} \n")

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        # Create directory if output directory specified by user does not yet exist
        os.makedirs(output_dir, exist_ok=True)
        # Ask for filename to use
        input_filename = input("Enter filename for results: ")
        poi.to_file(os.path.join(output_dir, f"{input_filename}.gpkg"), driver="GPKG")
        print("Done")

    if plot_aoi:
        # Get the total bounds of the poi file
        poi_total_bounds_wgs = poi.to_crs("EPSG:4326").total_bounds
        # Calculate the center coordinates
        center_coords = [(poi_total_bounds_wgs[1] + poi_total_bounds_wgs[3]) / 2, (poi_total_bounds_wgs[0] + poi_total_bounds_wgs[2]) / 2]
        # Create GeoJSON layers from the GeoDataFrames
        poi_column_names = list(filter(lambda col: col != 'geometry', poi.columns))

        if isinstance(point_of_interest, gpd.GeoDataFrame) or os.path.isfile(point_of_interest):
            # Create a base map
            map = folium.Map(location=center_coords, zoom_start=10)
            folium.GeoJson(poi.to_crs("EPSG:4326"),
                        name="PoI",
                        tooltip=folium.features.GeoJsonTooltip(fields=poi_column_names),
                        style_function=lambda feature: {'color': 'black'}).add_to(map)
            folium.GeoJson(aoi_gdf.to_crs("EPSG:4326"),
                        name="Buffer zones",
                        style_function=lambda feature: {'fillColor': 'blue', 'color': 'blue', 'fillOpacity': 0.3}).add_to(map)
            map_title = 'Areas of interest used for mean NDVI calculation'
        else:
            # Create color map for NDVI values in map
            linear = cmp.LinearColormap(['white', 'lightgreen', 'forestgreen', 'darkgreen'], vmin=0, vmax=1, caption='Mean NDVI')
            # Create a base map
            map = folium.Map(location=center_coords, zoom_start=15, tiles="cartodb positron")
            folium.GeoJson(poi.to_crs("EPSG:4326"),
                        name="Streets",
                        tooltip=folium.features.GeoJsonTooltip(fields=poi_column_names),
                        style_function=lambda feature: {'color': linear(feature['properties']['mean_NDVI']), 'weight': 5, 'fillColor': linear(feature['properties']['mean_NDVI']), 'fillOpacity': 0.8}).add_to(map)
            map_title = f'Mean NDVI for street network of {point_of_interest}'
            # Add the color legend to the map
            linear.caption = "Mean NDVI"
            linear.add_to(map)
        # Add layer control to the map
        folium.LayerControl().add_to(map)
        # Set the title
        map.get_root().html.add_child(folium.Element(f'<h3 style="text-align:center">{map_title}</h3>'))
        # Display map
        display(map)
    
    return poi

def get_landcover_percentages(point_of_interest, landcover_raster_file=None, crs_epsg=None, polygon_type="neighbourhood",
                              buffer_type=None, buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, plot_aoi=True,
                              write_to_file=True, save_lulc=True, output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
    if isinstance(point_of_interest, gpd.GeoDataFrame):
        poi = point_of_interest
    else:
        poi = gpd.read_file(point_of_interest)

    # Make sure that geometries in poi file are either all provided using point geometries or all using polygon geometries
    if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon') or all(poi['geometry'].geom_type == 'MultiPolygon'):
        geom_type = poi.iloc[0]['geometry'].geom_type
    else:
        raise TypeError("Please make sure all geometries are of 'Point' type, all geometries are of 'Polygon' type or all geometries are of 'MultiPolygon' type and re-run the function")

    # In case of house polygons, transform to centroids
    if geom_type == "Polygon" or geom_type == "MultiPolygon":
        if polygon_type not in ["neighbourhood", "house"]:
            raise ValueError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
        if polygon_type == "house":
            print("Changing geometry type to Point by computing polygon centroids...")
            poi['geometry'] = poi['geometry'].centroid
            geom_type = poi.iloc[0]['geometry'].geom_type
            print("Done \n")

    # Make sure buffer distance and type are set in case of point geometries
    if geom_type == "Point":
        if buffer_type not in ["euclidean", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # Make sure CRS is projected rather than geographic
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

    # Make sure the buffer_type argument has a valid value if not None
    if buffer_type is not None and buffer_type not in ["euclidean", "network"]:
        raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # If buffer type is set to euclidean, make sure that the buffer distance is set
    if buffer_type == "euclidean":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise TypeError("Please make sure that the buffer_dist argument is set to a positive integer")  
    
    # If buffer type is set to network, make sure that either the buffer distance is set or both trip_time and travel_speed are set
    if buffer_type == "network":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            if not isinstance(travel_speed, (float, int)) or (not travel_speed > 0) or (not isinstance(trip_time, int) or (not trip_time > 0)):
                raise ValueError("Please make sure that either the buffer_dist argument is set to a positive integer or both the travel_speed and trip_time are set to positive numbers")
            else:
                # Convert km per hour to m per minute
                meters_per_minute = travel_speed * 1000 / 60 
                # Calculate max distance that can be travelled based on argument specified by user and add 25% to account for edge effects
                buffer_dist = trip_time * meters_per_minute
        else:
            # Buffer_dist and combination of travel_speed and trip_time cannot be set at same time
            if isinstance(travel_speed, (float, int)) and travel_speed > 0 and isinstance(trip_time, int) and trip_time > 0:
                raise ValueError("Please make sure that one of the following requirements is met:\
                                \n1. If buffer_dist is set, travel_speed and trip_time should not be set\
                                \n2. If travel_speed and trip_time are set, buffer_dist shoud not be set")
            
    # Create polygon in which all pois are located to extract data from PC/OSM, incl. buffer if specified
    if buffer_dist is None:
        poi_polygon = sg.box(*poi.total_bounds)
    else:
        poi_polygon = sg.box(*poi.total_bounds).buffer(buffer_dist*1.10) # Add 10% to account for edge effects

    if landcover_raster_file is None:
        # Create epsg transformer to use planetary computer
        epsg_transformer = pyproj.Transformer.from_crs(f"epsg:{epsg}", "epsg:4326")
        print("Retrieving landcover class raster through planetary computer...")
        start_landcover_retrieval = time()
        # transform CRS to comply with planetary computer requirements
        bounding_box_pc = transform(epsg_transformer.transform, poi_polygon).bounds  
        # Swap coords order to match with planetary computer format
        bounding_box_pc = [bounding_box_pc[1], bounding_box_pc[0], bounding_box_pc[3], bounding_box_pc[2]] 

        # Query planetary computer
        catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
        search = catalog.search(
            collections=["esa-worldcover"],
            bbox=bounding_box_pc,
        )
        # Retrieve the items and select the first, most recent one
        items = search.item_collection()
        selected_item = items[0]
        # Extract landcover classes and store in dictionary to use in later stage
        class_list = selected_item.assets["map"].extra_fields["classification:classes"]
        classmap = {
            c["value"]: c["description"]
            for c in class_list
        }

        # Load raster using rioxarray
        landcover = rioxarray.open_rasterio(selected_item.assets["map"].href)
        # Clip raster to bounds of geometries in poi file
        landcover_clip = landcover.rio.clip_box(*bounding_box_pc)
        # Reproject to original poi file CRS
        landcover_src = landcover_clip.rio.reproject(f"EPSG:{epsg}", resampling= Resampling.nearest)
        
        # Provide landcover image information to user
        print(f"Information on the land cover image retrieved from planetary computer:\
              \n   Image description: {selected_item.properties['description']}\
              \n   Image timeframe: {selected_item.properties['start_datetime']} - {selected_item.properties['end_datetime']}")

        if save_lulc:
            # Create directory if the one specified by user does not yet exist
            os.makedirs(output_dir, exist_ok=True)
            # Ask for filename to use
            input_filename = input("Enter filename for landcover raster: ")
            # Write landcover raster to file
            landcover_src.rio.to_raster(os.path.join(output_dir, f"{input_filename}.tif"))
            print("Landcover image successfully saved to raster file")
        end_landcover_retrieval = time()
        elapsed_landcover_retrieval = end_landcover_retrieval - start_landcover_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_landcover_retrieval))} \n")
    else:    
        landcover_src = rioxarray.open_rasterio(landcover_raster_file)
        # Make sure landcover raster has same CRS as poi file
        if not landcover_src.rio.crs.to_epsg() == epsg:
            print("Adjusting CRS of land cover file to match with Point of Interest CRS...")
            landcover_src.rio.write_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Make sure all points of interest are within or do at least intersect (in case of polygons) the NDVI raster provided
        if not all(geom.within(sg.box(*landcover_src.rio.bounds())) for geom in poi['geometry']):
            if geom_type == "Point":
                raise ValueError("Not all points of interest are within the landcover file provided, please make sure they are and re-run the function")
            else:
                if not all(geom.intersects(sg.box(*landcover_src.rio.bounds())) for geom in poi['geometry']):
                    raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the landcover file provided, please make sure they are/do and re-run the function")
                else:
                    print("Warning: Not all polygons of interest are completely within the area covered by the landcover file provided, results will be based on intersecting part of polygons involved \n")

    ### Step 2: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        # Buffer type == None implies that polygons in poi file serve as areas of interest
        aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
    else:
        # Make sure buffer_dist is set in case buffer_type set to euclidean
        if buffer_type == "euclidean":
            # Create area of interest based on euclidean distance
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:
            # Make sure network_type argument has valid value
            if network_type not in ["walk", "bike", "drive", "all"]:
                raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")

            # In case poi still contains polygon geometries, compute centroids so that isochrones can be created
            if geom_type == "Polygon" or geom_type == "MultiPolygon":
                print("Changing geometry type to Point by computing polygon centroids so that isochrones can be retrieved...")
                poi['geometry'] = poi['geometry'].centroid
                print("Done \n")  
            
            print("Retrieving network within total bounds of point(s) of interest, extended by buffer distance as specified...")
            start_network_retrieval = time()
            # Transform bounds polygon of poi file to 4326 for OSM
            polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326")
            # Extract polygon in EPSG 4326    
            wgs_polygon = polygon_gdf_wgs['geometry'].values[0]     
            # Retrieve street network for desired network type
            graph = ox.graph_from_polygon(wgs_polygon, network_type=network_type) 
            # Project street network graph back to original poi CRS
            graph_projected = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") 
            end_network_retrieval = time()
            elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
            print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")         

            # Compose area of interest based on isochrones
            aoi_geometry = []
            for geom in tqdm(poi['geometry'], desc='Retrieving isochrone for point(s) of interest'):
                # Find node which is closest to point location as base for next steps
                center_node = ox.distance.nearest_nodes(graph_projected, geom.x, geom.y) 
                # Create sub graph of the street network which contains only parts which can be reached within specified travel parameters
                subgraph = nx.ego_graph(graph_projected, center_node, radius=buffer_dist, distance="length") 
                # Compute isochrones, see separate function for line by line explanation
                isochrone_poly = make_iso_poly(graph_projected=graph_projected, subgraph=subgraph) 
                aoi_geometry.append(isochrone_poly)

            # Create dataframe of isochrone polygons
            aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
            print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")         

    ### Step 3: Perform calculations and write results to file
    print("Calculating landcover class percentages...")
    start_calc = time()
    # Check whether areas of interest were succesfully created for all PoIs
    if not all(geom is not None for geom in aoi_gdf['geometry']):
        print(f"Warning: Buffer zones could not be created for all {geom_type}s of Interest based on the current argument values, resulting in missing values for the landcover class percentages")
    # Check whether areas of interest, created in previous steps, are fully covered by the landcover raster, provide warning if not
    if not all(geom.within(sg.box(*landcover_src.rio.bounds())) for geom in aoi_gdf['geometry'] if geom is not None):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the landcover raster, note that results will be based on the intersecting part of the buffer zone")
       
    # apply the landcover percentage function to each geometry in the GeoDataFrame and create a new Pandas Series
    landcover_percentages_series = aoi_gdf.geometry.apply(lambda x: pd.Series(calculate_landcover_percentages(landcover_src=landcover_src, geometry=x)))
    # rename the columns with the landcover class values
    if landcover_raster_file is None:
        landcover_percentages_series = landcover_percentages_series.rename(columns=lambda x: str(classmap.get(x, x)))
    else:
        landcover_percentages_series.columns = ["class_" + str(col) for col in landcover_percentages_series.columns]
    # concatenate the new series to the original dataframe
    poi = pd.concat([poi, landcover_percentages_series], axis=1)
    end_calc = time()
    elapsed_calc = end_calc - start_calc
    print(f"Done, running time: {str(timedelta(seconds=elapsed_calc))} \n")

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        # Create output directory if the one specified by user does not yet exist
        os.makedirs(output_dir, exist_ok=True)
        # Ask for filename to use
        input_filename = input("Enter filename for results: ")
        poi.to_file(os.path.join(output_dir, f"{input_filename}.gpkg"), driver="GPKG")
        print("Done")

    if plot_aoi:
        # Get the total bounds of the poi file
        poi_total_bounds_wgs = poi.to_crs("EPSG:4326").total_bounds
        # Calculate the center coordinates
        center_coords = [(poi_total_bounds_wgs[1] + poi_total_bounds_wgs[3]) / 2, (poi_total_bounds_wgs[0] + poi_total_bounds_wgs[2]) / 2]
        # Create a base map
        map = folium.Map(location=center_coords, zoom_start=10)
        # Create GeoJSON layers from the GeoDataFrames
        poi_column_names = list(filter(lambda col: col != 'geometry', poi.columns))
        folium.GeoJson(poi.to_crs("EPSG:4326"),
                       name="PoI",
                       tooltip=folium.features.GeoJsonTooltip(fields=poi_column_names),
                       style_function=lambda feature: {'color': 'black'}).add_to(map)
        folium.GeoJson(aoi_gdf.to_crs("EPSG:4326"),
                       name="Buffer zones",
                       style_function=lambda feature: {'fillColor': 'blue', 'color': 'blue', 'fillOpacity': 0.3}).add_to(map)
        # Add layer control to the map
        folium.LayerControl().add_to(map)
        # Set the title
        map_title = 'Areas of interest used for landcover class percentage calculation'
        map.get_root().html.add_child(folium.Element(f'<h3 style="text-align:center">{map_title}</h3>'))
        # Display map
        display(map)

    return poi

def get_canopy_percentage(point_of_interest, canopy_vector_file, crs_epsg=None, polygon_type="neighbourhood", buffer_type=None, 
                          buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, plot_aoi=True,
                          write_to_file=True, output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
    if isinstance(point_of_interest, gpd.GeoDataFrame):
        poi = point_of_interest
    else:
        poi = gpd.read_file(point_of_interest)

    # Make sure geometries of poi file are either all provided using point geometries or all using polygon geometries
    if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon') or all(poi['geometry'].geom_type == 'MultiPolygon'):
        geom_type = poi.iloc[0]['geometry'].geom_type
    else:
        raise TypeError("Please make sure all geometries are of 'Point' type, all geometries are of 'Polygon' type or all geometries are of 'MultiPolygon' type and re-run the function")

    # In case of house polygons, transform to centroids
    if geom_type == "Polygon" or geom_type == "MultiPolygon":
        if polygon_type not in ["neighbourhood", "house"]:
            raise ValueError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
        if polygon_type == "house":
            print("Changing geometry type to Point by computing polygon centroids...")
            poi['geometry'] = poi['geometry'].centroid
            geom_type = poi.iloc[0]['geometry'].geom_type
            print("Done \n")

    # Make sure buffer distance and type are set in case of point geometries
    if geom_type == "Point":
        if buffer_type not in ["euclidean", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # Make sure CRS is projected rather than geographic
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

    # Retrieve tree canopy data
    canopy_src = gpd.read_file(canopy_vector_file)
    # Make sure geometries in canopy file are of polygon or multipolygon as areas need to be calculated
    if not (canopy_src['geometry'].geom_type.isin(['Polygon', 'MultiPolygon']).all()):
        raise TypeError("Please make sure all geometries of the tree canopy file are of 'Polygon' or 'MultiPolygon' type and re-run the function")

    # Make sure canopy file has same CRS as poi file
    if not canopy_src.crs.to_epsg() == epsg:
        print("Adjusting CRS of Greenspace file to match with Point of Interest CRS...")
        canopy_src.to_crs(f'EPSG:{epsg}', inplace=True)
        print("Done \n")

    # Make sure all points of interest are within or do at least intersect (in case of polygons) the tree canopy file provided
    if not all(geom.within(sg.box(*canopy_src.total_bounds)) for geom in poi['geometry']):
        if geom_type == "Point":
            raise ValueError("Not all points of interest are within the tree canopy file provided, please make sure they are and re-run the function")
        else:
            if not all(geom.intersects(sg.box(*canopy_src.total_bounds)) for geom in poi['geometry']):
                raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the tree canopy file provided, please make sure they are/do and re-run the function")
            else:
                print("Warning: Not all polygons of interest are completely within the area covered by the tree canopy file provided, results will be based on intersecting part of polygons involved \n")

    # Make sure the buffer_type argument has a valid value if not None
    if buffer_type is not None and buffer_type not in ["euclidean", "network"]:
        raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # If buffer type is set to euclidean, make sure that the buffer distance is set
    if buffer_type == "euclidean":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise ValueError("Please make sure that the buffer_dist argument is set to a positive integer")  
    
    # If buffer type is set to network, make sure that either the buffer distance is set or both trip_time and travel_speed are set
    if buffer_type == "network":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            if not isinstance(travel_speed, (float, int)) or (not travel_speed > 0) or (not isinstance(trip_time, int) or (not trip_time > 0)):
                raise ValueError("Please make sure that either the buffer_dist argument is set to a positive integer or both the travel_speed and trip_time are set to positive numbers")
            else:
                # Convert km per hour to m per minute
                meters_per_minute = travel_speed * 1000 / 60 
                # Calculate max distance that can be travelled based on argument specified by user and add 25% to account for edge effects
                buffer_dist = trip_time * meters_per_minute
        else:
            # Buffer_dist and combination of travel_speed and trip_time cannot be set at same time
            if isinstance(travel_speed, (float, int)) and travel_speed > 0 and isinstance(trip_time, int) and trip_time > 0:
                raise ValueError("Please make sure that one of the following requirements is met:\
                                \n1. If buffer_dist is set, travel_speed and trip_time should not be set\
                                \n2. If travel_speed and trip_time are set, buffer_dist shoud not be set")
    
    # Create polygon in which all pois are located to extract data from PC/OSM, incl. buffer if specified
    if buffer_dist is None:
        poi_polygon = sg.box(*poi.total_bounds)
    else:
        poi_polygon = sg.box(*poi.total_bounds).buffer(buffer_dist*1.10) # Add 10% to account for edge effects
    
    ### Step 2: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        # Buffer type == None implies that polygon geometries serve as areas of interest
        aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
    else:
        # Make sure buffer dist is set in case buffer type set to euclidean
        if buffer_type == "euclidean":
            # Create area of interest based on euclidean buffer
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:     
            # Make sure network_type has valid value
            if network_type not in ["walk", "bike", "drive", "all"]:
                raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
                   
            # In case poi still contain polygon geometries, compute centroids so that isochrones can be created
            if geom_type == "Polygon" or geom_type == "MultiPolygon":
                print("Changing geometry type to Point by computing polygon centroids so that isochrone can be retrieved...")
                poi['geometry'] = poi['geometry'].centroid
                print("Done \n") 
            
            print("Retrieving network within total bounds of point(s) of interest, extended by buffer distance as specified...")
            start_network_retrieval = time()
            # Transform bounds polygon of poi file to 4326 for OSM
            polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") 
            # Extract polygon in EPSG 4326       
            wgs_polygon = polygon_gdf_wgs['geometry'].values[0]  
            # Retrieve street network for desired network type
            graph = ox.graph_from_polygon(wgs_polygon, network_type=network_type) 
            # Project street network graph back to original poi CRS
            graph_projected = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") 
            end_network_retrieval = time()
            elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
            print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")

            aoi_geometry = []
            for geom in tqdm(poi['geometry'], desc='Retrieving isochrone for point(s) of interest'):
                # Find node which is closest to point location as base for next steps
                center_node = ox.distance.nearest_nodes(graph_projected, geom.x, geom.y) 
                # Create sub graph of the street network which contains only parts which can be reached within specified travel parameters
                subgraph = nx.ego_graph(graph_projected, center_node, radius=buffer_dist, distance="length") 
                # Compute isochrones, see separate function for line by line explanation
                isochrone_poly = make_iso_poly(graph_projected=graph_projected, subgraph=subgraph) 
                aoi_geometry.append(isochrone_poly)

            # Create dataframe of isochrone polygons
            aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
            print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")    
            

    ### Step 3: Perform calculations and write results to file
    print("Calculating percentage of tree canopy coverage...")
    start_calc = time()
    # Check whether areas of interest were succesfully created for all PoIs
    if not all(geom is not None for geom in aoi_gdf['geometry']):
        print(f"Warning: Buffer zones could not be created for all {geom_type}s of Interest based on the current argument values, resulting in missing values for the canopy cover")
    # Check whether areas of interest, created in previous steps, are fully covered by the tree canopy file, provide warning if not
    if not all(geom.within(sg.box(*canopy_src.total_bounds)) for geom in aoi_gdf['geometry'] if geom is not None):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the tree canopy file, note that results will be based on the intersecting part of the buffer zone")

    # Calculate percentage of tree canopy cover   
    poi['canopy_cover'] = aoi_gdf.apply(lambda row: np.nan if row.geometry is None else str(((canopy_src.clip(row.geometry).area.sum()/row.geometry.area)*100).round(2))+'%', axis=1)
    end_calc = time()
    elapsed_calc = end_calc - start_calc
    print(f"Done, running time: {str(timedelta(seconds=elapsed_calc))} \n")

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        # Create directory if the one specified by the user does not yet exist
        os.makedirs(output_dir, exist_ok=True)
        # Ask for filename to use
        input_filename = input("Enter filename for results: ")
        poi.to_file(os.path.join(output_dir, f"{input_filename}.gpkg"), driver="GPKG")
        print("Done")

    if plot_aoi:
        # Get the total bounds of the poi file
        poi_total_bounds_wgs = poi.to_crs("EPSG:4326").total_bounds
        # Calculate the center coordinates
        center_coords = [(poi_total_bounds_wgs[1] + poi_total_bounds_wgs[3]) / 2, (poi_total_bounds_wgs[0] + poi_total_bounds_wgs[2]) / 2]
        # Create a base map
        map = folium.Map(location=center_coords, zoom_start=10)
        # Create GeoJSON layers from the GeoDataFrames
        poi_column_names = list(filter(lambda col: col != 'geometry', poi.columns))
        folium.GeoJson(poi.to_crs("EPSG:4326"),
                       name="PoI",
                       tooltip=folium.features.GeoJsonTooltip(fields=poi_column_names),
                       style_function=lambda feature: {'color': 'black'}).add_to(map)
        folium.GeoJson(canopy_src.to_crs("EPSG:4326"),
                       name="Tree canopy",
                       style_function=lambda feature: {'fillColor': 'green', 'color': 'green', 'fillOpacity': 0.7}).add_to(map)
        folium.GeoJson(aoi_gdf.to_crs("EPSG:4326"),
                       name="Buffer zones",
                       style_function=lambda feature: {'fillColor': 'blue', 'color': 'blue', 'fillOpacity': 0.3}).add_to(map)
        # Add layer control to the map
        folium.LayerControl().add_to(map)
        # Set the title
        map_title = 'Areas of interest used for tree canopy cover percentage calculation'
        map.get_root().html.add_child(folium.Element(f'<h3 style="text-align:center">{map_title}</h3>'))
        # Display map
        display(map)

    return poi

def get_greenspace_percentage(point_of_interest, greenspace_vector_file=None, crs_epsg=None, polygon_type="neighbourhood", buffer_type=None, 
                        buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, plot_aoi=True, write_to_file=True, 
                        output_dir=os.getcwd()):
    ### Step 1: Read and process user input, check conditions
    if isinstance(point_of_interest, gpd.GeoDataFrame):
        poi = point_of_interest
    else:
        poi = gpd.read_file(point_of_interest)

    # Make sure geometries of poi file are either all provided using point geometries or all using polygon geometries
    if all(poi['geometry'].geom_type == 'Point') or all(poi['geometry'].geom_type == 'Polygon') or all(poi['geometry'].geom_type == 'MultiPolygon'):
        geom_type = poi.iloc[0]['geometry'].geom_type
    else:
        raise TypeError("Please make sure all geometries are of 'Point' type, all geometries are of 'Polygon' type or all geometries are of 'MultiPolygon' type and re-run the function")
    
    # In case of house polygons, transform to centroids
    if geom_type == "Polygon" or geom_type == "MultiPolygon":
        if polygon_type not in ["neighbourhood", "house"]:
            raise ValueError("Please make sure that the polygon_type argument is set to either 'neighbourhood' or 'house'")
        if polygon_type == "house":
            print("Changing geometry type to Point by computing polygon centroids...")
            poi['geometry'] = poi['geometry'].centroid
            geom_type = poi.iloc[0]['geometry'].geom_type
            print("Done \n")

    # Make sure buffer distance and type are set in case of point geometries
    if geom_type == "Point":
        if buffer_type not in ["euclidean", "network"]:
            raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # Make sure CRS is projected rather than geographic
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
    
    # Make sure the buffer_type argument has a valid value if not None
    if buffer_type is not None and buffer_type not in ["euclidean", "network"]:
        raise ValueError("Please make sure that the buffer_type argument is set to either 'euclidean' or 'network' and re-run the function")

    # If buffer type is set to euclidean, make sure that the buffer distance is set
    if buffer_type == "euclidean":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            raise ValueError("Please make sure that the buffer_dist argument is set to a positive integer")  
    
    # If buffer type is set to network, make sure that either the buffer distance is set or both trip_time and travel_speed are set
    if buffer_type == "network":
        if not isinstance(buffer_dist, int) or (not buffer_dist > 0):
            if not isinstance(travel_speed, (float, int)) or (not travel_speed > 0) or (not isinstance(trip_time, int) or (not trip_time > 0)):
                raise ValueError("Please make sure that either the buffer_dist argument is set to a positive integer or both the travel_speed and trip_time are set to positive numbers")
            else:
                # Convert km per hour to m per minute
                meters_per_minute = travel_speed * 1000 / 60 
                # Calculate max distance that can be travelled based on argument specified by user and add 25% to account for edge effects
                buffer_dist = trip_time * meters_per_minute
        else:
            # Buffer_dist and combination of travel_speed and trip_time cannot be set at same time
            if isinstance(travel_speed, (float, int)) and travel_speed > 0 and isinstance(trip_time, int) and trip_time > 0:
                raise ValueError("Please make sure that one of the following requirements is met:\
                                \n1. If buffer_dist is set, travel_speed and trip_time should not be set\
                                \n2. If travel_speed and trip_time are set, buffer_dist shoud not be set")
    
    # Create polygon in which all pois are located to extract data from PC/OSM, incl. buffer if specified
    if buffer_dist is None:
        poi_polygon = sg.box(*poi.total_bounds)
    else:
        poi_polygon = sg.box(*poi.total_bounds).buffer(buffer_dist*1.10) # Add 10% to account for edge effects
    # Transform to 4326 for OSM
    polygon_gdf_wgs = gpd.GeoDataFrame(geometry=[poi_polygon], crs=f"EPSG:{epsg}").to_crs("EPSG:4326") 
    # Extract polygon in EPSG 4326 
    wgs_polygon = polygon_gdf_wgs['geometry'].values[0] 

    ### Step 2: Read greenspace data, retrieve from OSM if not provided by user
    if greenspace_vector_file is None:
        print(f"Retrieving greenspaces within total bounds of {geom_type}(s) of interest, extended by buffer distance if specified...")
        start_greenspace_retrieval = time()
        # Tags seen as Urban Greenspace (UGS) require the following:
        # 1. Tag represent an area
        # 2. The area is outdoor
        # 3. The area is (semi-)publically available
        # 4. The area is likely to contain trees, grass and/or greenery
        # 5. The area can reasonable be used for walking or recreational activities
        greenspace_tags = {'landuse':['allotments','forest','greenfield','village_green'], 'leisure':['garden','fitness_station','nature_reserve','park','playground'],'natural':'grassland'}
        # Extract greenspaces from OpenStreetMap
        greenspace_src = ox.geometries_from_polygon(wgs_polygon, tags=greenspace_tags)
        # Change CRS to the same one as poi file
        greenspace_src.to_crs(f"EPSG:{epsg}", inplace=True)
        # Create a boolean mask to filter out polygons and multipolygons
        polygon_mask = greenspace_src['geometry'].apply(lambda geom: geom.geom_type in ['Polygon', 'MultiPolygon'])
        # Filter the GeoDataFrame to keep only polygons and multipolygons
        greenspace_src = greenspace_src.loc[polygon_mask]
        end_greenspace_retrieval = time()
        elapsed_greenspace_retrieval = end_greenspace_retrieval - start_greenspace_retrieval
        print(f"Done, running time: {str(timedelta(seconds=elapsed_greenspace_retrieval))} \n")
    else:
        greenspace_src = gpd.read_file(greenspace_vector_file)
        # Make sure geometries are all polygons or multipolygons as areas should be calculated
        if not (greenspace_src['geometry'].geom_type.isin(['Polygon', 'MultiPolygon']).all()):
            raise TypeError("Please make sure all geometries of the greenspace file are of 'Polygon' or 'MultiPolygon' type and re-run the function")
        
        # Make sure CRS of greenspace file is same as CRS of poi file
        if not greenspace_src.crs.to_epsg() == epsg:
            print("Adjusting CRS of greenspace file to match with Point of Interest CRS...")
            greenspace_src.to_crs(f'EPSG:{epsg}', inplace=True)
            print("Done \n")

        # Make sure all points of interest are within or do at least intersect (in case of polygons) the greenspace file provided
        if not all(geom.within(sg.box(*greenspace_src.total_bounds)) for geom in poi['geometry']):
            if geom_type == "Point":
                raise ValueError("Not all points of interest are within the greenspace file provided, please make sure they are and re-run the function")
            else:
                if not all(geom.intersects(sg.box(*greenspace_src.total_bounds)) for geom in poi['geometry']):
                    raise ValueError("Not all polygons of interest are within, or do at least partly intersect, with the area covered by the greenspace file provided, please make sure they are/do and re-run the function")
                else:
                    print("Warning: Not all polygons of interest are completely within the area covered by the greenspace file provided, results will be based on intersecting part of polygons involved \n")

    # merge the overlapping greenspace into one polygon
    s_ = gpd.GeoDataFrame(geometry=[greenspace_src.unary_union], crs=f"EPSG:{epsg}").explode(index_parts=False).reset_index(drop=True)
    s_ = gpd.sjoin(s_, greenspace_src, how='left')
    greenspace_src = s_.dissolve(s_.index, aggfunc='first').reset_index(drop=True)
    
    ### Step 3: Construct the Area of Interest based on the arguments as defined by user
    if buffer_type is None:
        # Buffer type == None implies that polygon geometries serve as areas of interest
        aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'])
    else:
        # Make sure buffer dist is set in case buffer type is euclidean
        if buffer_type == "euclidean":
            # Create area of interest based on euclidean buffer
            aoi_gdf = gpd.GeoDataFrame(geometry=poi['geometry'].buffer(buffer_dist))
        else:
            # Make sure network type has valid value
            if network_type not in ["walk", "bike", "drive", "all"]:
                raise ValueError("Please make sure that the network_type argument is set to either 'walk', 'bike, 'drive' or 'all', and re-run the function")
            
            # If poi still contains polygon geometries, compute centroids so that isochrones can be created
            if geom_type == "Polygon" or geom_type == "MultiPolygon":
                print("Changing geometry type to Point by computing polygon centroids so that isochrones can be retrieved...")
                poi['geometry'] = poi['geometry'].centroid
                print("Done \n") 
            
            print(f"Retrieving network within total bounds of {geom_type}(s) of interest, extended by buffer distance as specified...")
            start_network_retrieval = time()       
            # Retrieve street network for desired network type
            graph = ox.graph_from_polygon(wgs_polygon, network_type=network_type) 
            # Project street network graph back to original poi CRS
            graph_projected = ox.project_graph(graph, to_crs=f"EPSG:{epsg}") 
            end_network_retrieval = time()
            elapsed_network_retrieval = end_network_retrieval - start_network_retrieval
            print(f"Done, running time: {str(timedelta(seconds=elapsed_network_retrieval))} \n")
                
            aoi_geometry = []
            for geom in tqdm(poi['geometry'], desc='Retrieving isochrone for point(s) of interest'):
                # Find node which is closest to point location as base for next steps
                center_node = ox.distance.nearest_nodes(graph_projected, geom.x, geom.y) 
                # Create sub graph of the street network which contains only parts which can be reached within specified travel parameters
                subgraph = nx.ego_graph(graph_projected, center_node, radius=buffer_dist, distance="length") 
                # Compute isochrones, see separate function for line by line explanation
                isochrone_poly = make_iso_poly(graph_projected=graph_projected, subgraph=subgraph) 
                aoi_geometry.append(isochrone_poly)

            # Create dataframe with isochrone geometries
            aoi_gdf = gpd.GeoDataFrame(geometry=aoi_geometry, crs=f"EPSG:{epsg}")
            print("Note: creation of isochrones based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb \n")  

    ### Step 4: Perform calculations and write results to file
    print("Calculating percentage of greenspace area coverage...")
    start_calc = time()
    # Check whether areas of interest were succesfully created for all PoIs
    if not all(geom is not None for geom in aoi_gdf['geometry']):
        print(f"Warning: Buffer zones could not be created for all {geom_type}s of Interest based on the current argument values, resulting in missing values for the mean NDVI")
    # Check whether areas of interest, created in previous steps, are fully covered by the greenspace bounds, provide warning if not
    if not all(geom.within(sg.box(*greenspace_src.total_bounds)) for geom in aoi_gdf['geometry'] if geom is not None):
        print(f"Warning: Not all buffer zones for the {geom_type}s of Interest are completely within the area covered by the greenspace file, note that results will be based on the intersecting part of the buffer zone")

    # Calculate percentage of greenspace area cover   
    poi['greenspace_cover'] = aoi_gdf.apply(lambda row: np.nan if row.geometry is None else str(((greenspace_src.clip(row.geometry).area.sum()/row.geometry.area)*100).round(2))+'%', axis=1)
    end_calc = time()
    elapsed_calc = end_calc - start_calc
    print(f"Done, running time: {str(timedelta(seconds=elapsed_calc))} \n")

    if write_to_file:
        print("Writing results to new geopackage file in specified directory...")
        # Create output directory if the one specified by user does not yet exist
        os.makedirs(output_dir, exist_ok=True)
        # Ask for filename to use
        input_filename = input("Enter filename for results: ")
        poi.to_file(os.path.join(output_dir, f"{input_filename}.gpkg"), driver="GPKG")
        print("Done")

    if plot_aoi:
        # Get the total bounds of the poi file
        poi_total_bounds_wgs = poi.to_crs("EPSG:4326").total_bounds
        # Calculate the center coordinates
        center_coords = [(poi_total_bounds_wgs[1] + poi_total_bounds_wgs[3]) / 2, (poi_total_bounds_wgs[0] + poi_total_bounds_wgs[2]) / 2]
        # Create a base map
        map = folium.Map(location=center_coords, zoom_start=10)
        # Create GeoJSON layers from the GeoDataFrames
        poi_column_names = list(filter(lambda col: col != 'geometry', poi.columns))
        folium.GeoJson(poi.to_crs("EPSG:4326"),
                       name="PoI",
                       tooltip=folium.features.GeoJsonTooltip(fields=poi_column_names),
                       style_function=lambda feature: {'color': 'black'}).add_to(map)
        folium.GeoJson(greenspace_src.to_crs("EPSG:4326"),
                       name="Greenspaces",
                       style_function=lambda feature: {'fillColor': 'green', 'color': 'green', 'fillOpacity': 0.7}).add_to(map)
        folium.GeoJson(aoi_gdf.to_crs("EPSG:4326"),
                       name="Buffer zones",
                       style_function=lambda feature: {'fillColor': 'blue', 'color': 'blue', 'fillOpacity': 0.3}).add_to(map)
        # Add layer control to the map
        folium.LayerControl().add_to(map)
        # Set the title
        map_title = 'Areas of interest used for greenspace cover percentage calculation'
        map.get_root().html.add_child(folium.Element(f'<h3 style="text-align:center">{map_title}</h3>'))
        # Display map
        display(map)

    return poi

##### SUPPORTING FUNCTIONS
# Function to create isochrone polygon of network
def make_iso_poly(graph_projected, subgraph, edge_buff=25, node_buff=0):
    #Note: based on code by gboeing, source: https://github.com/gboeing/osmnx-examples/blob/main/notebooks/13-isolines-isochrones.ipynb
    node_points = [sg.Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)] # Create list of point geometries existing of x and y coordinates for each node in subgraph retrieved from previous step
    nodes_gdf = gpd.GeoDataFrame({"id": list(subgraph.nodes)}, geometry=node_points) # Create geodataframe containing data from previous step
    nodes_gdf = nodes_gdf.set_index("id") # Set index to node ID

    edge_lines = []
    for n_fr, n_to in subgraph.edges(): # Iterate over edges in subgraph
        f = nodes_gdf.loc[n_fr].geometry # Retrieve geometry of the 'from' node of the edge
        t = nodes_gdf.loc[n_to].geometry # Retrieve geometry of the 'to' node of the edge
        edge_lookup = graph_projected.get_edge_data(n_fr, n_to)[0].get("geometry", sg.LineString([f, t])) # Retrieve edge geometry between from and to nodes
        edge_lines.append(edge_lookup) # Append edge geometry to list of edge lines

    n = nodes_gdf.buffer(node_buff).geometry # Create buffer around the nodes
    e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry # Create buffer around the edges
    all_gs = list(n) + list(e) # Concatenate nodes and edges
    
    try:
        isochrone_poly = gpd.GeoSeries(all_gs).unary_union # Create polygon of the concatenated nodes and edges
        isochrone_poly = sg.Polygon(isochrone_poly.exterior) # try to fill in surrounded areas so shapes will appear solid and blocks without white space inside them
    except:
        isochrone_poly = None

    return isochrone_poly

# Function to calculate land cover percentages for a single geometry
def calculate_landcover_percentages(landcover_src, geometry):
    if geometry is None:
        percentages = {}
        return percentages
    
    # Clip landcover raster to area of interest
    clipped = landcover_src.rio.clip([geometry]).clip(min=0) 
    # Count the occurrences of all unique raster values
    unique, counts = np.unique(clipped.values, return_counts=True) 
    # Calculate total nr. of occurrences 
    total = counts.sum() 
    # Calculate percentages for each class
    percentages = {value: str((count / total * 100).round(3)) + "%" for value, count in zip(unique, counts)} 
    
    return percentages