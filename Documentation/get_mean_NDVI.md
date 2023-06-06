**get_mean_NDVI(point_of_interest_file, ndvi_raster_file=None, crs_epsg=None, polygon_type="neighbourhood", buffer_type=None, buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, year=datetime.now().year, write_to_file=True, save_ndvi=True, output_dir=os.getcwd())**

> Retrieve the mean Normalised Difference Vegetation Index (NDVI) for areas or points of interest.

> Function to calculate mean NDVI values surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries around and for which to compute mean NDVI values.

>> - ndvi_raster_file *(string)* – optional, the absolute or relative path to the file containing the NDVI values in raster format. If not provided, the NDVI raster will be obtained by using satellite images from the planetary computer.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest file has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of file is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - polygon_type *(string {"neighbourhood", "house"})* - to be defined in case point_of_interest_file contains polygon geometries.

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries and optional in case point_of_interest_file contains polygon geometries, the way in which the area of interest should be composed. If "euclidian", a straight line distance will be used based on the buffer distance as specified by the buffer_dist argument. If "network", isoschrone maps will be composed based on the buffer distance or combination of trip_time and travel_speed.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidian" and optional if buffer_type is set to "network". The distance in meters that will be used to compute the area of interest surrounding the geometries of the point_of_interest file. NOTE: In case buffer_type is set to "network", buffer_dist should not be specified if travel_speed and trip_time arguments are specified.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which the network needs to be retrieved from OpenStreetMap.

>> - trip_time *(int)* – may optionally be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type. The trip_time, in addition to the travel_speed, will be used to compose an isochrone map if no buffer distance is provided. NOTE: travel_speed and trip_time arguments should not be specified if buffer_dist is specified.

>> - travel_speed *(int)* – may optionally be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type. The travel_speed, in addition to the trip_time, will be used to compose an isochrone map if no buffer distance is provided. NOTE: travel_speed and trip_time arguments should not be specified if buffer_dist is specified.

>> - year *(int)* – optional, may be defined if no ndvi raster file is provided. The year for which to retrieve satellite images using the planetary computer. If not specified while raster file is not provided, the current year will be used by default.

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - save_ndvi *(bool {"TRUE", "FALSE"})* - whether or not to save the satellite image retrieved from planetary computer and NDVI raster that was created afterwards in case NDVI raster was not provided. Files will be written to a new file in the directory specified in the output_dir argument if set to TRUE. By default, argument is set to TRUE.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column for mean NDVI value(s). Dataframe will also be written to new file in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframe