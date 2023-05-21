**get_mean_NDVI(point_of_interest_file, ndvi_raster_file, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None, network_type=None, trip_time=None, travel_speed=None)**

> Retrieve the mean Normalised Difference Vegetation Index (NDVI) for areas or points of interest.

> Function to calculate mean NDVI values surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries around and for which to compute mean NDVI values.

>> - ndvi_raster_file *(string)* – optional, the absolute or relative path to the file containing the NDVI values in raster format. If not provided, the NDVI raster will be obtained by using satellite images from the planetary computer.  

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest file has geographic CRS rather than projected. CRS will be transformed to the one specified. In case crs_epsg is not specified and CRS of file is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries and optional in case point_of_interest_file contains polygon geometries, the way in which the buffer distance should be considered.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidian" or "network", surrounding distance in meters to consider for mean NDVI calculation.

>> - network_file *(string)* – may optionally be defined in case buffer_type is set to "network", the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified while both conditions are met, network will be retrieved through the OSMnx package.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which network buffer needs to be composed.

>> - trip_time *(int)* – to be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type.

>> - travel_speed *(int)* – to be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written. If not specified, the current working directory will serve as default.

>> - year *(int)* – optional, may be defined if no raster file is provided. The year for which to retrieve satellite images using the planetary computer. 

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column for mean NDVI value(s). Dataframe will also be written to new file in specified directory (see output_dir argument). 

>>Return type:	
>>> Geodataframe