**get_mean_NDVI(point_of_interest_file, ndvi_raster_file, buffer_dist=None, buffer_type=None, network_file=None, network_type=None, trip_time=None, travel_speed=None)**

> Retrieve the mean Normalised Difference Vegetation Index (NDVI) for an area or points of interest.

> Function to calculate mean NDVI values surrounding a geographic location or within a geographic area. The point surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) distance buffer. 

>> Parameters: 

>> - point_of_interest_file *(string)* – The absolute or relative path to the file containing point or polygon geometries around and for which to compute mean NDVI values.

>> - ndvi_raster_file *(string)* – The absolute or relative path to the file containing the NDVI values in raster format. 

>> - buffer_dist *(int)* – to be defined in case point_of_interest_file contains point geometries rather than polygon geometries, surrounding distance in meters to consider for mean NDVI calculation.

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries rather than polygon geometries, the way in which the distance buffer should be considered.

>> - network_file *(string)* – may optionally be defined in case point_of_interest_file contains point geometries rather than polygon geometries and buffer_type is set to "network", the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified while both conditions are met, network will be retrieved through the OSMnx package.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case point_of_interest_file contains point geometries rather than polygon geometries and buffer_type is set to "network", the travel mode for which network buffer needs to be composed.

>> - trip_time *(int)* – to be defined in case point_of_interest_file contains point geometries rather than polygon geometries and buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type.

>> - travel_speed *(int)* – to be defined in case point_of_interest_file contains point geometries rather than polygon geometries and buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column for mean NDVI value(s).

>>Return type:	
>>> Geodataframe