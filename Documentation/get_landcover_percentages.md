**get_landcover_percentages(point_of_interest_file, landcover_raster_file, buffer_type=None, buffer_dist=None,network_file=None, network_type=None, trip_time=None, travel_speed=None, output_dir=os.getcwd())**

> Retrieve the percentage of area covered by each land cover class for areas or points of interest.

> Function to calculate the percentage of area covered by each land cover class surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries around and for which to compute mean NDVI values.

>> - landcover_raster_file *(string)* – the absolute or relative path to the raster file containing a land cover classification map, where each pixel is assigned a land cover class.

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries and optional in case point_of_interest_file contains polygon geometries, the way in which the buffer distance should be considered.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidian" or "network", surrounding distance in meters to consider for land cover type percentage calculation.

>> - network_file *(string)* – may optionally be defined in case buffer_type is set to "network", the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified while both conditions are met, network will be retrieved through the OSMnx package.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which network buffer needs to be composed.

>> - trip_time *(int)* – to be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type.

>> - travel_speed *(int)* – to be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including columns for land cover class percentages. Dataframe will also be written to new file in specified directory (see output_dir argument). 

>>Return type:	
>>> Geodataframe
