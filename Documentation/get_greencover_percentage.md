**get_greencover_percentage(point_of_interest_file, greenspace_vector_file, buffer_type=None, buffer_dist=None,network_file=None, network_type=None, trip_time=None, travel_speed=None, output_dir=os.getcwd())**

> Retrieve the percentage of area covered by greenspace for areas or points of interest.

> Function to calculate the percentage of area covered by greenspace surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries around and for which to compute the greenspace cover percentage.

>> - greenspace_vector_file *(string)* – the absolute or relative path to the vector file containing polygons of greenspaces such as tree canopies and parks. 

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries and optional in case point_of_interest_file contains polygon geometries, the way in which the buffer distance should be considered.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidian" or "network", surrounding distance in meters to consider for land cover type percentage calculation.

>> - network_file *(string)* – may optionally be defined in case buffer_type is set to "network", the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified while both conditions are met, network will be retrieved through the OSMnx package.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which network buffer needs to be composed.

>> - trip_time *(int)* – to be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type.

>> - travel_speed *(int)* – to be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column for greenspace cover percentage. Dataframe will also be written to new file in specified directory (see output_dir argument). 

>>Return type:	
>>> Geodataframe
