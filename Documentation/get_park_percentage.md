**get_park_percentage(point_of_interest_file, park_vector_file=None, crs_epsg=None, buffer_type=None, buffer_dist=None, network_file=None, network_type=None, trip_time=None, travel_speed=None, write_to_file=True, output_dir=os.getcwd())**

> Retrieve the percentage of area covered by parks for areas or points of interest.

> Function to calculate the percentage of area covered by parks surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidian (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries around and for which to compute the percentage of park cover.

>> - park_vector_file *(string)* – optional, the absolute or relative path to the vector file containing tree canopy data. Note that geometries should be polygon or multipolygon. In case no file is provided, park data will be retrieved from OpenStreetMap.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest file has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of file is geographic, CRS will be transformed to EPSG 3395 by default.

>> - buffer_type *(string {"euclidian", "network"})* – to be defined in case point_of_interest_file contains point geometries and optional in case point_of_interest_file contains polygon geometries, the way in which the area of interest should be composed. If "euclidian", a straight line distance will be used based on the buffer distance as specified by the buffer_dist argument. If "network", isoschrone maps will be composed based on the additional arguments of trip_time and travel_speed.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidian" OR "network" while no network is provided. In case buffer_type is "euclidian", the buffer distance is used to define the area, surrouding the point(s)/polygon(s) of interest, for which the percentage of park cover should be calculated. In case buffer_type is "network", the buffer distance will be used for extracting data OpenStreetMap if the park vector file or network is not provided whereas the area of interest for which to perform the park coverage calculation will then be defined based on isochrones, following the trip_time and travel_mode arguments.

>> - network_file *(string)* – optional, may be defined in case buffer_type is set to "network", the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified while buffer_type is set to "network", network will be retrieved through OpenStreetMap considering the point/polygon's surrounding straight line distance as specified in the buffer_dist argument.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network" and no network file is provided, the travel mode for which the network needs to be retrieved.

>> - trip_time *(int)* – to be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type. The trip_time, as well as the travel_speed, will be used to compose an isochrone map.

>> - travel_speed *(int)* – to be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type. The travel_speed, as well as the trip_time, will be used to compose an isochrone map.

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column for park cover percentage. Dataframe will also be written to new file in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframe
