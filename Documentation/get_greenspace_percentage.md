**get_greenspace_percentage(point_of_interest, greenspace_vector_file=None, crs_epsg=None, polygon_type="neighbourhood", buffer_type=None, buffer_dist=None, network_type=None, trip_time=None, travel_speed=None, plot_aoi=True, write_to_file=True, output_dir=os.getcwd())**

> Retrieve the percentage of area covered by greenspaces for areas or points of interest.

> Function to calculate the percentage of area covered by greenspaces surrounding geographic locations or within geographic areas. The surroundings can be defined by a euclidean (straight-line) or network (travel distance considering transportation infrastructure) buffer distance. 

>> Parameters: 

>> - point_of_interest *(string)* – either the absolute/relative path to the file or the geodataframe containing point or (multi)polygon geometries around and for which to compute the percentage of greenspace cover.

>> - greenspace_vector_file *(string)* – optional, the absolute or relative path to the vector file containing greenspace data. Note that geometries should be polygon or multipolygon. In case no file is provided, greenspace data will be retrieved from OpenStreetMap.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest data has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of data is geographic, CRS will be transformed to EPSG 3395 by default.

>> - polygon_type *(string {"neighbourhood", "house"})* - to be defined in case point_of_interest contains polygon geometries. If set to "house", polygon geometries will be converted to point geometries.

>> - buffer_type *(string {"euclidean", "network"})* – to be defined in case point_of_interest contains point geometries or polygon geometries which represent houses and optional in case point_of_interest contains polygon geometries which represent neighbourhoods, the way in which the area of interest should be composed. If "euclidean", a straight line distance will be used based on the buffer distance as specified by the buffer_dist argument. If "network", isochrone maps will be composed based on the buffer distance or combination of trip_time and travel_speed.

>> - buffer_dist *(int)* – to be defined if buffer_type is set to "euclidean" and optional if buffer_type is set to "network". The distance in meters that will be used to compute the area of interest surrounding the geometries of the point_of_interest data. NOTE: In case buffer_type is set to "network", buffer_dist should not be specified if travel_speed and trip_time arguments are specified.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which the network needs to be retrieved from OpenStreetMap.

>> - trip_time *(int)* – may optionally be defined in case buffer_type is set to "network", trip time in minutes to consider for travel mode specified in network_type. The trip_time, in addition to the travel_speed, will be used to compose an isochrone map if no buffer distance is provided. NOTE: travel_speed and trip_time arguments should not be specified if buffer_dist is specified.

>> - travel_speed *(float, int)* – may optionally be defined in case buffer_type is set to "network", travel speed in km/h to consider for travel mode specified in network_type. The travel_speed, in addition to the trip_time, will be used to compose an isochrone map if no buffer distance is provided. NOTE: travel_speed and trip_time arguments should not be specified if buffer_dist is specified.

>> - plot_aoi *(bool {"TRUE", "FALSE"})* - whether or not to plot the areas of interest that have been used for the greenspace cover calculation. If set to TRUE (default), the plot will be shown as part of the function execution.

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest including column for greenspace cover percentage. Dataframe will also be written to new file in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframe
