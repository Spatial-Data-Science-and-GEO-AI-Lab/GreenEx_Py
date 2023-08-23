**get_shortest_distance_greenspace(point_of_interest, crs_epsg=None, target_dist=300, greenspace_vector_file=None, distance_type="euclidean", destination="centroids", network_type=None, min_greenspace_area=None, plot_aoi=True, write_to_file=True, output_dir=os.getcwd())**

> Assess whether or not greenspaces are located within a threshold distance for points of interest.

> Function to calculate the shortest distance from point locations to greenspaces in a euclidean or network manner. The greenspace destination points can be set to either "centroids" or "entrance". If "entrance", the distance will be calculated from the point location to the network nodes which are within 20 meters of the greenspace boundaries - also referred to as pseudo greenspace entry points. 

>> Parameters: 

>> - point_of_interest *(string)* – either the absolute/relative path to the file or the geodataframe containing point or (multi)polygon geometries for which to compute the shortest distance to greenspaces.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest data has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of data is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - target_dist *(int)* – threshold distance in meters surrounding the point locations in which greenspaces should be located. If no greenspace destination is detected within threshold distance for the current parameters, the distance to greenspace will be set to the given threshold distance. The actual distance to the nearest greenspace cannot always be retrieved properly due to edge effects concerning the area for which the network and greenspaces are retrieved through OpenStreetMap. A warning is also provided in this case, as it should be considered for further analysis. 

>> - greenspace_vector_file *(string)* – optional, the absolute or relative path to the vector file containing greenspace data. In case no file is provided, greenspace data will be retrieved from OpenStreetMap.

>> - distance_type *(string {"euclidean", "network"})* – the way in which the shortest distance to a greenspace should be calculated, straight line distance (euclidean) or network distance (network). By default, argument is set to euclidean. Note that if set to network, processing times might increase significantly in case many points/polygons of interest are provided in the point_of_interest data.

>> - destination *(string {"centroids", "entrance"})* – the destination points for the greenspaces. If "entrance", distances will be computed between the point location and the network nodes that are within 20 meters of the greenspace boundaries - therefore considered as pseudo-entry points.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which the network needs to be retrieved from OpenStreetMap.

>> - min_greenspace_area *(int)* – optional, if set to integer value, only greenspaces with an area greater than or equal to the min_greenspace_area in squared meters will be considered for the shortest distance calculation. 

>> - plot_aoi *(bool {"TRUE", "FALSE"})* - whether or not to plot the areas of interest that have been used for the shortest distance calculation. If set to TRUE (default), the plot will be shown as part of the function execution.

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest including column indicating whether or not a greenspace is within specified target distance. Dataframe will also be written to new file in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframe
