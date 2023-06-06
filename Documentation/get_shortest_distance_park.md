**get_shortest_distance_park(point_of_interest_file, crs_epsg=None, target_dist=300, park_vector_file=None, distance_type="euclidian", destination="centroids", network_type=None, write_to_file=True, output_dir=os.getcwd())**

> Assess whether or not parks are located within a threshold network distance for points of interest.

> Function to calculate the shortest network distance from point locations to parks. The park destination points can be set to either "centroids" or "entrance". The network distance will be calculated from the point location's nearest network node to the network nodes which are within 20 meters of the park boundaries - also referred to as fake park entry points. The euclidian distance between the point location and the nearest network node will be added to minimize the error rate. In case of destination = "centroids", the euclidian distance between the park's fake entry points and the park's centroid will also be considered.

>> Parameters: 

>> - point_of_interest_file *(string)* – the absolute or relative path to the file containing point or polygon geometries for which to compute the shortest distance to parks.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest file has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of file is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - target_dist *(int)* – threshold distance surrounding the point locations in which parks should be located.

>> - park_vector_file *(string)* – optional, the absolute or relative path to the vector file containing park data. In case no file is provided, park data will be retrieved from OpenStreetMap.

>> - distance_type *(string {"euclidian", "network"})* – the way in which the shortest distance to a park should be calculated, straight line distance (euclidian) or network distance (network). By default, argument is set to euclidian. Note that if set to network, processing times might increase significantly in case many points/polygons of interest are provided in the point_of_interest file.

>> - destination *(string {"centroids", "entrance"})* – the destination points for the parks. If "entrance", network distances will be computed between the point location's nearest network node and the network nodes that are within 20 meters of the park boundaries - therefore considered as (fake) entry points. Euclidian distance between point location and nearest network node will be added. If "centroids", the same procedure applies while also taking into account the euclidian distance between the park's centroid and fake entry points.

>> - network_type *(string {"walk", "bike", "drive", "all"})* – to be defined in case buffer_type is set to "network", the travel mode for which the network needs to be retrieved from OpenStreetMap.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest_file including column indicating whether or not a park is within specified target distance. Dataframe will also be written to new file in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframe
