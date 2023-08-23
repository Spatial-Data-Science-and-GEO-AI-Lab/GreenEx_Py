**get_streetview_GVI(point_of_interest, access_token=None, crs_epsg=None, polygon_type="neighbourhood", buffer_dist=None, workers=4, crop_by_road_centres=1, network_file=None, write_to_file=True, output_dir=os.getcwd())**

> Retrieve the average Greenness Visibility Index (GVI), based on a streetview analysis, for points and/or areas of interest. Note that this function is based on research conducted by [Ilse A. Vázquez Sánchez](https://github.com/Spatial-Data-Science-and-GEO-AI-Lab/StreetView-NatureVisibility).

> Function to calculate the average GVI based on a streetview analysis. The function accepts a file containing point or polygon geometries. Based on the user's specifications, GVI values will be calculated for multiple road network locations that are close to or within the original points/polygons of interest. In the end, these values will be averaged to come up with a final score. Also, the number of points upon which this score is based, will be included in the results of the function. To calculate the GVI, streetview images are obtained from Mapillary for which an API token is required.

>> Parameters: 

>> - point_of_interest *(string)* – either the absolute/relative path to the file or the geodataframe containing point or (multi)polygon geometries around and for which to compute the average GVI value.

>> - access_token *(string)* – Mapillary API token that is required to access the streetview images of Mapillary.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest data has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of data is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - polygon_type *(string {"neighbourhood", "house"})* - to be defined in case point_of_interest contains polygon geometries. In case set to "neighbourhood", buffer_dist argument is optional and if not specified, GVI values will be calculated for road network locations within each polygon geometry. If set to "house", network will be retrieved from OpenStreetMap and based on buffer_dist if not provided. 

>> - buffer_dist *(int)* – to be defined in case point_of_interest contains point geometries and optional in case point_of_interest contains polygon geometries, the point/polygon of interest surrounding distance in meters that should be considered when defining road network locations for calculating GVI values. 

>> - workers *(int)* – the maximum number of concurrent worker threads to use (i.e. simultaneous tasks that can be executed), providing control over the level of concurrency in the function.

>> - crop_by_road_centres *(bool {"TRUE", "FALSE"})* – the way in which panoramic streetview images should be cropped. If set to False, panoramic images will be cropped into 4 equal parts. If set to True (default), panoramic images are cropped based on road centres.

>> - network_file *(string)* – optional, the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified, the network will be retrieved through OpenStreetMap considering the point/polygon's surrounding straight line distance as specified in the buffer_dist argument.

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest including columns for mean GVI value(s) and the number of points upon which these are based. A second dataframe will be returned and contains the surrounding road network locations (of the original points/polygons of interest) that were used to calculate the GVI values. The ID column is used to identify the original point/polygon of interest. Both dataframes will also be written to new files in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframes