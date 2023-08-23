**get_viewshed_GVI(point_of_interest, greendata_raster_file, dtm_raster_file, dsm_raster_file, network_file=None, crs_epsg=None, polygon_type="neighbourhood", buffer_dist=None, viewing_dist=250, sample_dist=50, observer_height=1.7, write_to_file=True, output_dir=os.getcwd())**

> Retrieve the average Greenness Visibility Index (GVI), based on a viewshed analysis, for points and/or areas of interest. Note that this function is based on research conducted by [Jonny Huck & Labib Labib](https://github.com/jonnyhuck/green-visibility-index/tree/master).

> Function to calculate the average GVI based on a viewshed analysis. The function accepts a file containing point or polygon geometries. Based on the user's specifications, GVI values will be calculated for multiple road network locations that are close to or within the original points/polygons of interest. In the end, these values will be averaged to come up with a final score. Also, the number of points upon which this score is based, will be included in the results of the function. To calculate the GVI, a viewshed is created based on a Digital Surface Model (DSM) and Digital Terrain Model (DTM) for each road network location that was defined in the process. Using this viewshed, the number of visible green pixels, as well as the total number of pixels, is retrieved from the greenspace raster data provided by the user. 

>> Parameters: 

>> - point_of_interest *(string)* – either the absolute/relative path to the file or the geodataframe containing point or (multi)polygon geometries around and for which to compute the average GVI value.

>> - greendata_raster_file *(string)* – the absolute or relative path to the boolean raster file indicating which cells are considered greenspace and which are not.

>> - dtm_raster_file *(string)* – the absolute or relative path to the raster file containing the Digital Terrain Model.

>> - dsm_raster_file *(string)* – the absolute or relative path to the raster file containing the Digital Surface Model.

>> - network_file *(string)* – optional, the absolute or relative path to the file containing the network (transportation infrastructure) to consider. If not specified, the network will be retrieved through OpenStreetMap considering the point/polygon's surrounding straight line distance as specified in the buffer_dist argument.

>> - crs_epsg *(int)* - optional, to be defined in case provided point of interest data has geographic CRS rather than projected. CRS will be transformed to the projected CRS that is specified. In case crs_epsg is not specified and CRS of data is geographic, CRS will be transformed to EPSG 3395 by default. 

>> - polygon_type *(string {"neighbourhood", "house"})* - to be defined in case point_of_interest contains polygon geometries. In case set to "neighbourhood", buffer_dist argument is optional and if not specified, GVI values will be calculated for road network locations within each polygon geometry. If set to "house", network will be retrieved from OpenStreetMap and based on buffer_dist if not provided. 

>> - buffer_dist *(int)* – to be defined in case point_of_interest contains point geometries and optional in case point_of_interest contains polygon geometries, the point/polygon of interest surrounding distance in meters that should be considered when defining road network locations for calculating GVI values. 

>> - viewing_dist *(int)* - the viewing distance in meters to consider when composing the viewshed.

>> - sample_dist *(float, int)* - the interval in meters that is used to define road network locations. For each road that is within the area of interest, point locations will be generated using an interval of sample_dist meters. If a road is too short, the road's center location will be used. 

>> - observer_height *(float, int)* - a person's height in meters to consider for the creation of viewsheds and lines of sight. 

>> - write_to_file *(bool {"TRUE", "FALSE"})* - whether or not to write the results to a new file in the directory specified in the output_dir argument. By default, results will be written to file.

>> - output_dir *(string)* – the absolute or relative path to the directory in which the output file will be written in case write_to_file is set to TRUE. If not specified, the current working directory will serve as default.

>>Returns:	
>>> Dataframe as obtained from point_of_interest including columns for mean GVI value(s) and the number of points upon which these are based. A second dataframe will be returned and contains the surrounding road network locations (of the original points/polygons of interest) that were used to calculate the GVI values. The ID column is used to identify the original point/polygon of interest. Both dataframes will also be written to new files in specified directory (see output_dir argument) if write_to_file set to TRUE. 

>>Return type:	
>>> Geodataframes