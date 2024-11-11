import conversions
import contextily as ctx
import find_euler_path
import geopandas as gpd
import matplotlib.pyplot as plt


# scroll down past this function

def plot(geojson_file, time_delay=0, arrow_spacing=15):
    """
    Plots the lines from a GeoJSON file and adds arrows to represent the direction of the lines.
    Supports delay while graphing to better visualize the direction of path

    Parameters:
    - geojson_file (str): The path to the GeoJSON file.
    - time_delay (float): The delay between each line plot in seconds. Default is 0.
    - arrow_spacing (int): The spacing between arrows. Only every `arrow_spacing` line will have an arrow. Default is 15.

    Returns:
    None
    """
    # Load the GeoJSON file
    gdf = gpd.read_file(geojson_file)

    # Create a new plot
    fig, ax = plt.subplots()

    # Iterate over the GeoDataFrame
    for i, row in gdf.iterrows():
        # Get the line's start and end points
        if row['geometry'].geom_type == 'MultiLineString' or row['geometry'].geom_type == 'MultiPolygon':
            start = row['geometry'][0].coords[0][:2]
            end = row['geometry'][-1].coords[-1][:2]
        else:
            start = row['geometry'].coords[0][:2]
            end = row['geometry'].coords[-1][:2]

        # Plot the line
        ax.plot(*row['geometry'].xy, color='blue')

        # Add an arrow at the end of the line to represent the direction
        if i % arrow_spacing == 0 and time_delay == 0:
            ax.annotate('', xy=end, xytext=start, arrowprops=dict(facecolor='red', edgecolor='red'))

        # fig.canvas.draw()

        # Add a delay
        if time_delay > 0:
            plt.pause(time_delay)

    # Show the plot
    # ctx.add_basemap(ax)
    plt.show()


"""
The lines below are the important ones. Here, you can run all the code from the rest
of the repository.

convert_to_graph_road_edges() converts a GeoJSON file to a NetworkX graph with road edges.
If you use your own geojson data, and its not from forsyth county ga, you will need to verify
that the geojson data is in the correct format and enter the appropriate label for road type
and road name as keyword parameters below. If you are using data that doesn't have names for 
each road, make sure to set has_properties to False. You can also set weighted_by_road_type 
to false if you want to minimize the length of the whole path and don't want certain roads to 
be traversed more times than others. 

- If weighted_by_road_type is True, the multipliers will default to being setup as can be seen
in "conversions.py". You can modify the find_multiplier() method there according to your
original data's labeling. 

find_euler_path.modify_graph() takes a graphml file and modifies it to be eulerian. The parameters 
are defaulted so that it converts the file produced by convert_to_graph_road_edges(). There are three
options for the method parameter: "built_in" and "min_weights". "built_in" is faster and more practical
for most applications.
will produce a shorter eulerian path.
- this returns a list of three numbers in the from of 
[(total distance of eulerian path), (total distance of original path), (eulerization distance increase multiplier), (number of artificial edges added)]
this output is later printed to the console.

- You can change the method to find an euler circuit, although the differences are minimal.
    - built_in: uses .eulerian_circuit() function from networkX
    - "trotter": uses algorithm taught in class
- You can run this program as is, with no modifications to the parameters below, for an example eulerization and output visualized.
"""

conversions.convert_to_graph_road_edges('forsyth_major_bottom_left_roads.geojson',
                                        dest='forsyth_major_bottom_left_roads.graphml',
                                        formatted_road_name='FullStName',  # road name label
                                        formatted_road_type='RoadPosTyp',  # road type label
                                        has_properties=True,
                                        length_unit='miles',
                                        weighted_by_road_type=True)  # toggle for multiplying busy roads

attributes = find_euler_path.modify_graph(graphml_input='forsyth_major_bottom_left_roads.graphml',
                                          dest='euler_path_output.graphml',
                                          euler_form_method="built_in",     # method to produce graph capable of forming euler circuit
                                          euler_order_method="trotter",    # method to order euler circuit
                                          length_unit="miles")

conversions.convert_to_geojson('euler_path_output.graphml')
print(attributes)
print("starting pop-up mapping")
plot('output_geojson.geojson', time_delay=0.00000000000000000001)
print('finished pop-up mapping')
