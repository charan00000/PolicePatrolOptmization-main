import ast
import math
import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import json
import lxml
import osmnx as ox
from find_euler_path import calculate_distance_raw
from shapely.geometry import MultiLineString, LineString


def convert_to_graph_road_nodes(geojson_file, dest='new_graph.graphml'):
    """
    inactive
    """
    G = nx.Graph()
    with open(geojson_file, 'r') as f:
        geojson_data = json.load(f)

    node_count = 0
    for feature in geojson_data['features']:
        coordinates = feature['geometry']['coordinates']
        coordinates_str = ",".join(map(str, coordinates))
        G.add_node(node_count, coordinates=coordinates_str)
        node_count += 1

    nx.write_graphml(G, dest)


def convert_to_graph_road_edges(geojson_file, dest='new_graph.graphml',
                                formatted_road_name='FullStName',
                                formatted_road_type='MapClass',
                                has_properties=True,
                                length_unit="Miles",
                                weighted_by_road_type=True):
    """
    Converts a GeoJSON file containing road data into a NetworkX graph with road edges.

    Parameters:
    - geojson_file (str): The path to the GeoJSON file.
    - dest (str): The destination path to save the resulting graph file (default: 'new_graph.graphml').
    - formatted_road_name (str): The name of the road property in the GeoJSON file (default: 'FullStName').
    - has_properties (bool): Indicates whether the GeoJSON file has road properties (default: True).
    - length_unit (str): The unit of length for calculating road distances (default: 'Miles').

    Returns:
    - None

    """
    gdf = gpd.read_file(geojson_file)
    total_distance = 0
    G = nx.MultiGraph()
    for _, road in gdf.iterrows():
        geometry = road.geometry
        if isinstance(geometry, LineString):
            line = [geometry]
        elif isinstance(geometry, MultiLineString):
            raise TypeError("MultiLineString found in the data.")
        for linestring in line:
            for source, target in zip(list(linestring.coords[:-1]), list(linestring.coords[1:])):
                if has_properties:
                    rd_name = road[formatted_road_name]
                    rd_type = road[formatted_road_type]
                else:
                    rd_name = "unnamed"
                    rd_type = 'no_type'
                if rd_type is None:
                    rd_type = 'no_type'
                distance = calculate_distance_raw(source[0],
                                                  source[1],
                                                  target[0],
                                                  target[1],
                                                  in_init_length_unit=length_unit)
                total_distance += distance
                multiplier = find_multiplier(rd_type, formatted_road_type)
                for _ in range(multiplier):
                    G.add_edge(source,
                               target,
                               name=rd_name,
                               type=rd_type,
                               length=distance)
    G.graph['total_distance'] = total_distance
    nx.write_graphml(G, dest)


def convert_to_geojson(graphml_file, dest='output_geojson.geojson', formatted_road_name='FullStName'):
    """
    Convert a GraphML file to GeoJSON format.

    Args:
        graphml_file (str): The path to the GraphML file.
        dest (str, optional): The destination path for the GeoJSON output file. Defaults to 'output_geojson.geojson'.
        formatted_road_name (str, optional): The column name for the road name in the GeoJSON file. Defaults to 'FullStName'.

    Returns:
        None
    """
    G = nx.read_graphml(graphml_file)
    gdf = gpd.GeoDataFrame(columns=['order', formatted_road_name, "length", 'heading', 'road_type', 'geometry'])
    order = 0  # count for each road to be taken to follow eulerian path. Later used to label each road
    first_road = list(G.edges(data=True))[0]
    previous_road_name = first_road[2]['name']
    for source, target, data in G.edges(data=True):
        source = ast.literal_eval(source)
        target = ast.literal_eval(target)
        if data['name'] != previous_road_name:
            order += 1
        previous_road_name = data['name']
        heading = find_heading(source, target)
        new_road = gpd.GeoDataFrame({
            'order': [str(order) + ", "],
            formatted_road_name: [data['name']],
            # 'FullStName' is the column name for the road name in the geojson file
            'length': [data['length']],  # 'Miles' is the column name for the road length in the geojson file
            'heading': [heading],
            'road_type': [data['type']],
            'geometry': [LineString([source, target])]
        })
        gdf = pd.concat([gdf, new_road], ignore_index=True)
    gdf.to_file(dest, driver='GeoJSON')


def find_heading(source, target):
    """
    Calculates the bearing between two points.
    The formula used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    :Parameters:
      - `source: The tuple representing the longitude/latitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `target: The tuple representing the long/lat for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    """
    if (type(source) != tuple) or (type(target) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(source[1])
    lat2 = math.radians(target[1])

    diffLong = math.radians(target[0] - source[0])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
                                           * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # normalize the initial bearing to be from 0 to 360 degrees
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def find_multiplier(rd_type, rd_format):
    if rd_format == "MapClass":
        triple = ["Limited Access Freeway"]
        double = ["Major Rd"]
        if rd_type in triple:
            return 3
        if rd_type in double:
            return 2
    elif rd_format == "RoadPosTyp":
        triple = ["Highway"]
        double = ["Parkway", "Boulevard"]
        if rd_type in triple:
            return 3
        if rd_type in double:
            return 2
    return 1
