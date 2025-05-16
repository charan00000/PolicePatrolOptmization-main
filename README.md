# PolicePatrolOptimization

This repository is a supplement to our paper: https://docs.google.com/document/d/10-syP0GjfmwT1vPXm6TVYVB3fopaaWNEwXu1rFjYXdU/edit?usp=sharing. Check requirements.txt for all necessary imports. 

main_driver: this is where you can run all functions. Scroll down past plot(). You can convert a geojson file, as found online on many government platforms, into a graphml file that can be used to find euler circuits, and then convert back into a geojson file. You can visualize the new geojson file using aoftware such as QGIS. When visualizing modified geojsons, make sure to configure labeling so that order labels are displayed, allowing you to see direction in circuits. This label is found as a column in the geo data frame of the output geojson file. 

a graphml file is an XML file that is essentially structured like a graph, containing nodes and edges. THe nodes represent intersections between different roads, and the edges represent small segments of roads thesmelves. The nodes are labeled by their coordinates, which are needed when converting back into geojson. Each edge also stores data contianing the name of the road, and in the future additional properties (speed limits, distance, etc) found on many geojson files for each road segment can be added to the graphml file.

the find_euler_path.modify_graph() function is used to convert a graphml file into one that represnts an euler circuit. The function has three keyword arguments: graphml_input for the input graphml path. 

The third parameter is the method used to eulerize the graph and has three options:
  - "fleury" (default): uses Networkx's builtin eulerize function, which connects vertices of odd degree together with no consideration of weight or distance between odd degreed nodes using the fleury algorithm. In the case the graph is already euler, you can leave this parameter at "base".
  - "min_weights": connects degreed vertices, this time attempting to find, for each odd-degree vertex, another one that's the shortest distance away to connect with. With this method, shortest = lowest number of other vertices crossed and doesn't consider lengths of each segment. Keep in mind that a segment is not an entire road, but a small portion of a road (not necessarily connecting a road intersection to another road intersection).
  - "dijkstra": similar approach to "min_weights" but uses djkstra's algorithm to factor in weights, which in the context of this project is lengths of road segments in (by default) miles. Close in outputted graph to fleury.

The best ways to visualize each graph are:
  - display it using QGIS with either both headings and order columns as labels, or as orders labeled and the geojson file pulled up alongside the QGIS window, reading the edges in order and paying attention to the headings.
    - Headings are in degrees (0-360), with 0 = due north, 90 = due east, 180 = due south, and 270 = due west
   
  - (easiest) use plot() function in main method, and set time_delay keyword parameter = around .01
      - This way currently doesn't show the underlying map of the earth, making it hard to visualize road routes if you're not familiar with the area. You can visualize the graph, without arrows, on the github repository website by opening output_geojson.geojson in order to help.

basic.geojson is a simple geojson representing a graph that's already in euler circuit form. When using it, make sure to set convert_to_graph_road_edges(has_road_names = False), as this graph doesnt have labeled road names.

basic_not_euler.geojson is similar to basic.geojson but with one leaf node. 

find_euler_circuit.modify() returns total distance of euelr circuit, which is used to determine eulerization effectiveness

things to know:
  - you can account for busier roads by duplicating edges represented by busy roads. This is done by convert_to_graph_road_edges()'s keyword parameter weighted_by_road_type. set it to True
      - If you supply your own geojson data, look through it and find what its keys are for road type and road name (find how it labels road name and road type). Then supply them in as keyword parameters to convert_to_graph_road_edges()
  - If you supply your own geojson data, make sure there are no multi line strings in it. These are not supported because they can result in edges not being connected.
      - If you supply your own data and it has multilines, one possible solution is to go back to where you downloaded the data and filter out neighborhood roads.
      - If you end up with a disconnected graph after filtering out multilines and get your geojson from a GIS website, try adjusting the area in which you factor in data before downloading it into a geojson.
