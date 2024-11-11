import ast
import networkx as nx
import matplotlib.pyplot as plt
import heapq
from itertools import combinations
from pyproj import Geod
from networkx.utils import pairwise


def modify_graph(graphml_input='new_graph.graphml',
                 dest='euler_path_output.graphml',
                 euler_form_method="built_in",
                 euler_order_method="built_in",
                 length_unit="miles"):
    """
    Modifies a graph by finding an Euler path and writing the modified graph to a GraphML file.

    Parameters:
    - graphml_input (str): Path to the input GraphML file.
    - dest (str): Path to the output GraphML file.
    - method (str): Method for Eulerization. Can be "built_in" or "min_weights" or "dijkstra".
    - length_unit (str): Unit of length for calculating distances.

    Returns:
    - list: A list containing the total distance of the Euler path and the number of artificial edges created.
    """
    G = nx.read_graphml(graphml_input)
    if euler_form_method == "min_weights":
        euler_G = eulerize_minimize_weights(G)
    elif euler_form_method == "dijkstra":
        euler_G = eulerize_minimize_weights_dijkistra(G)
    elif euler_form_method == "built_in_weighted":
        euler_G = eulerize_built_in_weighted(G)
    else:
        euler_G = eulerize_built_in(G)
    nx.write_graphml(euler_G, 'eulerized_graph.graphml')
    if euler_order_method == "trotter":
        circuit = trotter(euler_G)
    else:
        circuit = list(nx.eulerian_circuit(euler_G))
    nx.write_graphml(nx.MultiDiGraph(circuit), dest)
    new_G = nx.MultiDiGraph()
    total_distance = 0
    artificial_edges = 0
    for source, target in circuit:
        edge_data = euler_G.get_edge_data(source, target)
        if edge_data is None or len(edge_data[0]) == 0 or len(edge_data) == 0 or 'name' not in edge_data[0]:
            road_name = "unnamed"
            artificial_edges += 1
            init_length = calculate_distance(source, target, init_length_unit=length_unit)
        else:
            road_name = edge_data[0]['name']
            init_length = edge_data[0]['length']
            init_type = edge_data[0]['type']
        new_G.add_edge(source, target, name=road_name, length=init_length, type=init_type)
        total_distance += init_length

    nx.write_graphml(new_G, dest)
    old_length = G.graph['total_distance']
    circuit_length_multiplier = total_distance / old_length
    return [total_distance, old_length, circuit_length_multiplier, "artificial edges: " + str(artificial_edges)]


def eulerize_built_in(G):
    """
    Eulerizes a graph by adding edges.
    Doesn't add an edge between nodes that dont already have a single edge between them

    Parameters:
    G (networkx.Graph): The input graph.

    Returns:
    networkx.Graph: The Eulerized graph.
    """
    return nx.eulerize(G)


def eulerize_built_in_weighted(G):
    """
    WORK IN PROGRESS NOT READY YET

    Mostly contains source code from networkx eulerize() function with few tweaks to account for road segment length.
    Doesn't add an edge between nodes that dont already have a single edge between them

    Parameters:
    G (networkx.Graph): The input graph.

    Returns:
    networkx.Graph: The Eulerized graph.

    References:
    - Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and
      function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008),
      Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008
    """
    if G.order() == 0:
        raise nx.NetworkXPointlessConcept("Cannot Eulerize null graph")
    if not nx.is_connected(G):
        raise nx.NetworkXError("G is not connected")
    odd_degree_nodes = [n for n, d in G.degree() if d % 2 == 1]
    G = nx.MultiGraph(G)
    if len(odd_degree_nodes) == 0:
        return G

    # get all shortest paths between vertices of odd degree
    odd_deg_pairs_paths = []
    for m, n in combinations(odd_degree_nodes, 2):
        path = nx.dijkstra_path(G, source=m, target=n, weight=(lambda u,v, data:  G[u][v][0].get('weight', 1)))
        odd_deg_pairs_paths.append((m, {n: path}))

    # use the number of vertices in a graph + 1 as an upper bound on
    # the maximum length of a path in G
    upper_bound_on_max_path_length = len(G) + 1

    # use "len(G) + 1 - len(P)",
    # where P is a shortest path between vertices n and m,
    # as edge-weights in a new graph
    # store the paths in the graph for easy indexing later
    Gp = nx.Graph()
    for n, Ps in odd_deg_pairs_paths:
        for m, P in Ps.items():
            if n != m:
                Gp.add_edge(
                    m, n, weight=upper_bound_on_max_path_length - len(P), path=P
                )

    # find the minimum weight matching of edges in the weighted graph
    best_matching = nx.Graph(list(nx.max_weight_matching(Gp)))

    # duplicate each edge along each path in the set of paths in Gp
    for m, n in best_matching.edges():
        path = Gp[m][n]["path"]
        G.add_edges_from(nx.utils.pairwise(path))
    return G


def eulerize_minimize_weights(old_G):
    """
    Eulerize the given graph by adding edges between pairs of odd-degree nodes
    to minimize the weights of the resulting Eulerian circuit.

    Parameters:
    old_G (networkx.Graph): The input graph.

    Returns:
    networkx.Graph: The Eulerized graph.
    """
    G = old_G.copy()
    print("started eulerize_minimize_weights")
    odd_degree_nodes = [node for node, degree in G.degree() if degree % 2 == 1]
    shortest_paths = dict(nx.all_pairs_shortest_path(G))
    # Create a priority queue of pairs of odd-degree nodes, with distances as priorities
    pair_queue = []
    for i, node_A in enumerate(odd_degree_nodes):
        for node_B in odd_degree_nodes[i + 1:]:
            path_length = len(shortest_paths[node_A][node_B])
            pair_queue.append((path_length, node_A, node_B))
        heapq.heapify(pair_queue)
    # While there are nodes with odd degree
    while pair_queue:
        # Pop the pair with the shortest distance
        _, node_A, node_B = heapq.heappop(pair_queue)

        if node_A in odd_degree_nodes and node_B in odd_degree_nodes:
            # Duplicate all edges in the shortest path between node1 and node2
            path = shortest_paths[node_A][node_B]
            for i in range(len(path) - 1):
                G.add_edge(path[i], path[i + 1])
            # Remove node1 and node2 from the list of nodes with odd degree
            odd_degree_nodes.remove(node_A)
            odd_degree_nodes.remove(node_B)
    odd_degree_nodes = [node for node, degree in G.degree() if degree % 2 == 1]
    print("odd degree nodes: ", odd_degree_nodes)
    return nx.eulerize(G)


def eulerize_minimize_weights_dijkistra(old_G):
    """
    Eulerize the given graph by adding edges between pairs of odd-degree nodes
    to minimize the weights of the resulting Eulerian circuit.

    Parameters:
    old_G (networkx.Graph): The input graph.

    Returns:
    networkx.Graph: The Eulerized graph.
    """
    G = old_G.copy()
    print("started eulerize_minimize_weights")
    odd_degree_nodes = [node for node, degree in G.degree() if degree % 2 == 1]

    shortest_paths = dict(nx.all_pairs_dijkstra_path(G, weight='length'))
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G, weight='length'))

    # Create a priority queue of pairs of odd-degree nodes, with distances as priorities
    pair_queue = [(path_lengths[node_A][node_B], node_A, node_B) for i, node_A in enumerate(odd_degree_nodes) for node_B
                  in odd_degree_nodes[i + 1:]]
    heapq.heapify(pair_queue)
    while pair_queue:

        # Pop the pair with the shortest distance
        _, node_A, node_B = heapq.heappop(pair_queue)
        if node_A in odd_degree_nodes and node_B in odd_degree_nodes:

            # Duplicate all edges in the shortest path between node1 and node2
            path = shortest_paths[node_A][node_B]
            for i in range(len(path) - 1):
                G.add_edge(path[i], path[i + 1], length=G[path[i]][path[i + 1]][0]['length'])
                print(G[path[i]][path[i + 1]])
            odd_degree_nodes.remove(node_A)
            odd_degree_nodes.remove(node_B)
    odd_degree_nodes = [node for node, degree in G.degree() if degree % 2 == 1]
    print("odd degree nodes: ", odd_degree_nodes)
    return nx.eulerize(G)


def calculate_distance_raw(lon1, lat1, lon2, lat2, in_init_length_unit="miles"):
    """
    Calculate the distance between two points on the Earth's surface using longitude and latitude coordinates.

    Parameters:
    lon1 (float): The longitude of the first point.
    lat1 (float): The latitude of the first point.
    lon2 (float): The longitude of the second point.
    lat2 (float): The latitude of the second point.
    in_init_length_unit (str, optional): The unit of length for the calculated distance. Default is "miles".

    Returns:
    float: The calculated distance between the two points.

    """
    geod = Geod(ellps='WGS84')
    _, _, distance = geod.inv(lon1, lat1, lon2, lat2)
    if in_init_length_unit == "kilometers":
        return distance / 1000
    return distance / 1609.344


def calculate_distance(source, target, init_length_unit="miles"):
    """
    Calculate the distance between two points.

    Args:
        source (str): The coordinates of the source point in the format "(latitude, longitude)".
        target (str): The coordinates of the target point in the format "(latitude, longitude)".
        init_length_unit (str, optional): The initial length unit. Defaults to "miles".

    Returns:
        float: The calculated distance between the source and target points.
    """
    tup_source = ast.literal_eval(source)
    tup_target = ast.literal_eval(target)
    return calculate_distance_raw(tup_source[0],
                                  tup_source[1],
                                  tup_target[0],
                                  tup_target[1],
                                  in_init_length_unit=init_length_unit)


def trotter(G):
    """
    Trotter's algorithm from the lecture videos for finding an Euler path in a graph.

    Parameters:
    G (networkx.Graph): The input graph.

    Returns:
    list: A list containing the nodes in the Euler path.
    """
    G = G.copy()
    current_node = list(G.nodes)[0]
    total_circuit = [current_node]
    current_path = [current_node]

    while current_path:
        if G[current_node]:
            # Select the neighbor with the lowest number and remove edge between them
            next_node = min(G.neighbors(current_node))
            G.remove_edge(current_node, next_node)
            current_node = next_node
            current_path.append(current_node)
        else:
            # If the current node has no neighbors, add it to the Euler circuit
            total_circuit.insert(total_circuit.index(current_path[0]) + 1, current_node)
            current_node = current_path.pop()

    euler_circuit_edge_list = list(pairwise(total_circuit))
    return euler_circuit_edge_list