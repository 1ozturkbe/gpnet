import os

def read_inp(file_directory):
    file = open(file_directory, 'r')
    lines = file.readlines()
    file.close()
    # DATA TYPES
    data = {}
    elevation = {}
    head = {}
    sinks = {}
    coordinates = {}
    topology_dict = {}
    L = {}
    D = {}
    roughness = {}
    # MARKERS
    j = '[JUNCTIONS]\n'
    r = '[RESERVOIRS]\n'
    p = '[PIPES]\n'
    c = '[COORDINATES]\n'
    markers = [j,r,p,c]
    marker_indices = [lines.index(marker) for marker in markers]
    # GRABBING DATA
    for n, idx in enumerate(marker_indices):
        count = 2
        relevant_data = []
        while lines[idx+count] != '\n':
            relevant_data.append(lines[idx+count].split())
            count += 1
        data[markers[n]] = relevant_data
    # CLEANING DATA
    for marker, relevant_data in data.iteritems():
        if marker == j: # junctions
            for i in range(len(relevant_data)):
                junction_index, elev, demand = relevant_data[i][0:3]
                elevation[float(junction_index)] = float(elev)
                sinks[float(junction_index)] = float(demand)
        elif marker == r:
            for i in range(len(relevant_data)):
                junction_index, headpress = relevant_data[i][0:2]
                head[float(junction_index)] = float(headpress)
                elevation[float(junction_index)] = float(headpress)
        elif marker == c:
            for i in range(len(relevant_data)):
                junction_index, x, y = relevant_data[i][0:3]
                coordinates[float(junction_index)] = [float(x), float(y)]
        elif marker == p:
            for i in range(len(relevant_data)):
                pipe_ID, node1, node2, length, diameter, rough = relevant_data[i][0:6]
                pipe_ID = float(pipe_ID)
                topology_dict[pipe_ID] = [int(node1), int(node2)]
                L[pipe_ID] = float(length)
                D[pipe_ID] = float(diameter)
                roughness[pipe_ID] = float(rough)
    if len(coordinates) != len(elevation):
        raise ValueError('Dimensionality mismatch in coordinate and elevation data!')
    elif len(sinks) + len(head) != len(coordinates):
        raise ValueError('Problem is over/underdefined in demand and head pressure!')
    return elevation, head, sinks, coordinates, topology_dict, L, D, roughness

