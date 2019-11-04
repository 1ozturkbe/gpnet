import os

def read_inp(file_directory):
    file = open(file_directory + '.txt', 'r')
    lines = file.readlines()
    file.close()
    # DATA TYPES
    data = {}
    pipe_data = {}
    junction_data = {}
    coordinates = {}
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
        data[marker[n]] = relevant_data
    # CLEANING DATA
    for marker, relevant_data in data.iteritems():
        if marker == j: # junctions
            for i in range(len(relevant_data)):
                junction_index, elevation, demand = relevant_data[i][0:3]
                junction_data[float(junction_index)] = [float(elevation), float(demand)]
        elif marker == r:
            for i in range(len(relevant_data)):
                junction_index, head = relevant_data[i][0:2]
                junction_data[float(junction_index)] = [float(head), 0]
        elif marker == c:
            for i in range(len(relevant_data)):
                junction_index, x, y = relevant_data[i][0:3]
                coordinates[float(junction_index)] = [float(x), float(y)]
        elif marker == p:
            continue
    del data[marker]
    return data

