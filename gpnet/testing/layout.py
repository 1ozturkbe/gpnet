import numpy as np
from known_topology_network_model import HW_KT_FND, DW_KT_FND
from generate_data import define_length

def hanoi(friction='DW'):
    # Hanoi layout problem from University of Exeter Centre for Water Systems
    # http://emps.exeter.ac.uk/engineering/research/cws/downloads/benchmarks/layout/
    N = 32
    sinks = np.array([0, 890, 850, 130, 725, 1005, 1350, 550, 525,
             525, 500, 560, 940, 615, 280, 310,
             865, 1345, 60, 1275, 930, 485, 1045, 820, 170,
             900, 370, 290, 360, 360, 105, 805])/3600.0
    sources = [0 for i in range(N)]
    sources[0] = sum(sinks)
    topology_list = [[0, 1], [1, 2], [2, 3], [2, 19], [2, 18], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                     [10, 11], [11, 12], [9, 13], [13, 14], [14, 15], [16, 15], [17, 16], [18, 17], [15, 26], [26, 25],
                     [25, 24], [23, 24], [22, 23], [22, 27], [27, 28], [28, 29], [29, 30], [24, 31], [31, 30], [19, 20],
                     [20, 21], [19, 22]]
    coordinates = {1: (5021.20, 1582.17), 2: (5025.20, 2585.42), 3: (5874.22, 2588.30), 4: (6873.11, 2588.30),
                   5: (8103.51, 2585.42), 6: (8103.51, 3234.67), 7: (8106.66, 4179.28), 8: (8106.66, 5133.78),
                   9: (7318.64, 5133.78), 10: (7319.94, 5831.65), 11: (7319.94, 6671.19), 12: (5636.76, 6676.24),
                   13: (6530.63, 5133.78), 14: (5676.02, 5133.78), 15: (5021.20, 5133.78), 16: (5021.20, 4412.36),
                   17: (5021.20, 3868.52), 18: (5021.20, 3191.49), 19: (3587.87, 2588.30), 20: (3587.87, 1300.84),
                   21: (3587.87, 901.29), 22: (1978.55, 2588.30), 23: (1975.58, 4084.35), 24: (1980.46, 5137.63),
                   25: (3077.46, 5137.63), 26: (3933.52, 5133.78), 27: (846.04, 2588.20), 28: (-552.41, 2588.20),
                   29: (-552.38, 4369.06), 30: (-549.36, 5137.63), 31: (536.45, 5137.63), 0: (5360.71, 1354.05)}
    # Finding lengths
    L_all = define_length(coordinates)
    L = [L_all[pipe[0], pipe[1]] for pipe in topology_list]
    # It seems like lengths claimed by the data are wrong...(?)
    # L = [100, 1350, 900, 1150, 1450, 450, 850, 850, 800, 950, 1200, 3500, 800, 500, 550, 2730, 1750, 800, 400, 2200,
    #      1500, 500, 2650, 1230, 1300, 850, 300, 750, 1500, 2000, 1600, 150, 860, 950]

    # Defining roughness and minimum head (m)
    h_min = [30 for _ in range(N)]

    if friction == 'DW':
        m = DW_KT_FND(N, topology_list)
        roughness = 0.26e-6
    elif friction == 'HW':
        m = HW_KT_FND(N, topology_list)
        roughness = 130
    else:
        print('Friction model %s is not yet supported.' % friction)

    m.substitutions.update({
        "L": L,
        "\dot{V}_+": sources,
        "\dot{V}_-": sinks,
        "\\epsilon": roughness,
        "H_{min}": h_min,
        m["H"][0]: 100,
        "D_{max}": 1.016,
        "D_{min}": 0.3048,
        "F_{max}": 1e20,
    })
    # Substitutions depending on friction model
    if friction == 'DW':
        m.substitutions.update({
            "\\rho": 1000,
            "\\mu": 8.9e-4,
            "g": 9.81,
        })

    m.cost = m['C']
    m.coordinates = coordinates
    m.topology_list = topology_list
    return m

def small_graph():
    N = 5
    sinks = [0, 0.89, 0.85, 0.130, 0.725]
    sources = [sum(sinks), 0, 0, 0, 0]
    topology_list = [[0, 1], [1, 2], [1, 3], [2, 4] ,[3,4]]
    coordinates = {0: (0, 1000), 1: (0, 0), 2: (1000, 0), 3: (-1000, 0), 4: (1000, -1000)}
    return sources, sinks, topology_list, coordinates
