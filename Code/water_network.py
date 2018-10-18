from general_incompressible_network_model import IncompressibleFluidNetworkDistribution
from known_topology_network_model import KTFND
from generate_data import define_topology, define_length

if __name__ == '__main__':
    # Somewhat large problem
    # N = 32
    # sinks = [0, 890, 850, 130, 725, 1005, 1350, 550, 525, 525, 500, 560, 940, 615, 280, 310,
    #          865, 1345, 60, 1275, 930, 485, 1045, 820, 170, 900, 370, 290, 360, 360, 105, 805]
    # sources = [0 for i in xrange(N)]  # not correct like this
    # sources[0] = sum(sinks)
    # topology_list = [[0, 1], [1, 2], [2, 3], [2, 19], [2, 18], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
    #                  [10, 11], [11, 12], [9, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [15, 26], [26, 25],
    #                  [25, 24], [23, 24], [22, 23], [22, 27], [27, 28], [28, 29], [29, 30], [24, 31], [31, 30], [19, 20],
    #                  [20, 21], [19, 22]]
    # coordinates = {2: (5021.20, 1582.17), 3: (5025.20, 2585.42), 4: (5874.22, 2588.30), 5: (6873.11, 2588.30),
    #                6: (8103.51, 2585.42), 7: (8103.51, 3234.67), 8: (8106.66, 4179.28), 9: (8106.66, 5133.78),
    #                10: (7318.64, 5133.78), 11: (7319.94, 5831.65), 12: (7319.94, 6671.19), 13: (5636.76, 6676.24),
    #                14: (6530.63, 5133.78), 15: (5676.02, 5133.78), 16: (5021.20, 5133.78), 17: (5021.20, 4412.36),
    #                18: (5021.20, 3868.52), 19: (5021.20, 3191.49), 20: (3587.87, 2588.30), 21: (3587.87, 1300.84),
    #                22: (3587.87, 901.29), 23: (1978.55, 2588.30), 24: (1975.58, 4084.35), 25: (1980.46, 5137.63),
    #                26: (3077.46, 5137.63), 27: (3933.52, 5133.78), 28: (846.04, 2588.20), 29: (-552.41, 2588.20),
    #                30: (-552.38, 4369.06), 31: (-549.36, 5137.63), 32: (536.45, 5137.63), 1: (5360.71, 1354.05)}
    # Small problem
    # N = 5
    # sinks = [0,50,100,100,50]
    # sources = [sum(sinks), 0,0,0,0]
    # topology_list = [[0,1],[1,2],[1,3],[2,4],[3,4]]
    # coordinates = {0: (0,1000), 1: (0,0), 2: (1000,0), 3:(-1000,0), 4:(1000,-1000)}
    # Small tree problem
    N = 3
    sinks = [0,5,10]
    sources = [sum(sinks), 0,0]
    topology_list = [[0,1],[1,2]]
    coordinates = {0: (0,1000), 1: (0,0), 2:(-1000,0)}
    connect = define_topology(topology_list,N)
    L_all = define_length(coordinates)
    L = [L_all[pipe[0],pipe[1]] for pipe in topology_list]
    roughness = [[0.26e-6 for _ in xrange(N)] for _ in xrange(N)]
    h_min = [30 for _ in xrange(N)]

    water_distribution = KTFND(N, topology_list)

    water_distribution.substitutions.update({
        "L": L,
        "Sr": sources,
        "Sk": sinks,
        "\\epsilon": 0.26e-6,
        "H_{min}": h_min,
        "\\rho": 1000,
        "\\mu": 8.9e-4,
        "g": 9.81,
        "D_{max}": 1.016,
        "D_{min}": 0.3048,
        "F_{max}": 1e20,
        "C_s": 1,
    })

    water_distribution.cost = water_distribution['C']

    sol = water_distribution.localsolve(verbosity=4, reltol=1e-7, iteration_limit=150)

