from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
from generate_data import define_length, define_topology
import numpy as np


class IncompressibleFluidNetworkDistribution(Model):
    def setup(self, N):
        H = VectorVariable(N, "H", "m", "Head")
        H_part = VectorVariable([N, N], "H_p", "m", "Partial Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        H_s = Variable("H_{s}", "m", "Head Source")
        source = VectorVariable(N, "Sr", "m^3/s", "Source")
        sink = VectorVariable(N, "Sk", "m^3/s", "Sink")
        rough = VectorVariable([N, N], "\\epsilon", "m", "Pipe Roughness")
        relRough = VectorVariable([N, N], "\\epsilon/D", "-", "Relative Pipe Roughness")
        flowCost = VectorVariable([N, N], "C_f", "s/m^3",
                                  "Pipe Flow Cost")
        L = VectorVariable([N, N], "L", "m", "Pipe Length")
        D = VectorVariable([N, N], "D", "m", "Pipe Diameter")
        maxFlow = VectorVariable([N, N], "F_{max}", "m^3/s",
                                 'Maximum Flow Rate')
        connect = VectorVariable([N, N], "x", "-", "connectivity")  # Integer Variable
        flow = VectorVariable([N, N], "F", "m^3/s", "Flow Rate")
        V = VectorVariable([N, N], "V_f", "m/s", "Flow Velocity")
        H_loss = VectorVariable([N, N], "H_L", "m", "Head Loss")
        Re = VectorVariable([N, N], "Re", "-", "Reynold's Number")
        f = VectorVariable([N, N], "f", "-", "Friction Factor")
        slack_1 = VectorVariable(N, "S_1", "-", "First Slack")
        slack_2 = VectorVariable(N, "S_2", "-", "Second Slack")
        udr = Variable("Udr", "-", "Undirected Topology")
        slackCost = Variable("C_s", "-", "Slack Cost")
        totalCost = Variable("C", "-", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        conCost = Variable("C_c", "-", "Connection Cost")
        rho = Variable("\\rho", "kg/m^3", "Density")
        mu = Variable("\\mu", "kg/m/s", "Viscosity")
        g = Variable("g", "m/s^2", "Gravity")

        constraints = []

        with SignomialsEnabled():
            for i in range(0, N):
                # conservation of mass constraints
                constraints.extend([
                    Tight([sink[i] + sum(flow[i, :]) <= slack_1[i] * (source[i] + sum(flow[:, i]))]),
                    Tight([source[i] + sum(flow[:, i]) <= slack_2[i] * (sink[i] + sum(flow[i, :]))]),
                    Tight([slack_2[i] >= 1]), Tight([slack_1[i] >= 1]),
                ])
                # energy constraints
                constraints.extend([
                    H[i] >= sum(H_loss[i, :] * connect[i, :]) + sum(H_part[i, :] * connect[i, :]),
                    H[i] <= sum(H_part[i, :]),
                    H[i] >= H_min[i],
                ])

                for j in range(0, N):
                    # topology constraints
                    constraints += [
                        connect[i, j] <= 1,
                        connect[i, j] * connect[j, i] <= udr,
                        flow[i, j] <= connect[i, j] * maxFlow[i, j],
                        H_part[i, j] <= connect[i, j]*H_s,
                        D[i, j] <= D_max,
                        D[i, j] >= D_min*(connect[i, j] + connect[j, i]),
                    ]
                    # flow cost
                    constraints += [
                        flowCost[i, j] == 1.1 * D[i, j] ** 1.5 * L[i, j] * units.s / units.m ** 5.5,
                    ]
                    # Darcy-Weisbach Equations
                    constraints += [
                        H_loss[i, j] == f[i, j] * L[i, j] * V[i, j] ** 2 / (2 * D[i, j] * g),
                        V[i, j] == 4 * flow[i, j] / (np.pi * D[i, j] ** 2),
                        relRough[i, j] == rough[i, j] / D[i, j],
                        Re[i, j] == rho * V[i, j] * D[i, j] / mu,
                    ]
                    # friction factor posynomial approximation
                    constraints += [f[i, j] ** 2.39794 >= 3.26853e-06 * Re[i, j] ** 0.0574443 * relRough[
                        i, j] ** 0.364794 + 0.0001773 * Re[i, j] ** -0.529499 * relRough[
                                        i, j] ** -0.0810121
                                    + 0.00301918 * Re[i, j] ** -0.0220498 * relRough[i, j] ** 1.73526 + 0.0734922 * Re[
                                        i, j] ** -1.13629 * relRough[i, j] ** 0.0574655
                                    + 0.000214297 * Re[i, j] ** 0.00035242 * relRough[i, j] ** 0.823896]
                    constraints += [f[i, j] <= 10]

                for j in range(i + 1, N):
                    constraints.extend([
                                        D[i, j] == D[j, i],
                                        relRough[i, j] == relRough[j, i],
                                        L[i, j] == L[j, i]])
            constraints += [totalCost >= np.sum(flow * flowCost + conCost * connect) * (1 + slackCost * np.prod(slack_1) * np.prod(slack_2))]
            constraints += [H[0] == H_s]
        return constraints


if __name__ == '__main__':
    # Somewhat large problem
    """
    N = 32
    sinks = [0, 890 / 3600.0, 850 / 3600.0, 130 / 3600.0, 725 / 3600.0, 1005 / 3600.0, 1350 / 3600.0, 550 / 3600.0,
             525 / 3600.0,
             525 / 3600.0, 500 / 3600.0, 560 / 3600.0, 940 / 3600.0, 615 / 3600.0, 280 / 3600.0, 310 / 3600.0,
             865 / 3600.0, 1345 / 3600.0, 60 / 3600.0, 1275 / 3600.0, 930 / 3600.0, 485 / 3600.0, 1045 / 3600.0,
             820 / 3600.0, 170 / 3600.0,
             900 / 3600.0, 370 / 3600.0, 290 / 3600.0, 360 / 3600.0, 360 / 3600.0, 105 / 3600.0, 805 / 3600.0]
    sinks = [i / 10.0 for i in sinks]
    sources = [0 for i in xrange(N)]
    sources[0] = sum(sinks)
    topology_list = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11],
                     [11, 12],
                     [9, 13], [13, 14], [14, 15], [16, 15], [17, 16], [18, 17], [2, 18], [2, 19], [19, 20],
                     [20, 21],
                     [19, 22], [22, 23], [23, 24], [25, 24], [26, 25], [15, 26], [22, 27], [27, 28], [28, 29],
                     [29, 30],
                     [31, 30], [24, 31]]
    L = [100, 1350, 900, 1150, 1450, 450, 850, 850, 800, 950, 1200, 3500, 800, 500, 550, 2730, 1750, 800, 400, 2200,
         1500, 500, 2650, 1230, 1300, 850, 300, 750, 1500, 2000, 1600, 150, 860, 950]
    coordinates = {1: (5021.20, 1582.17), 2: (5025.20, 2585.42), 3: (5874.22, 2588.30), 4: (6873.11, 2588.30),
                   5: (8103.51, 2585.42), 6: (8103.51, 3234.67), 7: (8106.66, 4179.28), 8: (8106.66, 5133.78),
                   9: (7318.64, 5133.78), 10: (7319.94, 5831.65), 11: (7319.94, 6671.19), 12: (5636.76, 6676.24),
                   13: (6530.63, 5133.78), 14: (5676.02, 5133.78), 15: (5021.20, 5133.78), 16: (5021.20, 4412.36),
                   17: (5021.20, 3868.52), 18: (5021.20, 3191.49), 19: (3587.87, 2588.30), 20: (3587.87, 1300.84),
                   21: (3587.87, 901.29), 22: (1978.55, 2588.30), 23: (1975.58, 4084.35), 24: (1980.46, 5137.63),
                   25: (3077.46, 5137.63), 26: (3933.52, 5133.78), 27: (846.04, 2588.20), 28: (-552.41, 2588.20),
                   29: (-552.38, 4369.06), 30: (-549.36, 5137.63), 31: (536.45, 5137.63), 0: (5360.71, 1354.05)}

    # Small problem
    """
    """
    N = 5
    sinks = [0, 0.89, 0.85, 0.130, 0.725]
    sources = [sum(sinks), 0, 0, 0, 0]
    topology_list = [[0, 1], [1, 2], [1, 3], [2, 4] ,[3,4]]
    coordinates = {0: (0, 1000), 1: (0, 0), 2: (1000, 0), 3: (-1000, 0), 4: (1000, -1000)}
    """
    """
    # Small tree problem
    N = 3
    sinks = [0, 0.8, 0.65]
    sources = [sum(sinks), 0, 0]
    topology_list = [[0, 1], [1, 2]]
    coordinates = {0: (0, 1000), 1: (0, 0), 2: (-1000, 0)}
    """
    N = 2
    sinks = [0, 0.8]
    sources = [sum(sinks), 0]
    topology_list = [[0, 1]]
    coordinates = {0: (0, 1000), 1: (0, 0)}
    connect = define_topology(topology_list, N)

    L = [[0 for i in xrange(N)] for j in xrange(N)]
    for p1 in coordinates:
        for p2 in coordinates:
            L[p1][p2] = np.sqrt((coordinates[p1][0] - coordinates[p2][0]) ** 2 +
                                (coordinates[p1][1] - coordinates[p2][1]) ** 2)

    water_distribution = IncompressibleFluidNetworkDistribution(N)

    water_distribution.substitutions.update({
        "L": L,
        "Sr": sources,
        "Sk": sinks,
        "\\epsilon": 0.26e-6,
        "H_{min}": 30,
        "H_{s}": 2000,
        "\\rho": 1000,
        "\\mu": 8.9e-4,
        "g": 9.81,
        "D_{max}": 1.016,
        "D_{min}": 0.3048,
        "F_{max}": 1e20,
        "C_s": 1,
        "Udr": 0,
        "C_c": 1,
        "x": connect,
    })

    water_distribution.cost = water_distribution['C']
    sol = water_distribution.localsolve(verbosity=2, reltol=1e-2, iteration_limit=1500)

