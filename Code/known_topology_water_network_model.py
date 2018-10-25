from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
import numpy as np


# Hazen-Williams Known Topology Water Network Design
class HWKTWND(Model):
    def setup(self, N, topology_list):
        number_of_pipes = len(topology_list)
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        source = VectorVariable(N, "Sr", "m^3/s", "Source")
        sink = VectorVariable(N, "Sk", "m^3/s", "Sink")
        pipeCost = VectorVariable(number_of_pipes, "P_f", "-", "Pipe Cost")
        L = VectorVariable(number_of_pipes, "L", "m", "Pipe Length")
        D = VectorVariable(number_of_pipes, "D", "m", "Pipe Diameter")
        maxFlow = Variable("F_{max}", "m^3/s", 'Maximum Flow Rate')
        flow = VectorVariable(number_of_pipes, "F", "m^3/s", "Flow Rate")
        H_loss = VectorVariable(number_of_pipes, "H_L", "m", "Head Loss")
        slack_1 = VectorVariable(N, "S_1", "-", "First Slack")
        slack_2 = VectorVariable(N, "S_2", "-", "Second Slack")
        slackCost = Variable("C_s", "-", "Slack Cost")
        totalCost = Variable("C", "m^3/s", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        rough = Variable("\\epsilon", "-", "Pipe Roughness")

        constraints = []

        source_dict = {}

        with SignomialsEnabled():

            for i in range(0, N):
                first_flow_p = sink[i]
                second_flow_p = source[i]
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        first_flow_p += flow[pipe_index]
                    if pipe[1] == i:
                        second_flow_p += flow[pipe_index]

                constraints.extend([
                    Tight([first_flow_p <= slack_1[i] * second_flow_p]),
                    Tight([slack_1[i] >= 1]),
                    # Tight([slack_1[i] <= 10]),
                    Tight([second_flow_p <= slack_2[i] * first_flow_p]),
                    Tight([slack_2[i] >= 1]),
                    # Tight([slack_2[i] <= 10]),
                    H[i] >= H_min[i]
                ])
                left_flow_p = H[i]
                right_flow_p = 0
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        H_part = Variable("H_{%s,%s}" % (i, pipe[1]), "m")
                        if pipe[1] in source_dict:
                            source_dict[pipe[1]].append(H_part)
                        else:
                            source_dict[pipe[1]] = [H_part]
                        right_flow_p += H_part
                        right_flow_p += H_loss[pipe_index]
                if right_flow_p != 0:
                    constraints.extend(Tight([left_flow_p >= right_flow_p]))

            for key in source_dict.keys():
                constraints.extend(Tight([H[key] <= sum(source_dict[key])]))

            for pipe_index in xrange(number_of_pipes):
                constraints += [flow[pipe_index] <= maxFlow,
                                pipeCost[pipe_index] == 1.1 * D[pipe_index] ** 1.5 * L[pipe_index] / units.m ** 2.5,
                                H_loss[pipe_index] == L[pipe_index] * (flow[pipe_index] / rough) ** 1.8099 / (
                                            994.62 * (D[pipe_index]) ** 4.8099),
                                D[pipe_index] <= D_max,
                                D[pipe_index] >= D_min]

            constraints += [totalCost >= np.sum(flow * pipeCost) * (slackCost * np.prod(slack_1) * np.prod(slack_2))]
            constraints += [H[0] == 100 * units.m]
        return constraints


if __name__ == '__main__':
    N = 32
    sinks = [0, 890 / 3600.0, 850 / 3600.0, 130 / 3600.0, 725 / 3600.0, 1005 / 3600.0, 1350 / 3600.0, 550 / 3600.0,
             525 / 3600.0,
             525 / 3600.0, 500 / 3600.0, 560 / 3600.0, 940 / 3600.0, 615 / 3600.0, 280 / 3600.0, 310 / 3600.0,
             865 / 3600.0, 1345 / 3600.0, 60 / 3600.0, 1275 / 3600.0, 930 / 3600.0, 485 / 3600.0, 1045 / 3600.0,
             820 / 3600.0, 170 / 3600.0,
             900 / 3600.0, 370 / 3600.0, 290 / 3600.0, 360 / 3600.0, 360 / 3600.0, 105 / 3600.0, 805 / 3600.0]
    sources = [0 for i in xrange(N)]
    sources[0] = sum(sinks)
    topology_list = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11],
                     [11, 12],
                     [9, 13], [13, 14], [14, 15], [16, 15], [17, 16], [18, 17], [2, 18], [2, 19], [19, 20], [20, 21],
                     [19, 22], [22, 23], [23, 24], [25, 24], [26, 25], [15, 26], [22, 27], [27, 28], [28, 29], [29, 30],
                     [31, 30], [24, 31]]
    L = [100, 1350, 900, 1150, 1450, 450, 850, 850, 800, 950, 1200, 3500, 800, 500, 550, 2730, 1750, 800, 400, 2200,
         1500, 500, 2650, 1230, 1300, 850, 300, 750, 1500, 2000, 1600, 150, 860, 950]

    water_distribution = HWKTWND(N, topology_list)

    water_distribution.substitutions.update({
        "L": L,
        "Sr": sources,
        "Sk": sinks,
        "\\epsilon": 130,
        "H_{min}": 10,
        "D_{max}": 1.016,
        "D_{min}": 0.3048,
        "F_{max}": 1e20,
        "C_s": 1,
    })

    water_distribution.cost = water_distribution['C']
    sol = water_distribution.localsolve(verbosity=4, reltol=1e-4, iteration_limit=1500)
