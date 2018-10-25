from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
import numpy as np


# Known Topology Fluid Network Design
class KTFND(Model):
    def setup(self, N, topology_list):
        number_of_pipes = len(topology_list)
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        H_s = Variable("H_{s}", "m", "Head Source")
        source = VectorVariable(N, "Sr", "m^3/s", "Source")
        sink = VectorVariable(N, "Sk", "m^3/s", "Sink")
        rough = Variable("\\epsilon", "m", "Pipe Roughness")
        relRough = VectorVariable(number_of_pipes, "\\epsilon/D", "-", "Relative Pipe Roughness")
        pipeCost = VectorVariable(number_of_pipes, "P_f", "-",
                                  "Pipe Cost")
        L = VectorVariable(number_of_pipes, "L", "m", "Pipe Length")
        D = VectorVariable(number_of_pipes, "D", "m", "Pipe Diameter")
        maxFlow = Variable("F_{max}", "m^3/s", 'Maximum Flow Rate')
        flow = VectorVariable(number_of_pipes, "F", "m^3/s", "Flow Rate")
        V = VectorVariable(number_of_pipes, "V_f", "m/s", "Flow Velocity")
        H_loss = VectorVariable(number_of_pipes, "H_L", "m", "Head Loss")
        Re = VectorVariable(number_of_pipes, "Re", "-", "Reynold's Number")
        f = VectorVariable(number_of_pipes, "f", "-", "Friction Factor")
        slack_1 = VectorVariable(N, "S_1", "-", "First Slack")
        slack_2 = VectorVariable(N, "S_2", "-", "Second Slack")
        slackCost = Variable("C_s", "-", "Slack Cost")
        totalCost = Variable("C", "m^3/s", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        rho = Variable("\\rho", "kg/m^3", "Density")
        mu = Variable("\\mu", "kg/m/s", "Viscosity")
        g = Variable("g", "m/s^2", "Gravity")

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
                    Tight([second_flow_p <= slack_2[i] * first_flow_p]),
                    Tight([slack_2[i] >= 1]), Tight([slack_1[i] >= 1]),
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
                                H_loss[pipe_index] == f[pipe_index] * L[pipe_index] * V[pipe_index] ** 2 / (2 * D[pipe_index] * g),
                                V[pipe_index] == 4 * flow[pipe_index] / (np.pi * D[pipe_index] ** 2),
                                relRough[pipe_index] == rough / D[pipe_index],
                                Re[pipe_index] == rho * V[pipe_index] * D[pipe_index] / mu,
                                D[pipe_index] <= D_max,
                                D[pipe_index] >= D_min]

                constraints += [f[pipe_index] ** 2.39794 >= 3.26853e-06 * Re[pipe_index] ** 0.0574443 *
                                relRough[pipe_index] ** 0.364794 + 0.0001773 * Re[pipe_index] ** -0.529499 *
                                relRough[pipe_index] ** -0.0810121 + 0.00301918 * Re[pipe_index] ** -0.0220498 *
                                relRough[pipe_index] ** 1.73526 + 0.0734922 * Re[pipe_index] ** -1.13629 *
                                relRough[pipe_index] ** 0.0574655 + 0.000214297 * Re[pipe_index] ** 0.00035242 *
                                relRough[pipe_index] ** 0.823896]

            constraints += [totalCost >= np.sum(flow * pipeCost) * (slackCost * np.prod(slack_1) * np.prod(slack_2))]
            constraints += [H[0] == H_s]
        return constraints


if __name__ == '__main__':
    pass
