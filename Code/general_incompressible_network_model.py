from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
import numpy as np


class IncompressibleFluidNetworkDistribution(Model):
    def setup(self, N):
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
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
        slackCost = Variable("C_s", "-", "Slack Cost")
        totalCost = Variable("C", "-", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        conCost = Variable("C_c", "-", "Connection Cost")
        rho = Variable("\\rho", "kg/m^3", "Density")
        mu = Variable("\\mu", "kg/m/s", "Viscosity")
        g = Variable("g", "m/s^2", "Gravity")
        # flowCostFactor = Variable("FC_f", "s/m^5.5")

        constraints = []

        with SignomialsEnabled():

            for i in range(0, N):
                constraints.extend([
                    Tight([sink[i] + sum(flow[i, :]) <= slack_1[i] * (source[i] + sum(flow[:, i]))]),
                    Tight([slack_1[i] >= 1]),
                    Tight([source[i] + sum(flow[:, i]) <= slack_2[i] * (sink[i] + sum(flow[i, :]))]),
                    Tight([slack_2[i] >= 1]),
                    H[i] >= H_min[i]
                ])
                if i > 0:
                    constraints.extend([sum(H * connect[:, i]) >= H[i] + sum(H_loss[:, i] * connect[:, i])])
                for j in range(0, N):
                    constraints += [flow[i, j] <= connect[i, j] * maxFlow[i, j],
                                    connect[i, j] <= 1,
                                    flowCost[i, j] == 1.1 * D[i, j] ** 1.5 * L[i, j] * units.s / units.m ** 5.5,
                                    H_loss[i, j] == f[i, j] * L[i, j] * V[i, j] ** 2 / (2 * D[i, j] * g),
                                    V[i, j] == 4 * flow[i, j] / (np.pi * D[i, j] ** 2),
                                    relRough[i, j] == rough[i, j] / D[i, j],
                                    Re[i, j] == rho * V[i, j] * D[i, j] / mu,
                                    D[i, j] <= D_max,
                                    D[i, j] >= D_min*(connect[i, j] + connect[j, i])
                                    ]

                    constraints += [f[i, j] ** 2.39794 >= 3.26853e-06 * Re[i, j] ** 0.0574443 * relRough[
                        i, j] ** 0.364794 + 0.0001773 * Re[i, j] ** -0.529499 * relRough[
                                        i, j] ** -0.0810121
                                    + 0.00301918 * Re[i, j] ** -0.0220498 * relRough[i, j] ** 1.73526 + 0.0734922 * Re[
                                        i, j] ** -1.13629 * relRough[i, j] ** 0.0574655
                                    + 0.000214297 * Re[i, j] ** 0.00035242 * relRough[i, j] ** 0.823896]
                    constraints += [f[i,j] <= 10]
            for i in range(0, N):
                for j in range(i + 1, N):
                    constraints.extend([connect[i, j] * connect[j, i] <= 1e-20,
                                        D[i, j] == D[j, i],
                                        relRough[i, j] == relRough[j, i],
                                        L[i, j] == L[j, i]])
            constraints += [totalCost >= np.sum(flow * flowCost + conCost * connect) * (
                    1 + slackCost * np.prod(slack_1) * np.prod(slack_2))]
            constraints += [H[0] == 100 * units.m]
        return constraints


if __name__ == '__main__':
    pass
