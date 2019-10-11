from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
import numpy as np

# Known Topology Fluid Network Design
class DWKTFND(Model):
    # Darcy-Weisbach Known Topology Water Network Design
    def setup(self, N, topology_list, penalty=10.):
        number_of_pipes = len(topology_list)
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        source = VectorVariable(N, "\dot{V}_+", "m^3/s", "Source")
        sink = VectorVariable(N, "\dot{V}_-", "m^3/s", "Sink")
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
        slack_out = VectorVariable(N, "S_{out}", "-", "Outflow Slack")
        slack_in = VectorVariable(N, "S_{in}", "-", "Inflow Slack")
        slack_h = VectorVariable(number_of_pipes, "S_h", "-", "Head Slack")
        totalCost = Variable("C", "m^3/s", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        rho = Variable("\\rho", "kg/m^3", "Density")
        mu = Variable("\\mu", "kg/m/s", "Viscosity")
        g = Variable("g", "m/s^2", "Gravity")

        constraints = []

        with SignomialsEnabled():

            for i in range(0, N):
                flow_in = sink[i]
                flow_out = source[i]
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        flow_in += flow[pipe_index]
                    if pipe[1] == i:
                        flow_out += flow[pipe_index]

                constraints.extend([
                    Tight([flow_in <= slack_out[i] * flow_out]),
                    Tight([flow_out <= slack_in[i] * flow_in]),
                    Tight([slack_in[i] >= 1]), Tight([slack_out[i] >= 1]),
                    H[i] >= H_min[i]
                ])
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        constraints.extend([
                            Tight([H[i] >= H_loss[pipe_index] + H[pipe[1]]]),
                            Tight([H[i] <= slack_h[pipe_index]*(H_loss[pipe_index] + H[pipe[1]])]),
                            Tight([slack_h[pipe_index] >= 1]),
                        ])
            for pipe_index in xrange(number_of_pipes):
                constraints += [flow[pipe_index] <= maxFlow,
                                pipeCost[pipe_index] == 1.1 * D[pipe_index] ** 1.5 * L[pipe_index] / units.m ** 2.5,
                                H_loss[pipe_index] == f[pipe_index] * L[pipe_index] * V[pipe_index] ** 2 / (
                                            2 * D[pipe_index] * g),
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
                constraints += [f[pipe_index] <= 1]

            constraints += [totalCost >= np.sum(flow * pipeCost) *
                            (np.prod(slack_out) * np.prod(slack_in) * np.prod(slack_h)**penalty)]
        return constraints

class HWKTWND(Model):
    # Hazen-Williams Known Topology Water Network Design
    def setup(self, N, topology_list, penalty=10.):
        number_of_pipes = len(topology_list)
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        source = VectorVariable(N, "\dot{V}_+", "m^3/s", "Source")
        sink = VectorVariable(N, "\dot{V}_-", "m^3/s", "Sink")
        pipeCost = VectorVariable(number_of_pipes, "P_f", "-", "Pipe Cost")
        L = VectorVariable(number_of_pipes, "L", "m", "Pipe Length")
        D = VectorVariable(number_of_pipes, "D", "m", "Pipe Diameter")
        maxFlow = Variable("F_{max}", "m^3/s", 'Maximum Flow Rate')
        flow = VectorVariable(number_of_pipes, "F", "m^3/s", "Flow Rate")
        H_loss = VectorVariable(number_of_pipes, "H_L", "m", "Head Loss")
        slack_out = VectorVariable(N, "S_{out}", "-", "Outflow Slack")
        slack_in = VectorVariable(N, "S_{in}", "-", "Inflow Slack")
        slack_h = VectorVariable(number_of_pipes, "S_h", "-", "Head Slack")
        slackCost = Variable("C_s", "-", "Slack Cost")
        totalCost = Variable("C", "m^3/s", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        rough = Variable("\\epsilon", "-", "Pipe Roughness")

        constraints = []

        with SignomialsEnabled():

            for i in range(0, N):
                flow_in = sink[i]
                flow_out = source[i]
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        flow_in += flow[pipe_index]
                    if pipe[1] == i:
                        flow_out += flow[pipe_index]

                constraints.extend([
                    Tight([flow_in <= slack_out[i] * flow_out]),
                    Tight([slack_out[i] >= 1]),
                    # Tight([slack_out[i] <= 10]),
                    Tight([flow_out <= slack_in[i] * flow_in]),
                    Tight([slack_in[i] >= 1]),
                    # Tight([slack_in[i] <= 10]),
                    H[i] >= H_min[i]
                ])
                for pipe_index, pipe in enumerate(topology_list):
                    if pipe[0] == i:
                        constraints.extend([
                            Tight([H[i] >= H_loss[pipe_index] + H[pipe[1]]]),
                            Tight([H[i] <= slack_h[pipe_index]*(H_loss[pipe_index] + H[pipe[1]])]),
                            Tight([slack_h[pipe_index] >= 1]),
                        ])

            for pipe_index in xrange(number_of_pipes):
                constraints += [flow[pipe_index] <= maxFlow,
                                pipeCost[pipe_index] == 1.1 * D[pipe_index] ** 1.5 * L[pipe_index] / units.m ** 2.5,
                                # H_loss[pipe_index] == L[pipe_index] * (flow[pipe_index] / rough) ** 1.8099 / (
                                #             994.62 * (D[pipe_index]) ** 4.8099),
                                H_loss[pipe_index] == 10.67*L[pipe_index] * (flow[pipe_index]/units('m^3/s')) ** 1.852 /
                                                        (rough**1.852*(D[pipe_index]/units('m'))**4.8704),
                                # S (hydraulic slope H/L) = 10.67*Q^1.852 (volumetric flow rate) /
                                # C^1.852 (pipe roughness) /d^4.8704 (pipe diameter)
                                D[pipe_index] <= D_max,
                                D[pipe_index] >= D_min]

            constraints += [totalCost >= np.sum(flow * pipeCost) *
                                        (np.prod(slack_out) * np.prod(slack_in) * np.prod(slack_h)**penalty)]
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

    m = DWKTWND(N, topology_list)

    m.substitutions.update({
        "L": L,
        "\dot{V}_+": sources,
        "\dot{V}_-": sinks,
        "\\epsilon": 130,
        "H_{min}": 10,
        "D_{max}": 1.016,
        "D_{min}": 0.3048,
        "F_{max}": 1e20,
        "C_s": 1,
    })

    m.cost = m['C']
    sol = m.localsolve(verbosity=4, reltol=1e-4, iteration_limit=1500)
