
from builtins import range
from gpkit import Model, Variable, VectorVariable, SignomialsEnabled, units
from gpkit.constraints.tight import Tight
import numpy as np

# Known Topology Fluid Network Design
class KT_FND(Model):
    # Known Topology Fluid Network Design
    def setup(self, N, topology_dict, friction='DW', penalty=10.):
        n_p = len(topology_dict) # number of pipes
        H = VectorVariable(N, "H", "m", "Head")
        H_min = VectorVariable(N, "H_{min}", "m", "Minimal Head Required")
        source = VectorVariable(N, "\dot{V}_+", "m^3/s", "Source", pr=20)
        sink = VectorVariable(N, "\dot{V}_-", "m^3/s", "Sink", pr=20)
        rough = VectorVariable(n_p, "\\epsilon", "m", "Pipe Roughness")
        pipeCost = VectorVariable(n_p, "c_p", "-",
                                  "Pipe Cost")
        L = VectorVariable(n_p, "L", "m", "Pipe Length")
        D = VectorVariable(n_p, "D", "m", "Pipe Diameter")
        flow = VectorVariable(n_p, "q", "m^3/s", "Flow Rate")
        V = VectorVariable(n_p, "v_f", "m/s", "Flow Velocity")
        maxV = VectorVariable(n_p, "v_{max}", 1e20*np.ones(n_p), "m/s", 'Maximum Flow Velocity')
        H_loss = VectorVariable(n_p, "H_L", "m", "Head Loss")
        slack_out = VectorVariable(N, "S_{out}", "-", "Outflow Slack")
        slack_in = VectorVariable(N, "S_{in}", "-", "Inflow Slack")
        slack_h = VectorVariable(n_p, "S_h", "-", "Head Slack")
        totalCost = Variable("C", "-", "Total Cost")
        D_max = Variable("D_{max}", "m", "Maximum Diameter")
        D_min = Variable("D_{min}", "m", "Minimum Diameter")
        rho = Variable("\\rho", 1000, "kg/m^3", "Density")
        mu = Variable("\\mu", 8.9e-4, "kg/m/s", "Viscosity")
        g = Variable("g", 9.81, "m/s^2", "Gravity")

        if friction == 'DW':
            relRough = VectorVariable(n_p, "\\bar{\\epsilon}", "-", "Relative Pipe Roughness")
            Re = VectorVariable(n_p, "Re", "-", "Reynold's Number")
            f = VectorVariable(n_p, "f", "-", "Friction Factor")

        constraints = []

        with SignomialsEnabled():

            for i in range(0, N):
                flow_in = sink[i]
                flow_out = source[i]
                for pipe_index, node in list(topology_dict.items()):
                    if node[0] == i:
                        flow_in += flow[pipe_index]
                    if node[1] == i:
                        flow_out += flow[pipe_index]

                constraints.extend([
                    Tight([flow_in <= slack_out[i] * flow_out]),
                    Tight([flow_out <= slack_in[i] * flow_in]),
                    Tight([slack_in[i] >= 1]), Tight([slack_out[i] >= 1]),
                    H[i] >= H_min[i]
                ])
                # Head loss constraints
                for pipe_index, node in list(topology_dict.items()):
                    if node[0] == i:
                        constraints.extend([
                            Tight([H[node[0]] >= H_loss[pipe_index] + H[node[1]]]),
                            Tight([H[node[0]] <= slack_h[pipe_index]*(H_loss[pipe_index] + H[node[1]])]),
                            Tight([slack_h[pipe_index] >= 1]),
                        ])

            for pipe_index in range(n_p):
                constraints += [V[pipe_index] <= maxV,
                                pipeCost[pipe_index] == 1.1 * D[pipe_index] ** 1.5 * L[pipe_index] / units.m ** 2.5,
                                V[pipe_index] == 4 * flow[pipe_index] / (np.pi * D[pipe_index] ** 2),
                                D[pipe_index] <= D_max,
                                D[pipe_index] >= D_min]

                if friction == "HW":
                    constraints += [H_loss[pipe_index] == 10.67*L[pipe_index] * (flow[pipe_index]/units('m^3/s')) ** 1.852 /
                                                        ((rough[pipe_index]/units('m'))**1.852*(D[pipe_index]/units('m'))**4.8704)]
                                # S (hydraulic slope H/L) = 10.67*Q^1.852 (volumetric flow rate) /
                                # C^1.852 (pipe roughness) /d^4.8704 (pipe diameter)

                if friction == 'DW':
                    constraints += [H_loss[pipe_index] == f[pipe_index] * L[pipe_index] * V[pipe_index] ** 2 / (
                                            2 * D[pipe_index] * g),
                                relRough[pipe_index] == rough[pipe_index] / D[pipe_index],
                                Re[pipe_index] == rho * V[pipe_index] * D[pipe_index] / mu,
                    # From frictionFactorFitting.py
                                f[pipe_index] ** 2.39794 >= 3.26853e-06 * Re[pipe_index] ** 0.0574443 *
                                relRough[pipe_index] ** 0.364794 + 0.0001773 * Re[pipe_index] ** -0.529499 *
                                relRough[pipe_index] ** -0.0810121 + 0.00301918 * Re[pipe_index] ** -0.0220498 *
                                relRough[pipe_index] ** 1.73526 + 0.0734922 * Re[pipe_index] ** -1.13629 *
                                relRough[pipe_index] ** 0.0574655 + 0.000214297 * Re[pipe_index] ** 0.00035242 *
                                relRough[pipe_index] ** 0.823896,
                                    f[pipe_index] <= 1]

            constraints += [totalCost >= np.sum(pipeCost) *
                            (np.prod(slack_out) * np.prod(slack_in) * np.prod(slack_h)**penalty)]
        return constraints

def subs_with_dict(m, varkey, dict):
    for idx, val in list(dict.items()):
        m.substitutions.update({varkey[idx]:val})

