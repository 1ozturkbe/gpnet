from gpfit.fit import fit
import numpy as np
from numpy import logspace, log, meshgrid
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt


# fixed initial guess for fitting


def generate_friction_factor_data(Re, relRough):
    Re, relRough = meshgrid(Re, relRough)
    Re = Re.flatten()
    relRough = relRough.flatten()

    f = 0.25 / log(relRough / 3.7 + 5.74 / Re ** 0.9) ** 2

    x = [Re, relRough]
    x = np.log(x)
    y = np.log(f)
    return x, y


def fit_friction_factor_data(x, y, k, method):
    c, error = fit(x, y, k, method)
    return c, error


def plot_friction_factor_fit_k_equal_5(Re, relRough):
    u_1, u_2 = meshgrid(Re, relRough)

    fig1 = plt.figure()
    ax = Axes3D(fig1)
    w = 0.25 / log(u_2 / 3.7 + 5.74 / u_1 ** 0.9) ** 2
    ax.plot_surface(np.log(u_1), np.log(u_2), np.log(w))

    plt.figure()
    plt.contour(np.log(u_1), np.log(u_2), np.log(w))

    fig2 = plt.figure()
    ax2 = Axes3D(fig2)
    w_app = (3.26853e-06 * u_1 ** 0.0574443 * u_2 ** 0.364794 + 0.0001773 * u_1 ** -0.529499 * u_2 ** -0.0810121
             + 0.00301918 * u_1 ** -0.0220498 * u_2 ** 1.73526 + 0.0734922 * u_1 ** -1.13629 * u_2 ** 0.0574655
             + 0.000214297 * u_1 ** 0.00035242 * u_2 ** 0.823896) ** (1 / 2.39794)
    ax2.plot_surface(np.log(u_1), np.log(u_2), np.log(w_app))

    plt.figure()
    plt.contour(np.log(u_1), np.log(u_2), np.log(w_app))


if __name__ == '__main__':
    Re = logspace(3, 8, 50)
    relRough = logspace(-5, -1, 50)
    xData, yData = generate_friction_factor_data(Re, relRough)
    c, e = fit_friction_factor_data(xData, yData, 5, 'SMA')
    plot_friction_factor_fit_k_equal_5(Re, relRough)
