# The part of <possionSolver2D.py> 
# This script to plot all the required graphs.
# also this script will polt the delaunay triangulation.

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.spatial import Delaunay
from matplotlib.colors import LogNorm


def delunay_plot(points):
    tri = Delaunay(points)
    plt.triplot(points[:, 0], points[:, 1], tri.simplices.copy())
    plt.plot(points[:, 0], points[:, 1], 'o')
    plt.title('Delanauy Trianglaution of the domain $\Omega = [-1,1] ^2$')
    plt.show()


def poltall(x, y, zt, mat):
    # to plot the mass matrix 
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('The Stiffness Matrix')
    ax1.matshow(mat)
    ax2.spy(mat)
    plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # f = func(p[0,:], p[1,:])
    # Plot the surface.
    surf = ax.plot_surface(x, y, zt, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.set_zlim(-1.33, 4)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.7, aspect=15)
    plt.title('The solution $u$')
    plt.show()

    x_list = x
    y_list = y
    z_list = zt
    N = int(len(zt) ** .5)
    z = zt.reshape(N, N)
    plt.title('The solution $u$')
    plt.imshow(z, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)), norm=LogNorm(),
               aspect='auto')
    plt.colorbar()
    plt.show()
