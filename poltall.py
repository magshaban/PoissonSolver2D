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


    U = np.zeros(len(zt))

    for i in range(len(zt)):
        U[i] = zt[i][0]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    z =U
    surf = ax.plot_trisurf(x, y, z, cmap=cm.Spectral,linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('  %.03f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.7, aspect=9)
    plt.title('the solution $u_h(x)$')
    plt.show()   

#    x_list = x
#    y_list = y
#    z_list = zt
#    N = int(len(zt) ** .5)
#    z = zt.reshape(N, N)
#    plt.title('The solution $u$')
#    plt.imshow(z, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)), norm=LogNorm(),
#               aspect='auto')
#    plt.colorbar()
#    plt.show()
