import numpy
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.special
from scipy.special import sph_harm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import ListedColormap

import warnings
# Ignore all warnings
warnings.filterwarnings("ignore")

def hydrogen_wf(n, l, m, X, Y, Z):
    R = numpy.sqrt(X**2 + Y**2 + Z**2)
    Theta = numpy.arccos(Z/R)
    Phi = numpy.arctan2(Y, X)

    rho = 2.0 * R / n
    s_harm = sph_harm(m, l, Phi, Theta)
    l_poly = scipy.special.genlaguerre(n - l - 1, 2 * l + 1)(rho)

    prefactor = numpy.sqrt((2.0 / n)**3 * math.factorial(n - l - 1) / (2.0 * n * math.factorial(n + l)))
    wf = prefactor * numpy.exp(-rho / 2.0) * rho**l * s_harm * l_poly
    wf = numpy.nan_to_num(wf)
    return wf

def update_plot(val):
    n = int(slider_n.val)
    l = int(slider_l.val)
    m = int(slider_m.val)

    #print(n_init, l, m)

    data = hydrogen_wf(n, l, m, X, Y, Z)
    data = abs(data)**2

    index = int((slider_y.val - zmin) / dz)
    im.set_data(data[index, :, :])
    #im.set_clim(vmin=0, vmax=numpy.max(data))
    ax.set_title("Hydrogen Orbital xz Slice (y=" + str("%.2f" % slider_y.val) +
                 "): n=" + str(n) + ", l=" + str(l) + ", m=" + str(m))


def plot_Horbs(n=4):
    
    global fig, ax, im, slider_n, slider_l, slider_m, slider_y, X, Y, Z, zmin, dz

    dz = 0.5
    zmin = -15
    zmax = 15
    x = numpy.arange(zmin, zmax, dz)
    y = numpy.arange(zmin, zmax, dz)
    z = numpy.arange(zmin, zmax, dz)
    X, Y, Z = numpy.meshgrid(x, y, z)

    n_init = n
    l_init = n-1
    m_init = 0

    data = hydrogen_wf(n_init, l_init, m_init, X, Y, Z)
    data = abs(data)**2

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.15, bottom=0.25)

    im = plt.imshow(data[int((0 - zmin) / dz), :, :], 
            vmin=0, vmax=numpy.max(data), 
            extent=[zmin, zmax, zmin, zmax], 
            cmap='Purples')
    plt.colorbar()

    slider_n = Slider(plt.axes([0.25, 0.1, 0.65, 0.03]), "n", n_init, n_init, valinit=n_init, valstep=1)
    slider_l = Slider(plt.axes([0.25, 0.05, 0.65, 0.03]), "l", 0, n_init-1, valinit=l_init, valstep=1)
    slider_m = Slider(plt.axes([0.25, 0.0, 0.65, 0.03]), "m", -(n_init-1), (n_init-1), valinit=m_init, valstep=1)
    slider_y = Slider(plt.axes([0.25, 0.15, 0.65, 0.03]), "Y", z[0], z[len(z) - 1], valinit=0)

    ax.set_title("Hydrogen Orbital xz Slice (y=" + str("%.2f" % slider_y.val) +
             "): n=" + str(n_init) + ", l=" + str(l_init) + ", m=" + str(m_init))


    slider_n.on_changed(update_plot)
    slider_l.on_changed(update_plot)
    slider_m.on_changed(update_plot)
    slider_y.on_changed(update_plot)

    plt.show()

