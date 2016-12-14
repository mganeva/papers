import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np
import bornagain as ba

# experimental data
datafile = "sc11929_27.dat"

# integration range
x_range = (5, 230)
y_range = (115, 135)

# noise accounts for the constant background noise
noise = 0.2

# settings for plotting
settings_hexagonal = {'title': "Hexagonal Lattice",
                      'filename': "results/sim2016-12-14T13:09:28.dat",
                      'f': 1.0}
settings_fcc = {'title': r"FCC Lattice with $\beta=64.55^{\circ}$",
                'filename': "results/sim2016-12-14T13:27:46.dat",
                'f': 1.0}


def load_real_data(filename=datafile, rot90=0):
    """
    Fill histogram representing our detector with intensity data from tif file.
    """
    data = np.loadtxt(filename)
    hist = ba.Histogram2D(np.rot90(data, rot90))
    return hist


def sum_over_xrange(data, pixel_range=(0, 30)):
    arr = data.getArray()
    result = data.projectionY(120)
    n = result.getTotalNumberOfBins()
    for i in range(n):
        result.setBinContent(i, np.mean(arr[n-i-1, pixel_range[0]:pixel_range[1]]))
    return result


def sum_over_yrange(data, pixel_range=(0, 30)):
    arr = data.getArray()
    result = data.projectionX(120)
    n = result.getTotalNumberOfBins()
    for i in range(n):
        result.setBinContent(i, np.mean(arr[pixel_range[0]:pixel_range[1], n-i-1]))
    return result


def plot_data(settings):
    """
    Load data and plot results
    """
    expdata = load_real_data(rot90=1)       # experimental data
    simdata = load_real_data(filename=settings['filename'])     # simulated data

    # showing the result
    # =======
    # experimental data
    # =======
    plt.suptitle(settings['title'], fontsize=18)
    ax = plt.subplot(2, 2, 1)
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    im = plt.imshow(expdata.getArray() + 1.0,
                    norm=matplotlib.colors.LogNorm(1.0, 501.0),
                    extent=[expdata.getXmin(), expdata.getXmax(), expdata.getYmin(), expdata.getYmax()],
                    aspect='auto')
    # show the integration area over x range
    ax.add_patch(
        patches.Rectangle(
            (x_range[0], 2),  # (x,y)
            x_range[1] - x_range[0],  # width
            237,  # height
            fill=False, edgecolor="red", linewidth=2
        )
    )
    # show the integration area over y range
    ax.add_patch(
        patches.Rectangle(
            (2, y_range[0]),  # (x,y)
            230,  # width
            y_range[1] - y_range[0],  # height
            fill=False, edgecolor="green", linewidth=2
        )
    )
    cb = plt.colorbar(im)
    cb.set_label(r'Intensity (arb. u.)', size=16)
    plt.xlabel(r'channel$_y$', fontsize=16)
    plt.ylabel(r'channel$_z$', fontsize=16)
    plt.title("Experiment")

    # =======
    # simulation result
    # =======
    plt.subplot(2, 2, 2)
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    im = plt.imshow(settings['f']*simdata.getArray() + noise,
                    norm=matplotlib.colors.LogNorm(1.0, 501.0),
                    extent=[simdata.getXmin(), simdata.getXmax(), simdata.getYmin(), simdata.getYmax()],
                    aspect='auto')
    cb = plt.colorbar(im)
    cb.set_label(r'Intensity (arb. u.)', size=16)
    plt.xlabel(r'channel$_y$', fontsize=16)
    plt.ylabel(r'channel$_z$', fontsize=16)
    plt.title("Simulation")

    # =======
    # slice along qy
    # =======
    plt.subplot(2, 2, 3)
    expslice = sum_over_yrange(expdata, y_range)
    simslice = sum_over_yrange(simdata, y_range)
    plt.semilogy(simslice.getBinCenters(), settings['f']*simslice.getBinValues() + noise, color='g', linewidth=2)
    plt.semilogy(expslice.getBinCenters(), expslice.getBinValues(), color='k', marker='.', linestyle='None')
    plt.xlim(expslice.getXmin(), expslice.getXmax())
    plt.ylim(1e-3, 200.0)
    plt.xlabel(r'channel$_y$', fontsize=16)
    plt.ylabel(r'Intensity (arb. u.)', fontsize=16)

    # =======
    # slice along qz
    # =======
    plt.subplot(2, 2, 4)
    expslice = sum_over_xrange(expdata, x_range)
    simslice = sum_over_xrange(simdata, x_range)
    plt.semilogy(simslice.getBinCenters(), settings['f']*simslice.getBinValues() + noise, color='r', linewidth=2)
    plt.semilogy(expslice.getBinCenters(), expslice.getBinValues(), color='k', marker='.', linestyle='None')
    plt.xlim(expslice.getXmin(), expslice.getXmax())
    plt.ylim(1e-3, 200.0)
    plt.xlabel(r'channel$_z$', fontsize=16)
    plt.ylabel(r'Intensity (arb. u.)', fontsize=16)

    plt.show()

if __name__ == '__main__':
    plot_data(settings=settings_hexagonal)
    # plot_data(settings=settings_fcc)
