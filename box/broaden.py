import numpy as np


def gaussian_peak(grid, x0, FWHM):
    """ Return a normalized gaussian peak with FWHM=width
        in the given grid. """
    a = np.sqrt(4*np.log(2)/np.pi)/FWHM
    b = 4*np.log(2)/FWHM**2
    return a*np.exp( - b*(grid-x0)**2 )


def broaden(grid, values, width, weights=None):
    """ Broaden given values with Gaussian distributions over
        given grid. Width is the FWHM of the peak, and values
        can be provided with different weights (default value
        is that the weights are 1 for all values). Each peak
        integrates to the corresponding weight."""
    if weights == None:
        weights = np.ones_like(values)
    ret = np.zeros_like(grid)
    for val, weight in zip(values, weights):
        ret += weight * gaussian_peak(grid, val, width)
    return ret


def make_cumulative_plot(grid, values, width, weights_list, labels=None, colors=None):
    """ Make a cumulative plot from the values that are broadened with
        gaussian distributions of FWHM=width. The weights_list is an
        array of the shape [:,len(values)]. """
    import pylab
    if labels == None:
        labels = ['_nolegend_' for i in range(len(weights_list))]
    if colors == None:
        colors = np.random.random((len(values), 4))
        colors[:,3] = 1
    low = np.zeros_like(grid)
    for weights, color, label in zip(weights_list, colors, labels):
        up = low + broaden(grid, values, width, weights)
        ## MS: incompatibility issue with matplotlib>=3.1
#        x, y = pylab.poly_between(grid, low, up)
#        pylab.fill(x, y, facecolor=color, edgecolor='none', label=label)
        pylab.fill(np.append(grid,low), np.append(up, low),
                   facecolor=color, edgecolor='none', label=label)
        low = up

