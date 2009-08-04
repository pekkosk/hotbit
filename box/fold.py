import numpy as nu


def gaussian_peak(grid, x0, width):
    """ Return a gaussian peak with FWHM=width in the given grid. """
    return nu.exp( - ( (grid-x0) / (2*width) )**2 )


def fold(grid, values, width, weights=None):
    """ Fold given values with gaussian distributions over
        given grid. Width is the FWHM of the peak, and values
        can be provided with different weights (default value
        is that the weights are 1 for all values)."""
    if weights == None:
        weights = nu.ones_like(values)
    ret = nu.zeros_like(grid)
    for val, weight in zip(values, weights):
        ret += weight * gaussian_peak(grid, val, width)
    return ret


def make_cumulative_plot(grid, values, width, weights_list, labels=None, colors=None):
    """ Make a cumulative plot from the values that are folded with
        gaussian distributions of FWHM=width. The weights_list is an
        array of the shape [:,len(values)]. """
    import pylab
    if labels == None:
        labels = ['_nolegend_' for i in range(len(weights_list))]
    if colors == None:
        colors = nu.random.random((len(values), 4))
        colors[:,3] = 1
    low = nu.zeros_like(grid)
    for weights, color, label in zip(weights_list, colors, labels):
        up = low + fold(grid, values, width, weights)
        x, y = pylab.poly_between(grid, low, up)
        pylab.fill(x, y, facecolor=color, edgecolor='none', label=label)
        low = up

