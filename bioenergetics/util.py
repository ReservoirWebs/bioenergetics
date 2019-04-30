from scipy.interpolate import interp1d
from scipy.integrate import trapz
import numpy as np
from csv import DictReader


def select_rows(csvfile, site, month, year):
    """From an input csv file, select rows for a given site and season."""
    with open(csvfile) as fid:
        reader = DictReader(fid)
        rows = [r for r in reader if (r['Site'] == site
                                      and r['Month'] == month
                                      and r['Year'] == year)]
        x = [float(r['Depth']) for r in rows]
        y = [float(r['Total Daphnia']) for r in rows]
        return (x, y)


def compute_curves(depths, counts, sum_counts=None):
    """Interpolate depth-count data and compute AUC.

    If sum_counts is given, the values in counts are scaled such that
    the AUC = sum_counts.
    """
    surface_count = counts[np.argmin(depths)]

    auc = trapz(counts, depths)
    if sum_counts:
        counts = counts / auc * sum_counts
        auc = trapz(counts, depths)

    return (interp1d(depths, counts,
                     bounds_error=False, fill_value=surface_count),
            auc)


def interpolated_function(x, y, clip_max=None, clip_min=None):
    """Wrapper for an interpolated function

    Given a set of x and y values, return a callable object that
    returns interpolated y values for novel x values.

    Optional parameters clip_min and clip_max may be given to exclude
    datapoints where the function value exceeds a specified range.

    Example:

    depths = [0,1,2,3,4,5]
    temperatures = [25,23,21,20,19,17.5]
    temp_fn = interpolated_function(depths, temperatures, clip_max=4)
    t3_5 = temp_fn(3.5)

    """

    idx = np.argsort(x)
    x = np.sort(x)
    y = np.array(y)[idx]

    # clip values
    idx = np.repeat(True, y.size)
    if clip_max:
        idx = np.logical_and(idx, (y <= clip_max))
    if clip_min:
        idx = np.logical_and(idx, (y >= clip_min))

    if not np.any(idx):
        raise ValueError('Clip boundaries exclude all datapoints')
    else:
        x = x[idx]
        y = y[idx]

    return interp1d(x, y, bounds_error=False, fill_value=(y[0], y[-1]))
