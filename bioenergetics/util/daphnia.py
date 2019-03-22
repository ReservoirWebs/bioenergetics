from scipy.interpolate import interp1d
from scipy.integrate import trapz
import numpy as np
from csv import DictReader

def select_rows(csvfile, site, month, year):
    """Select rows for a given site and season."""
    with open(csvfile) as fid:
        reader = DictReader(fid)
        rows = [r for r in reader if (r['Site'] == site
                                      and r['Month'] == month
                                      and r['Year'] == year)]
        x = [float(r['Depth']) for r in rows]
        y = [float(r['Total Daphnia']) for r in rows]
        return (x,y)

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
