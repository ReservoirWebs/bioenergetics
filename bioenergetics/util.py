from csv import DictReader, DictWriter
import json
import os

from scipy.interpolate import interp1d
from scipy.integrate import trapz
import numpy as np
from matplotlib import pyplot


def select_rows(csvfile, x_key, y_key, site=None, month=None, year=None):
    """From an input csv file, select rows for a given site and season."""

    site = site and str(site)
    month = month and str(month)
    year = year and str(year)
    with open(csvfile) as fid:
        reader = DictReader(fid)
        rows = [
            r
            for r in reader
            if (
                (site is None or r["site"].lower() == site.lower())
                and (month is None or r["month"].lower() == month.lower())
                and (year is None or r["year"] == year)
            )
        ]
        x = [float(r[x_key]) for r in rows]
        y = [float(r[y_key]) for r in rows]
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

    return (interp1d(depths, counts, bounds_error=False, fill_value=surface_count), auc)


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
        raise ValueError("Clip boundaries exclude all datapoints")
    else:
        x = x[idx]
        y = y[idx]

    return interp1d(x, y, bounds_error=False, fill_value=(y[0], y[-1]))


def transpose_dict(d):
    """Convert a dictionary of lists to a list of dictionaries."""
    ks = d.keys()
    vs = list(d.values())
    n = len(vs[0])
    transposed = [dict() for i in range(n)]
    for i in range(n):
        for k in ks:
            transposed[i][k] = d[k][i]
    return transposed


def export_results(results, filename, fmt=None, extra_columns=None):
    """Write model outputs to file

    Args:
        results: A dictionary containing the output from
            `bioenergetics.model.Model.run`.
        filename: A string containing the path to the output file.
        fmt: An optional string specifying the format of the output file.
            Should be either 'json' or 'csv'. If not given, will try to
            guess from the filename extension, defaulting to 'csv' for
            unrecognized extensions.

    Raises:
        AssertionError: when an unrecognized value is passed to `fmt`
        TypeError: extra_columns is not a dictionary or None
    """

    if not fmt:
        path, ext = os.path.splitext(filename)
        if ext == "json":
            fmt = "json"
        else:
            fmt = "csv"
    else:
        assert fmt.lower() in ["json", "csv"], 'fmt must be either "json" or "csv"'

    rs = results.copy()
    if extra_columns is not None:
        rs.update(extra_columns)
    rs = transpose_dict(rs)

    print("exporting %d rows to %s" % (len(rs), filename))
    with open(filename, "w") as fid:
        if ext.lower() == "json":
            json.dump(rs, fid)
        else:
            writer = DictWriter(fid, rs[0].keys())
            writer.writeheader()
            writer.writerows(rs)


def plot_results(results, filename=None, title=None):
    fig = pyplot.figure(facecolor="#c8e9b1")
    if title is None:
        title = "Juvenile Spring Chinook"
    fig.suptitle(title, fontsize=20)
    mass_plot = fig.add_subplot(221)
    mass_plot.plot(results["mass"], label="Mass (g)")
    mass_plot.set_ylabel("Mass (g)")
    mass_plot.set_xlabel("Day of Month")

    growth_plot = fig.add_subplot(222)
    growth_plot.plot(results["growth"])
    growth_plot.set_ylabel("Growth (g/g/d)")
    growth_plot.set_xlabel("Day of Month")

    day_depth_plot = fig.add_subplot(223)
    day_depth_plot.plot(results["day_depth"], "black", label="Day Depth (m)")
    day_depth_plot.set_ylabel("Day Depth (m)")
    day_depth_plot.set_xlabel("Day of Month")
    day_depth_plot.set_ylim(35, 0)
    day_depth_plot.yticklabels = np.arange(0, 35, 5)

    night_depth_plot = fig.add_subplot(224)
    night_depth_plot.set_ylabel("Night Depth (m)")
    night_depth_plot.set_xlabel("Day of Month")
    night_depth_plot.plot(results["night_depth"], "black", label="Night Depth (m)")
    night_depth_plot.yticklabels = np.arange(0, 35, 5)
    night_depth_plot.set_ylim(35, 0)

    pyplot.subplots_adjust(top=0.3)
    fig.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)

    if filename:
        pyplot.savefig(filename, facecolor=fig.get_facecolor(), edgecolor="lightblue")
