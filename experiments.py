from importlib import reload
from csv import DictReader, DictWriter
import os

from matplotlib import pyplot
import numpy as np

from bioenergetics.model import Model, InterpolatedFunction
from bioenergetics.prey import DaphniaData
from bioenergetics import util

reload(util)

with open("tests/data/monthly-values-max.csv") as fid:
    reader = DictReader(fid)
    monthly_values = [r for r in reader]


def get_values(site, year, month):
    for row in monthly_values:
        if (
            row["site"] == site
            and row["year"] == year
            and row["month"] == str(month)
        ):

            try:
                return (
                    float(row["daphnia density"]),
                    float(row["daphnia size"]),
                    float(row["light extinction"]),
                )
            except ValueError:
                return None
    return None


def get_daphnia_dist(site, month, year, daphnia_density, daphnia_size):
    depths, total = util.select_rows(
        "tests/data/daphnia-vd.csv", "depth", "total daphnia", site, month, year
    )
    return DaphniaData(daphnia_density, depths, total, size=daphnia_size)


def get_temp_fn(site, month, year):
    depths, temps = util.select_rows(
        "tests/data/temperatures.csv", "depth", "temperature", site, month, year
    )
    return InterpolatedFunction(depths, temps)


daylengths = [
    12,
    12,
    12,
    11.83,
    13.4,
    14.73,
    15.42,
    15.12,
    13.97,
    12.45,
    12,
    12,
    12,
]
days_in_month = [-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
name_template = "%s_%s_%s"


def merge_results(a, b):
    c = {}
    for k in b.keys():
        if k in a:
            c[k] = a[k] + b[k]
        else:
            c[k] = b[k]
    return c


def write_csv(data, filename):
    with open(filename, "w") as fid:
        writer = DictWriter(fid, data[0].keys())
        writer.writeheader()
        writer.writerows(data)


def run_sensitivity(site="fall creek", year="2016"):
    elev, area = util.select_rows(
        "tests/data/bathymetry.csv", "elevation (m)", "2d_area (m2)", site=site
    )
    bathymetry = InterpolatedFunction(elev, area)
    k_results = []
    w_results = []
    ld_results = []
    ad_results = []
    for month in [4, 5, 6, 7, 8, 9]:
        print(month)
        values = get_values(site, year, month)
        daphnia_density, daphnia_size, light_extinction = values

        daphnia = get_daphnia_dist(
            site, month, year, daphnia_density, daphnia_size
        )
        temp_fn = get_temp_fn(site, month, year)
        starting_mass = 10
        for k in np.geomspace(0.01, 10, 15):
            batch = Model(
                starting_mass,
                daphnia,
                temp_fn,
                bathymetry,
                daylengths[month],
                k,
                day_light=39350,
                allow_dvm=True,
                allow_functional_response=False,
                max_P=2.0,
            )
            results = batch.run(n_days=1)[0]
            k_results.append(
                {"month": month, "growth": results["growth"][0], "k": k}
            )
        for w in np.linspace(0.01, 200, 15):
            batch = Model(
                w,
                daphnia,
                temp_fn,
                bathymetry,
                daylengths[month],
                k,
                day_light=39350,
                allow_dvm=True,
                allow_functional_response=False,
                max_P=2.0,
            )
            results = batch.run(n_days=1)[0]
            w_results.append(
                {
                    "month": month,
                    "growth": results["growth"][0],
                    "starting_mass": w,
                }
            )
        for ld in np.geomspace(0.1, 4, 15):
            daphnia = get_daphnia_dist(site, month, year, daphnia_density, ld)
            batch = Model(
                starting_mass,
                daphnia,
                temp_fn,
                bathymetry,
                daylengths[month],
                light_extinction,
                day_light=39350,
                allow_dvm=True,
                allow_functional_response=False,
                max_P=2.0,
            )
            results = batch.run(n_days=1)[0]
            ld_results.append(
                {
                    "month": month,
                    "growth": results["growth"][0],
                    "daphnia_size": ld,
                }
            )

        for ad in np.geomspace(1, 50000, 15):
            daphnia = get_daphnia_dist(site, month, year, ad, daphnia_size)
            batch = Model(
                starting_mass,
                daphnia,
                temp_fn,
                bathymetry,
                daylengths[month],
                light_extinction,
                day_light=39350,
                allow_dvm=True,
                allow_functional_response=False,
                max_P=2.0,
            )
            results = batch.run(n_days=1)[0]
            ad_results.append(
                {
                    "month": month,
                    "growth": results["growth"][0],
                    "daphnia_density": ad,
                }
            )
    write_csv(k_results, "k_sensitivity.csv")
    write_csv(w_results, "starting_mass_sensitivity.csv")
    write_csv(ld_results, "daphnia_size_sensitivity.csv")
    write_csv(ad_results, "daphnia_density_sensitivity.csv")


def run(allow_dvm, max_P, write_out=True):
    pyplot.close("all")
    for site in ["fall creek", "hills creek", "lookout point"]:
        elev, area = util.select_rows(
            "tests/data/bathymetry.csv",
            "elevation (m)",
            "2d_area (m2)",
            site=site,
        )
        bathymetry = InterpolatedFunction(elev, area)
        for year in ["2015", "2016"]:
            starting_mass = 0.3
            site_results = {}
            for month in [3, 4, 5, 6, 7, 8, 9]:
                n_days = days_in_month[month]
                values = get_values(site, year, month)
                if values is None:
                    continue
                print(site, year, month)
                daphnia_density, daphnia_size, light_extinction = values
                depths, total = util.select_rows(
                    "tests/data/daphnia-vd.csv",
                    "depth",
                    "total daphnia",
                    site,
                    month,
                    year,
                )
                daphnia = DaphniaData(
                    daphnia_density, depths, total, size=daphnia_size
                )

                depths, temps = util.select_rows(
                    "tests/data/temperatures.csv",
                    "depth",
                    "temperature",
                    site,
                    month,
                    year,
                )

                temp_fn = InterpolatedFunction(depths, temps)

                start_date = (int(year), month, 1)
                batch = Model(
                    starting_mass,
                    daphnia,
                    temp_fn,
                    bathymetry,
                    daylengths[month],
                    light_extinction,
                    day_light=39350,
                    allow_dvm=allow_dvm,
                    max_P=max_P,
                )
                results = batch.run(n_days=n_days, start_date=start_date)[0]

                starting_mass = max(0.3, results["mass"][-1])
                site_results = merge_results(site_results, results)

            if write_out:
                label = ""
                if not allow_dvm:
                    label += "nodvm"
                if max_P != 1:
                    label += "+maxP%d" % max_P
                label = label.lstrip("+")
                dirname = "experiments/%s/" % (label or "control")
                if not os.path.isdir(dirname):
                    os.makedirs(dirname)
                name = name_template % (site, year, label)
                name = name.rstrip("_")
                extras = {"site": [site] * len(site_results["date"])}
                util.export_results(
                    site_results, dirname + name + ".csv", extra_columns=extras
                )
            print("final length: %0.3f mm" % results["length"][-1])


run(True, 2.0, False)
# run(True, 2.0)
# run(True, 1.0)
# run(False, 1.0)
# run(False, 2.0)

# run_sensitivity()
