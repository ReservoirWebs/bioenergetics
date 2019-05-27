import json
from datetime import date
import pytest
from scipy.integrate import trapz
import numpy as np
from bioenergetics import util
from bioenergetics.model import Model, InterpolatedFunction
from bioenergetics.prey import DaphniaData

year = '2015'
month = 6
site = 'Fall Creek'
surface_elevation = 833.3
light_extinction = 0.792

daylengths = [12,12,12,11.83,13.4,14.73,15.42,15.12,13.97,12.45,12,12,12]


@pytest.fixture
def daphnia_data():
    total_daphnia = 5020.65
    daphnia_size = 1.26
    depth, counts = util.select_rows('tests/data/daphnia-vd.csv',
                                     'depth', 'total daphnia',
                                     site, month, year)
    return DaphniaData(total_daphnia, depth, counts, size=daphnia_size)


@pytest.fixture
def bathymetry():
    elev, area = util.select_rows('tests/data/bathymetry.csv',
                                  'elevation (m)', '2d_area (m2)',
                                  site=site)
    return InterpolatedFunction(elev, area)


@pytest.fixture
def temperature():
    depths, temperatures = util.select_rows('tests/data/temperatures.csv',
                                            'depth', 'temperature',
                                            site, month, year)
    return InterpolatedFunction(depths, temperatures)


def regression_test(results, filename):
    try:
        with open(filename) as fid:
            artifact = json.load(fid)
    except FileNotFoundError:
        artifact = None

    if artifact:
        assert results == artifact
    else:
        with open(filename, 'w') as fid:
            json.dump(results, fid, default=lambda x: int(x))


def safe_round(values, decimals):
    try:
        return list(np.round(values, decimals))
    except Exception:
        return values


def truncate(results, decimals=3):
    return {k: safe_round(v, decimals) for k, v in results.items()}


def test_prey_depth_profile(daphnia_data):
    xs = np.arange(daphnia_data.min_depth, daphnia_data.max_depth, 0.5)
    ys = list(map(daphnia_data.prey_count, xs))
    auc = trapz(ys, xs)
    assert auc == 5020.65
    assert (daphnia_data.prey_count(10) <
            daphnia_data.prey_count(0) <
            daphnia_data.prey_count(6))


def test_batch(temperature, bathymetry, daphnia_data):
    starting_mass = 3.0
    batch = Model(starting_mass, daphnia_data, temperature, bathymetry,
                  daylengths[month], light_extinction)
    results = truncate(batch.run(start_date=(2015, 6, 1))[0])
    results_dt = truncate(batch.run(start_date=date(2015,6,1))[0])
    regression_test(results, 'tests/data/fc-201506-regression.json')
    regression_test(results_dt, 'tests/data/fc-201506-regression.json')


def test_no_dvm(temperature, bathymetry, daphnia_data):
    starting_mass = 3.0
    batch = Model(starting_mass, daphnia_data, temperature, bathymetry,
                  daylengths[month], light_extinction, allow_dvm=False)
    results = truncate(batch.run()[0])
    regression_test(results, 'tests/data/fc-2015506-nodvm-regression.json')
