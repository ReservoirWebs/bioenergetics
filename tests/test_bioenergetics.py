import json
from csv import DictReader
import pytest
from scipy.integrate import trapz
import numpy as np
from bioenergetics import util
from bioenergetics.model import Model, InterpolatedFunction
from bioenergetics.prey import DaphniaData

year = '2015'
month = 'June'
site = 'Fall Creek'
surface_elevation = 833.3
light_extinction = 0.792

daylengths = {'March': 11.83, 'April': 13.4, 'May': 14.73, 'June': 15.42,
              'July': 15.12, 'August': 13.97, 'September': 12.45}


@pytest.fixture
def daphnia_data():
    total_daphnia = 5020.65
    daphnia_size = 1.26
    x, y = util.select_rows('tests/data/daphnia-vd.csv', site, month, year)
    return DaphniaData(total_daphnia, x, y, size=daphnia_size)


@pytest.fixture
def bathymetry():
    elev = []
    area = []
    with open('tests/data/fall-creek-bath.csv') as fid:
        reader = DictReader(fid)
        for row in reader:
            elev.append(int(row['elevation (m)']))
            area.append(float(row['2d_area (m2)']))
    return InterpolatedFunction(elev, area)


@pytest.fixture
def temperature():
    with open('tests/data/fall-creek-temp-june-2015.csv') as fid:
        reader = DictReader(fid)
        temperatures = []
        depths = []
        for row in reader:
            temperatures.append(float(row['temp']))
            depths.append(float(row['depth']))
    return InterpolatedFunction(depths, temperatures)


@pytest.fixture
def regression_artifact():
    with open('tests/data/fc-201506-regression.json') as fid:
        return json.load(fid)


def truncate(results, decimals=3):
    return {k: list(np.round(v, decimals)) for k, v in results.items()}


def test_prey_depth_profile(daphnia_data):
    xs = np.arange(daphnia_data.min_depth, daphnia_data.max_depth, 0.5)
    ys = list(map(daphnia_data.prey_count, xs))
    print(ys)
    auc = trapz(ys, xs)
    assert auc == 5020.65
    assert (daphnia_data.prey_count(10) <
            daphnia_data.prey_count(0) <
            daphnia_data.prey_count(6))


def test_batch(temperature, bathymetry, daphnia_data, regression_artifact):
    starting_mass = 3.0
    batch = Model(starting_mass, daphnia_data, temperature, bathymetry,
                  daylengths[month], light_extinction)
    results = truncate(batch.run()[0])
    assert results == regression_artifact
