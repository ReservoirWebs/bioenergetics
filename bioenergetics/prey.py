""" Classes for encapsulating prey data.

    Specific prey species should inherit from the PreyData base class.

    This file is part of GrowChinook.

    GrowChinook is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    GrowChinook is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GrowChinook.  If not, see <https://www.gnu.org/licenses/>.
"""


import numpy as np
from scipy.integrate import trapz
from scipy.interpolate import interp1d


class PreyData(object):
    def __init__(self, abundance, depths, counts, energy=None,
                 indigestibility=None, size=None, weight=None, *args, **kwargs):
        self.abundance = abundance
        self.depths = depths
        self.min_depth = np.min(depths)
        self.max_depth = np.max(depths)
        self.counts = counts
        self.energy = energy
        self.indigestibility = indigestibility
        self._init_extra(*args, **kwargs)

        self.compute_depth_profile()

        if weight:
            self.weight = weight
        elif size:
            self.weight = self.weight_from_size(size)

    def _init_extra(self):
        """
        Subclasses can override this method to support extra
        __init__ arguments.
        """

        pass

    def compute_depth_profile(self):
        x = self.depths
        y = self.counts
        surface_count = y[np.argmin(x)]
        auc = trapz(y, x)
        y = y / auc*self.abundance
        self.depth_fn = interp1d(x, y, bounds_error=False,
                                 fill_value=surface_count)

    def prey_count(self, d):
        return self.depth_fn(d)

    def weight_from_size(self, size):
        pass


class DaphniaData(PreyData):
    def _init_extra(self, weight_eq='smirnov'):
        if weight_eq.lower() in ['smirnov', 'cornell', 'pechen']:
            self.weight_eq = weight_eq
        else:
            raise ValueError('kind must be one of: ' +
                             '"smirnov", "cornell", "pechen"')

        # From Luecke and Brandt 22.7 overall, 23.3 kJ/g for unfrozen
        # Daphnia (dry weight) 1.62 kJ/g wet weight
        if self.energy is None:
            self.energy = 22700

        # Noue and Choubert 1985 suggest Daphnia are 82.6% digestible
        # by Rainbow Trout
        if self.indigestibility is None:
            self.indigestibility = 0.174

    def weight_from_size(self, size):
        """Compute weight for daphnia from size in mm.

        Three weight relationships are available, selected by the
        'weight_eq' parameter
        """
        if self.weight_eq == 'smirnov':
            # Wet weight from Smirnov 2014 (g from mg)
            return (0.075 * size ** 2.925) / 1000
        elif self.weight_eq == 'cornell':
            # Based off Cornell equation (g from ug)
            return (np.exp(1.468 + 2.83 * np.log(size))) / 1000000
        elif self.weight_eq == 'pechen':
            # Using Pechen 1965 fresh weight / length relationship
            # reported in Dumont for D. magna
            return (0.052 * size ** 3.012) / 1000
