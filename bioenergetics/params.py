"""Fish species parameters

    Specific fish species should inherit from the
    `bioenergetics.params.FishParameters` base class.

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

from csv import DictReader, DictWriter, QUOTE_NONNUMERIC

KEYS = ['c_eq', 'CA', 'CB', 'CQ', 'CTO', 'CTM', 'CTL', 'CK1', 'CK4', 'respeq',
        'RA', 'RB', 'RQ', 'RTO', 'RTM', 'RTL', 'RK1', 'RK4', 'ACT', 'BACT',
        'SDA', 'egexeq', 'FA', 'FB', 'FG', 'UA', 'UB', 'UG', 'prededeq',
        'energydensity', 'AlphaI', 'BetaI', 'cutoff', 'AlphaII', 'BetaII',
        'swim_speed']
"""Key names for wrapped parameters"""


class FishParameters(object):
    """Fish species modeling parameters

    This class contains parameter values for the Wisconsin Bioenergetics
    model for various fish species. In general, multiple equations are
    available for each bioenergetics process, and not all parameters will
    be used for a particular species. Refer to Hanson 1997 for the
    definitions of the parameters.

    Additionally, this class contains the function describing the
    relationship between weight and fork length.
    """

    def __init__(self, params={}, *args, **kwargs):
        """Constructor method

        Args:
            params: A dictionary containing the keys found in
                `bioenergetics.params.KEYS`.

        Returns:
            The initialized FishParameters instance.
        """
        self.params = params
        self._init_extra(*args, **kwargs)

        assert int(self.prededeq) in [1, 2], \
            "Predator energy density equation (prededeq) must be 1 or 2."
        assert int(self.egexeq) in [1, 2, 3], \
            "Egestion/excretion equation (exegeq) must be 1, 2, or 3."
        assert int(self.c_eq) in [1, 2, 3], \
            "Consumption equation (c_eq) must be 1, 2, or 3."

    def _init_extra(self):
        """
        Subclasses can override this method to support extra
        __init__ arguments.
        """

        pass

    def __getitem__(self, key):
        return self.params[key]

    def __getattr__(self, name):
        return self.params[name]

    def length_from_weight(self, w):
        """Length in millimeters from weight in grams.
        Should be implemented by subclasses.
        """
        pass

    def to_csv(self, filename):
        """Utility function to export parameters to a csv file.

        Args:
            filename: a string containing the path to the output csv file.

        Returns:
            Nothing
        """
        with open(filename, 'w') as fid:
            writer = DictWriter(fid, KEYS)
            writer.writeheader()
            writer.writerow(self.params)

    def from_csv(self, filename):
        """Utility function to import parameters from a csv file."""
        with open(filename) as fid:
            reader = DictReader(fid, quoting=QUOTE_NONNUMERIC)
            self.params = next(reader)


CHINOOK_DEFAULTS = {'c_eq': 3.0, 'CA': 0.303, 'CB': -0.275, 'CQ': 5.0,
                    'CTO': 15.0, 'CTM': 20.93, 'CTL': 24.05, 'CK1': 0.36,
                    'CK4': 0.53, 'respeq': 1.0, 'RA': 0.00264,
                    'RB': -0.217, 'RQ': 0.06818, 'RTO': 0.0234, 'RTM': 0.0,
                    'RTL': 25.0, 'RK1': 1.0, 'RK4': 0.13, 'ACT': 9.7,
                    'BACT': 0.0405, 'SDA': 0.172, 'egexeq': 3.0, 'FA': 0.212,
                    'FB': -0.222, 'FG': 0.631, 'UA': 0.0314, 'UB': 0.58,
                    'UG': -0.299, 'prededeq': 2.0, 'energydensity': 0.0,
                    'AlphaI': 5764.0, 'BetaI': 0.9862, 'cutoff': 4000.0,
                    'AlphaII': 7602.0, 'BetaII': 0.5266,'swim_speed':2.0}
"""Default bioenergetics parameters for Oncorhynchus tshawytscha."""


class Chinook(FishParameters):
    """Parameters for Oncorhynchus tshawytscha.

    Any bioeenergetics parameters not passed to the constructor will
    be filled with default values.

    Args:
        length_eq: A string value of either "murphy" or "macfarlane",
            specifying which weight to length relationship to use.
            Defaults to "murphy"
    """

    def _init_extra(self, length_eq='murphy'):
        length_eq = length_eq.lower()
        if length_eq in ['murphy', 'macfarlane']:
            self.length_eq = length_eq
        else:
            raise ValueError('length_eq must be "murphy" or "macfarlane"')
        params = CHINOOK_DEFAULTS.copy()
        params.update(self.params)
        self.params = params

    def length_from_weight(self, w):
        """Return the fork length given a weight.

           The constructor argument `length_eq` specifies which
           relationship is used.

           Args:
               w: Weight in grams

           Returns:
               The length in millimeters as a float.
        """
        if self.length_eq == 'murphy':
            # From LP and FC screw trap data (R2 = 0.9933)
            return (w / 0.000004) ** (1 / 3.1776)
        elif self.length_eq == 'macfarlane':
            # weight to fork length (MacFarlane and Norton 2008)
            return (w / 0.0003) ** (1 / 2.217)
