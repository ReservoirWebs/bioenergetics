import numpy as np
from scipy.integrate import trapz
from scipy.interpolate import interp1d

class PreyData(object):
    def __init__(abundance, weight, depths, counts):
        self.abundance = abundance
        self.weight = weight
        self.depths = depths
        self.counts = counts
        self.compute_depth_profile()

    def compute_depth_profile(self):
        x = self.depths
        y = self.counts
        surface_count = y[np.argmin(x)]
        auc = trapz(x,y)
        y = y / auc*self.abundance
        self.depth_fn = interp1d(x, y, bounds_error=False,
                                 fill_value=surface_count)

    def prey_count(self,d):
        return self.depth_fn(d)

def daphnia_weight(size, kind='smirnov'):
    """Compute weight for daphnia from size in mm.

    Three weight relationships are available, selected by the 'kind' parameter
    """
    if kind == 'smirnov':
        #Wet weight from Smirnov 2014 (g from mg)
        return (0.075 * size ** 2.925) / 1000
    elif kind == 'cornell':
        return (np.exp(1.468 + 2.83 * np.log(size))) / 1000000

class DaphniaData(PreyData):
    def __init__(abundance, size, depths, counts):
        super.__init__(abundance,size,depths,counts)


        self.weight =

        # Based off Cornell equation (g from ug)


        # Using Pechen 1965 fresh weight / length relationship
        # reported in Dumont for D. magna
        # self.daphnia_weight = (0.052* self.daphnia_size ** 3.012) / 1000
