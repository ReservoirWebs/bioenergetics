"""Linked Bioenergetics and foraging model

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


from datetime import date, timedelta
from collections import defaultdict
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize, brute
from scipy import stats

from bioenergetics import params

O2CONV = 13560
"""J/gram of O2 in respiration conversions (Elliot and Davison 1975)."""


class InterpolatedFunction(object):
    """Wrapper for an interpolated function

    Given a set of x and y values, create an interpolated function
    which may be accessed by calling the class instance. The class
    will also contain read-only properties that return the original
    data points as well as their domain and range.

    Optional parameters clip_min and clip_max may be given to exclude
    datapoints where the function value exceeds a specified range.

    Example:
    ```
    depths = [0,1,2,3,4,5]
    temperatures = [25,23,21,20,19,17.5]
    temp_fn = InterpolatedDepthFunction(depths, temperatures, clip_max=4)
    t3_5 = temp_fn(3.5)
    ```

    Attributes:
        x: The x-values used to initialize the instance
        y: The y-values used to initialize the instance
        domain: A two-element tuple containing the minimum and maximum x-values
        range: A two-element tuple contianing the minimum and maximum y-values
    """

    def __init__(self, x, y, clip_max=None, clip_min=None):
        """Constructor method

        Args:
            x: A collection of numeric values
            y: A collection of numeric values
            clip_max: An optional numeric value. If given, all x,y pairs
                where y > clip_max will be ignored.
            clip_min: An optional numeric value. If given, all x,y pairs
                where y < clip_min will be ignored

        Returns:
            The initialized instance
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

        self._x = x
        self._y = y
        self._domain = (np.min(x), np.max(x))
        self._range = (np.min(y), np.max(y))

        self.fn = interp1d(x, y, bounds_error=False, fill_value=(y[0], y[-1]))

    def __call__(self, d):
        return self.fn(d)

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def domain(self):
        return self._domain

    @property
    def range(self):
        return self._range


class Model:
    """Linked bioenergetics and foraging model

    This class links the Wisconsin bioenergetics model with a visual
    foraging model to simulate optimal growth for a fish species. As
    inputs, it requires parameters describing temperature depth
    profiles, bathymetry, prey abundance, and species-specific
    bioenergetics parameters. See the documentation for
    `bioenergetics.model.Model.__init__` for more details.

    Once instantiated, the `bioenergetics.model.Model.run`
    method will simulate a period of several days (default 30),
    determining the day- and night-time depths which
    optimize the fish's growth during that period.

    """

    def __init__(
        self,
        starting_mass,
        prey_data,
        temp_fn,
        bathymetry_fn,
        day_hours,
        light_extinction,
        params=params.Chinook(),
        depth_max=10000,
        depth_min=-1,
        day_light=39350,
        night_light=0.1,
        surface_elevation=100000,
        allow_dvm=True,
        max_P=1.0,
        allow_functional_response=True,
        force_depth=None,
        infection_stats=True,
    ):
        """Constructor method

        Args:
            starting_mass: The starting mass of the fish in grams.
            prey_data: A `bioenergetics.prey.PreyData` instance.
            temp_fn: A `bioenergetics.model.InterpolatedFunction`
                instance representing temperature in degrees C as a function of
                depth in meters.
            bathymetry_fn: An instance of
                `bioenergetics.model.InterpolatedFunction`
                 representing surface area in meters^2 as a function of
                elevation in meters.
            day_hours: A numeric value specifying the day length in hours.
            light_extinction: The light extinction coefficient.
            params: An optional `bioenergetics.params.FishParameters` instance
                containing the parameters for the fish species being modeled.
                Defaults to an instance of `bioenergetics.params.Chinook`.
            depth_max: An optional numeric value that caps the maximum
                simulated depth. If no value is given, this is equal the
                maximum `domain` value of `temp_fn`
            depth_min: An optional numeric value that caps the minimum
                simulated depth. If no value is given, this is equal the
                minimum `domain` value of `temp_fn`
            day_light: An optional numeric value specifying the surface light
                intensity during day time, in lux. Defaults to 39350.
            night_light: An optional numeric value specifying the surface light
                intensity during night time, in lux. Defaults to 0.1.
            surface_elevation: An optional numeric value specifying the
                elevation of the water's surface in meters. If not given,
                this is derived from the the `domain` value of `bathymetry_fn`.
            allow_dvm: An optional boolean parameter that allows the fish to
                perform diel vertical migration, selecting different depths for
                day and night time periods. Defaults to true.
            allow_functional_response: An optional boolean parameter that
                limits the fish's consumption by functional response.
                Defaults to true.
            max_P: An optional numeric parameter that caps the proportion of
                maximum consumption that is allowed. Defaults to 1.0
            force_depth: An optional parameter which, if given, disables
                optimization for depth and forces the fish to select the
                given depth value when computing growth. May be a pair of
                values [day_depth, night_depth] or a scalar value which will
                be used for both day and night. Defaults to None
            infection_stats: Whether to compute copepod infection stats.
                Defaults to True.
        """
        self.prey_data = prey_data
        self.light_extinction = light_extinction
        self.starting_mass = starting_mass
        self.depth_min = max(depth_min, temp_fn.domain[0])
        self.depth_max = max(min(depth_max, temp_fn.domain[1]), self.depth_min + 0.2)
        self.day_hours = day_hours
        self.depths = []
        self.surface_elevation = min(surface_elevation, bathymetry_fn.domain[1])
        self.params = params

        self.day_light = day_light
        self.night_light = night_light

        self.starting_length = params.length_from_weight(self.starting_mass)
        self.temp_from_depth = temp_fn
        self.area_from_elevation = bathymetry_fn
        self.allow_dvm = allow_dvm
        self.allow_functional_response = allow_functional_response
        self.max_P = max_P
        self.force_depth = force_depth
        self.infection_stats = infection_stats

    def compute_foragingbydepth(self, length, mass, surface_light, depth):
        """Specific encounter rate for a fish of given size at a given depth.

        Foraging from Beauchamps paper. Current reaction distance is from
        Gregory and Northcote 1992.

        Args:
            length: Length of the fish in millimeters.
            mass: Mass of the fish in grams.
            surface_light: Intensity of light at the surface in lux.
            depth: Depth in meters.

        Returns:
            The specific encounter rate in grams/gram/hour.
        """
        light = surface_light * np.exp((-self.light_extinction) * depth)
        depth = depth
        # prey per cc
        prey_density = self.prey_data.prey_count(depth) / 1000000
        lightenergy = light / 51.2

        # Note that reaction distance is in cm
        suspendedsediment = -((np.log10(lightenergy) - 1.045) / (0.0108))
        if suspendedsediment <= 0:
            reactiondistance = 31.64
        else:
            turbidity = 0.96 * np.log10(suspendedsediment + 1) - 0.002
            reactiondistance = 31.64 - 13.31 * turbidity

        # ~1.1 from this paper, 8 based on kokanee (is ~ the
        # median observed for this Chinook study)
        reactiondistance = max(reactiondistance, 1.1)
        swim_speed = self.params.swim_speed * length / 10
        searchvolume = np.pi * (reactiondistance ** 2) * swim_speed
        # prey per hour
        encounter_rate = searchvolume * prey_density * 60 * 60

        # Capping ER based on 2017 Haskell et al.  Haskell equation is
        # in L, daphnia are currently per cc and was per min, convert
        # to hr. Haskell may underestimate maximum, note the high
        # density corresponds to a pt ~48 that is not represented by
        # the 29.858 cap
        if self.allow_functional_response:
            max_er = (
                29.858
                * (prey_density * 1000)
                * ((4.271 + prey_density * 1000) ** (-1))
                * 60
            )
            # print(encounter_rate, max_er)
            if encounter_rate > max_er:
                encounter_rate = max_er

        # use if want to further restrict capture
        # encounter_rate = 0.9 * encounter_rate
        gramsER = encounter_rate * self.prey_data.wet_weight
        return gramsER / mass

    def compute_ft(self, temperature):
        """Temperature dependence of consumption.

        Args:
            temperature: Temperature in degrees C.

        Returns:
            A numeric value representing the temperature-dependent
            coefficient of consumption.
        """
        CQ = self.params["CQ"]
        CTL = self.params["CTL"]
        CTM = self.params["CTM"]
        CTO = self.params["CTO"]
        CK1 = self.params["CK1"]
        CK4 = self.params["CK4"]
        eq = self.params["c_eq"]
        if eq == 1:
            return np.exp(CQ * temperature)
        elif eq == 2:
            V = (CTM - temperature) / (CTM - CTO)
            Z = np.log(CQ) * (CTM - CTO)
            Y = np.log(CQ) * (CTM - CTO + 2)
            X = (Z ** 2 * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400
            return (V ** X) * np.exp(X * (1 - V))
        elif eq == 3:
            G1 = (1 / (CTO - CQ)) * np.log((0.98 * (1 - CK1)) / (CK1 * 0.002))
            G2 = (1 / (CTL - CTM)) * np.log((0.98 * (1 - CK4)) / (CK4 * 0.02))
            L1 = np.exp(G1 * (temperature - CQ))
            L2 = np.exp(G2 * (CTL - temperature))
            K_A = (CK1 * L1) / (1 + CK1 * (L1 - 1))
            K_B = (CK4 * L2) / (1 + CK4 * (L2 - 1))
            return K_A * K_B
        else:
            raise ValueError("Unknown consumption equation type: " + eq)

    def compute_cmax(self, W):
        """Maximum consumption for a given weight.

        Args:
            W: mass in grams.

        Returns:
            Maximum specific feeding rate in grams/gram/day
        """
        CA = self.params["CA"]
        CB = self.params["CB"]
        return CA * (W ** CB)

    def compute_waste(self, consumption, P, temperature):
        """Waste losses (egestion and excretion)

        Args:
            consumption: specific consumption in grams/gram/day.
            P: feeding level.
            temperature: temperature in degrees C.

        Returns:
            A two-element tuple containing specific egestion (fecal waste),
            and excretion (nitrogenous waste) rates in grams/gram/day.
        """
        # Units are g/g/d
        FA = self.params["FA"]
        FB = self.params["FB"]
        FG = self.params["FG"]
        UA = self.params["UA"]
        UB = self.params["UB"]
        UG = self.params["UG"]
        eq = self.params["egexeq"]
        if eq == 1:
            egestion = FA * consumption
            excretion = UA * (consumption - egestion)
            return (egestion, excretion)
        elif eq == 2:
            egestion = FA * (temperature ** FB) * np.exp(FG * P) * consumption
            excretion = (
                UA * (temperature ** UB) * np.exp(UG * P) * (consumption - egestion)
            )
            return (egestion, excretion)
        elif eq == 3:
            if self.prey_data.indigestibility is None:
                raise ValueError("Prey indigestibility not defined")
            PFF = self.prey_data.indigestibility
            PE = FA * (temperature ** FB) * np.exp(FG * P)
            PF = ((PE - 0.1) / 0.9) * (1 - PFF) + PFF
            egestion = PF * consumption
            excretion = (
                UA * (temperature ** UB) * np.exp(UG * P) * (consumption - egestion)
            )
            return (egestion, excretion)
        else:
            raise ValueError("Unknown egestion/excretion equation type: " + eq)

    def compute_respiration(self, W, temperature):
        """Respiration losses

        Args:
            W: mass in grams
            temperature: temperature in degrees C

        Returns:
            Specific respiration rate in grams/gram/day
        """
        RA = self.params["RA"]
        RB = self.params["RB"]
        RQ = self.params["RQ"]
        RTO = self.params["RTO"]
        RTM = self.params["RTM"]
        RTL = self.params["RTL"]
        RK1 = self.params["RK1"]
        RK4 = self.params["RK4"]
        ACT = self.params["ACT"]
        BACT = self.params["BACT"]
        eq = self.params["respeq"]
        if eq == 1:
            if temperature > RTL:
                VEL = RK1 * W ** RK4
                print(
                    "SOME OF THE INCLUDED TEMPERATURES ARE LETHAL, "
                    "PLEASE MODIFY THE TEMPERATURE TO EXCLUDE "
                    "TEMPERATURES OVER 25C!"
                )
            else:
                VEL = ACT * (W ** RK4) * np.exp(BACT * temperature)
            FTmetabolism = np.exp(RQ * temperature)
            activity = np.exp(RTO * VEL)
        elif eq == 2:
            Vresp = (RTM - temperature) / (RTM - RTO)
            Zresp = np.log(RQ) * (RTM - RTO)
            Yresp = np.log(RQ) * (RTM - RTO + 2)
            Xresp = (((Zresp ** 2) * (1 + (1 + 40 / Yresp) ** 0.5)) ** 2) / 400
            FTmetabolism = (Vresp ** Xresp) * np.exp(Xresp * (1 - Vresp))
            activity = ACT
        else:
            raise ValueError("Unknown respiration equation type: " + eq)
        return RA * (W ** RB) * FTmetabolism * activity

    def energy_density(self, W):
        """Predator energy density as a function of mass.

        Args:
            W: mass of fish in grams

        Returns:
            Energy density in joules/gram (wet mass).
        """
        if self.params["prededeq"] == 1:
            return self.energy_density
        elif self.params["prededeq"] == 2:
            if W < self.params["cutoff"]:
                alpha = self.params["AlphaI"]
                beta = self.params["BetaI"]
            else:
                alpha = self.params["AlphaII"]
                beta = self.params["BetaII"]
            return alpha + beta * W

    def compute_growth(self, w, egain):
        """Growth as a function of mass and energy gain.

        Args:
            w: mass in grams
            egain: energy gain in joules

        Returns:
            Amount of growth in grams. May be negative.
        """
        if self.params["prededeq"] == 1:
            return w + egain / self.params["energydensity"]
        elif self.params["prededeq"] == 2:
            if w < self.params["cutoff"]:
                alpha = self.params["AlphaI"]
                beta = self.params["BetaI"]
            elif w >= self.params["cutoff"]:
                alpha = self.params["AlphaII"]
                beta = self.params["BetaII"]

            if beta != 0:
                w_new = (
                    -alpha
                    + np.sqrt(alpha ** 2 + 4 * beta * (egain + w * (alpha + beta * w)))
                ) / (2 * beta)
            else:
                w_new = w + egain / alpha

            crossed_cutoff = (
                w < self.params["cutoff"] < w_new or w_new < self.params["cutoff"] < w
            )
            if crossed_cutoff:
                egain_cutoff = self.energy_density(
                    self.params["cutoff"]
                ) - self.energy_density(w)
                return self.compute_growth(self.params["cutoff"], egain_cutoff) + (
                    self.params["cutoff"] - w
                )
            else:
                return w_new - w

    def compute_bioenergetics(self, W, temp, P):
        """Compute components of the Wisconsin bioenergetics model.

        Args:
            W: mass in grams
            temp: temperature in degrees C
            P: feeding level

        Returns:
            A five-element tuple
            `(consumption, egestion, excretion, respiration, SDAction)`
        """
        cmax = self.compute_cmax(W)
        ft = self.compute_ft(temp)
        consumption = cmax * P * ft
        (egestion, excretion) = self.compute_waste(consumption, P, temp)
        respiration = self.compute_respiration(W, temp)
        SDAction = self.params["SDA"] * (consumption - egestion)
        return (consumption, egestion, excretion, respiration, SDAction)

    def best_depth(self, length, mass, x0=None):
        """Solve for optimal day- and night- time depths

        Each day is divided into a day and night period, which entails
        differing light levels affecting foraging efficiency. For each
        period, find the depth that maximizes growth given the
        light-dependent foraging level and metabolic rates derived
        from the temperature at that depth.

        Args:
            length: Length in millimeters
            mass: Mass in grams
            x0: An optional parameter giving an initial guess for the
                optimizer. A two-element tuple containing
                `(day_depth, night_depth)` in meters. If not given, a grid
                search will first be performed to provide the initial guess.

        Returns:
            A two-element tuple `(best_depths, best_results)`
            best_depths: A two element tuple containing the day- and
                night-time depths which produced the optimal growth.
            best_results: Growth and bioenergetics rates corresponding
                to the optimized depths, identical to the output of
                `bioenergetics.model.Model.growth_fn`.

        """
        day_hours = self.day_hours
        night_hours = 24.0 - day_hours

        if self.force_depth:
            if np.isscalar(self.force_depth):
                best_depths = (self.force_depth, self.force_depth)
            elif len(self.force_depth) == 2:
                best_depths = self.force_depth
            else:
                raise ValueError("Invalid value for force_depth")
        else:
            if self.allow_dvm:
                depth_bounds = (
                    (self.depth_min, self.depth_max),
                    (self.depth_min, self.depth_max),
                )

                def objective(depth_pair):
                    (day_depth, night_depth) = depth_pair
                    res = self.growth_fn(
                        day_depth, night_depth, length, mass, day_hours, night_hours
                    )
                    return -res[0]

            else:
                depth_bounds = ((self.depth_min, self.depth_max),)

                def objective(d):
                    res = self.growth_fn(d, d, length, mass, day_hours, night_hours)
                    return -res[0]

            if x0 is None:
                # find an initial guess via grid search
                x0 = brute(objective, depth_bounds)
            elif not self.allow_dvm:
                x0 = [x0[0]]

            res = minimize(
                objective,
                x0=x0,
                method="L-BFGS-B",
                bounds=depth_bounds,
                jac="2-point",
                options={"eps": 1e-3},
            )
            if self.allow_dvm:
                best_depths = res.x
            else:
                best_depths = (res.x[0], res.x[0])

        (day_depth, night_depth) = best_depths
        best_results = self.growth_fn(
            day_depth, night_depth, length, mass, day_hours, night_hours
        )
        return best_depths, best_results

    def growth_fn(self, day_depth, night_depth, length, mass, day_hours, night_hours):
        """Growth over a 24-hour period.

        This function computes growth for a fish of given size,
        located at given day- and night-time depths. It is used as an
        optimization target in the `bioenergetics.model.Model.best_depth`
        method.

        Args:
            day_depth: daytime depth in meters
            night_depth: nighttime depth in meters
            length: length in millimeters
            mass: mass in grams
            day_hours: daytime duration in hours
            night_hours: nighttime duration in hours

        Returns:
            A 9-element tuple `(growth, consumption, egestion, excretion,
            respiration, SDAction, P, day_P, night_P)`

            growth: growth in grams
            P: feeding level for entire 24-hour period
            day_P: daytime feeding level
            night_P: nighttime feeding level

        """
        day_temp = self.temp_from_depth(day_depth)
        night_temp = self.temp_from_depth(night_depth)
        cmax = self.compute_cmax(mass)
        day_foraging = self.compute_foragingbydepth(
            length, mass, self.day_light, day_depth
        )
        night_foraging = self.compute_foragingbydepth(
            length, mass, self.night_light, night_depth
        )
        if day_foraging > night_foraging:
            day_foraging *= day_hours
            day_P = min(day_foraging / cmax, self.max_P)
            night_P = min(self.max_P - day_P, night_foraging * night_hours)
        else:
            night_foraging *= night_hours
            night_P = min(night_foraging / cmax, self.max_P)
            day_P = min(self.max_P - night_P, day_foraging * day_hours)
        print(day_foraging * mass / self.prey_data.wet_weight)

        day_bioe = self.compute_bioenergetics(mass, day_temp, day_P)
        night_bioe = self.compute_bioenergetics(mass, night_temp, night_P)
        day_bioe = np.array(day_bioe) * day_hours / 24.0
        night_bioe = np.array(night_bioe) * night_hours / 24.0
        (consumption, egestion, excretion, respiration, SDAction) = (
            day_bioe + night_bioe
        )
        P = day_P + night_P
        consumptionjoules = consumption * self.prey_data.energy
        egain = (
            consumptionjoules
            - (
                (egestion + excretion + SDAction) * self.prey_data.energy
                + respiration * O2CONV
            )
        ) * mass
        growth = self.compute_growth(mass, egain)
        return (
            growth,
            consumption,
            egestion,
            excretion,
            respiration,
            SDAction,
            P,
            day_P,
            night_P,
        )

    def sustainability_estimate(self, depth, consumed):
        """Estimate the number of fish which can grow optimally at a given
        consumption rate, constrained by the prey abundance at the
        given depth.

        Args:
            depth: depth in meters
            consumed: number of individual prey consumed per day

        Returns:
            Estimated number of fish which can grow optimally.
        """
        elevation = self.surface_elevation - depth
        total_prey = self.prey_data.prey_count(depth)
        area = self.area_from_elevation(elevation)
        consumable = area * total_prey * 0.58
        pop_est = consumable / (consumed * 4)
        return pop_est

    def with_infection_stats(self, results):
        # compute time at surface or 16-17 degrees
        infect_hours = 0
        relative_infect_hours = 0
        passthroughs = 0
        dts = results["day_temperature"]
        nts = results["night_temperature"]
        dds = results["day_depth"]
        nds = results["night_depth"]
        for i, daytemp in enumerate(results["day_temperature"]):
            nighttemp = nts[i]
            daydepth = dds[i]
            nightdepth = nds[1]
            if 16 <= daytemp <= 17 or daydepth <= 1:
                infect_hours += self.day_hours
            if 16 <= nighttemp <= 17 or nightdepth <= 1:
                infect_hours += 24 - self.day_hours
            day_convert = self.day_hours/np.round(stats.norm.pdf(16.5, loc=16.5, scale=1), decimals=3)
            day_prob = stats.norm.pdf(daytemp, loc=16.5, scale=1)
            day_prob = np.round(day_prob, decimals=3)
            relative_infect_hours += np.round(day_prob * day_convert, decimals=3)

            night_convert = (24 - self.day_hours)/np.round(stats.norm.pdf(16.5, loc=16.5, scale=1), decimals=3)
            night_prob = np.round(stats.norm.pdf(nighttemp, loc=16.5, scale=1), decimals=3)
            relative_infect_hours += np.round(night_prob * night_convert, decimals=3)

            if daytemp > 16.5 and nighttemp <= 16.5:
                passthroughs += 1
            if i > 1 and daytemp > 16.5 and nts[i - 1] <= 16.5:
                passthroughs += 1

        night_hours = 24.0 - self.day_hours
        atus = sum(dts) * (self.day_hours / 24.0) + sum(nts) * (night_hours / 24.0)

        size_risk = [1.01 ** (length - 60) for length in results["length"]]

        results["infection_stats"] = {
            "atus": atus,
            "passthroughs": passthroughs,
            "infect_hours": infect_hours,
            "weighted_infect_hours": relative_infect_hours,
            "size_risk": size_risk,
            "size_risk_avg": np.average(size_risk),
        }
        return results

    def run(self, n_days=30, start_date=None):
        """Simulate fish growth over a period of several days.

        Beginning with a fish with a mass equal to
        `bioenergetics.model.Model.starting_mass`, call
        `bioenergetics.model.Model.best_depth` for each day, updating
        the fish's mass each time.

        Args:
            n_days: An optional parameter specifying the number of days
                to simulate. Defaults to 30.
            start_date: An optional parameter specifying the start date of
                the simulation. Accepts either a `datetime.date` object or a
                3-element tuple containing `(year, month, day)`

        Returns:
            A dictionary containing simulation outputs. Each entry is a
            list with `n_days` elements, corresponding to the values from
            each day. The output contains the following keys:

            date: If `start_date` is given, a date string in ISO format
                (YYYY-MM-DD), otherwise an integer index starting at 1.
            mass: The mass of the fish in grams at the end of each day.
            length: The length of fish in millimeters at the end of each day.
            day_depth: The daytime depth of the fish in meters.
            night_depth: The nighttime depth of the fish in meters.
            consumption: The specific consumption in grams/gram/day.
            respiration: The specific respiration in grams/gram/day.
            egestion: The specific egestion in grams/gram/day.
            excretion: The specific excretion in grams/gram/day.
            P: The combined feeding rate for the day.
            day_P: The daytime feeding rate.
            night_P: The nighttime feeding rate.
            day_temperature: The temperature corresponding to `day_depth`.
            night_temperature: The temperature corresponding to `night_depth`.
            daily_consumption: The number of individual prey items consumed

        Raises:
            AssertionError: n_days is not a positive integer

        """
        assert n_days >= 1, "n_days must be a positive integer"

        out = defaultdict(list)
        condition1 = float(
            100 * self.starting_mass * ((self.starting_length / 10) ** (-3.0))
        )
        last_best_depths = None
        length = self.starting_length
        mass = self.starting_mass

        use_date = start_date is not None
        if use_date and not isinstance(start_date, date):
            start_date = date(*start_date)
        ONE_DAY = timedelta(days=1)

        for d in range(n_days):
            best_depths, best_results = self.best_depth(length, mass, last_best_depths)
            (day_depth, night_depth) = best_depths
            last_best_depths = best_depths
            (
                growth,
                consumption,
                egestion,
                excretion,
                respiration,
                SDAction,
                P,
                day_P,
                night_P,
            ) = best_results
            day_temp = self.temp_from_depth(day_depth)
            night_temp = self.temp_from_depth(night_depth)
            pw = self.prey_data.wet_weight
            dailyconsume = (consumption * mass) / pw
            mass += growth
            if growth > 0:
                length = self.params.length_from_weight(mass)

            if use_date:
                datestamp = (start_date + ONE_DAY * d).isoformat()
                out["date"].append(datestamp)
            else:
                out["date"].append(int(d + 1))
            out["day_depth"].append(day_depth)
            out["night_depth"].append(night_depth)
            out["growth"].append(growth)
            out["mass"].append(mass)
            out["length"].append(length)
            out["egestion"].append(egestion)
            out["excretion"].append(excretion)
            out["consumption"].append(consumption)
            out["p"].append(P)
            out["day_temperature"].append(day_temp)
            out["night_temperature"].append(night_temp)
            out["day_P"].append(day_P)
            out["night_P"].append(night_P)
            out["daily_consumption"].append(dailyconsume)

            # early stopping condition
            if mass < 0:
                break

        PopEst = self.sustainability_estimate(day_depth, dailyconsume)
        condition = float(100 * (mass - self.starting_mass) * ((length / 10) ** (-3.0)))
        if self.infection_stats:
            out = self.with_infection_stats(out)
        return (out, dailyconsume, condition, condition1, PopEst, day_P, night_P)


# documentation overrides
__pdoc__ = {
    "InterpolatedFunction.x": False,
    "InterpolatedFunction.y": False,
    "InterpolatedFunction.domain": False,
    "InterpolatedFunction.range": False,
}
