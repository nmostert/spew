import tephra2io
import gridutils as grd
import os
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon, MultiPoint
from scipy.spatial import Voronoi  # , voronoi_plot_2d
import shapely.geometry
import shapely.ops
from matplotlib.colors import LogNorm
from matplotlib import rc


rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


class Eruption:
    """Describes a set of data for a single eruption."""

    def __init__(self, data=None, vent=None, test=False, config=None):
        """Construct Eruption Object.

        This aims to be the constructor of the new Eruption class.

        """
        if os.path.isfile(data):
            self.df, self.phi_labels = tephra2io.tephra2_to_df(data)
        elif isinstance(data, GeoDataFrame):
            self.df = data
        elif isinstance(data, pd.DataFrame):
            self.df = data
        else:
            print("ERROR: Invalid or no data. ")

        if test:
            self.df = self.df.head()

        if isinstance(vent, Point):
            self.vent = vent
            self.grid_vent = grd.nearest_grid_point(
                self.df, vent)
            self.vent_data = grd.nearest_data_point(self.df, vent)

            radii = self.df.apply(
                lambda row: self.vent.distance(row['geometry']), axis=1)
            self.df = self.df.assign(radius=radii)
        else:
            print("ERROR: No vent.")

        if isinstance(config, dict):
            self.vent_easting = config["vent_easting"]
            self.vent_northing = config["vent_northing"]
            if self.vent is None:
                self.vent = Point(config["vent_easting"],
                                  config["vent_northing"])
            self.vent_elevation = config["vent_elevation"]
            self.plume_height = config["plume_height"]
            self.alpha = config["alpha"]
            self.beta = config["beta"]
            self.eruption_mass = config["eruption_mass"]
            self.max_grainsize = config["max_grainsize"]
            self.min_grainsize = config["min_grainsize"]
            self.median_grainsize = config["median_grainsize"]
            self.std_grainsize = config["std_grainsize"]
            self.eddy_const = config["eddy_const"]
            self.diffusion_coefficient = config["diffusion_coefficient"]
            self.fall_time_threshold = config["fall_time_threshold"]
            self.lithic_density = config["lithic_density"]
            self.pumice_density = config["pumice_density"]
            self.col_steps = config["col_steps"]
            self.part_steps = config["part_steps"]
            self.plume_model = config["plume_model"]

    def sample(self, n, weight="Mass/Area", alpha=0.5, beta=0.1, filters=None,
               scorefun=None):
        """Sample the dataframe.

        Each point is assigned a uniform random number, which is weighted
        according to the score function. The n points with the highest assigned
        scores are returned.

        Keyword arguments:
        weight   -- the name of the column by which the sampling probability
                    should be weighted (default None)
        alpha    -- the sampling factor (default 1.0)
        filters  -- dictionary of the form {key:tuple} where key is the column
                    to be filtered by, and the tuple is the lower and upper
                    bounds of the accepted range. (default None)
                    Eg: {"radius": (1000, 500000)} will exclude all points with
                    radii not in the (1000,500000) range.
        scorefun -- user defined score function. Should be of the form
                    f(w, a, b), where w is the weight and a and b are
                    parameters.
                    Default w^(b*log(1+alpha)).

        """
        temp_df = self.df.copy()

        if filters is not None:
            for col in filters:
                temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                                  (temp_df[col] < filters[col][1])]
        if scorefun is None:
            def scorefun(w, alpha, beta):
                r = np.random.random_sample()
                return r * (w ** (beta * np.log(1 + alpha)))
        else:
            def scorefun(w, alpha, beta):
                r = np.random.random_sample()
                return r * scorefun(w, alpha, beta)

        temp_df['norm_weight'] = temp_df[weight].values / \
            np.max(self.df[weight].values)

        probs = temp_df.apply(lambda row: scorefun(
            row['norm_weight'], alpha, beta), axis=1)

        probs = pd.Series(probs, temp_df.index)

        probs.sort_values(inplace=True, ascending=False)

        sample = temp_df.loc[probs[:n].index, :]

        return sample

    def radial_sample(self, azimuth):
        """Sample the dataframe.

        Sample the dataframe radially

        """
        angle = 90 - azimuth  # alpha = 90 - theta

        # finding grid corners (and transforming origin to (0,0))
        x_max = self.df['Easting'].max() - self.vent.x
        x_min = self.df['Easting'].min() - self.vent.x
        y_max = self.df['Northing'].max() - self.vent.y
        y_min = self.df['Northing'].min() - self.vent.y

        def pythagoras(x, y):
            return np.sqrt(x**2 + y**2)

        # finding longest possible distance from vent on grid
        # i.e. distance from vent to furthest corner
        r = max([pythagoras(x_max, y_max), pythagoras(x_max, y_min),
                 pythagoras(x_min, y_min), pythagoras(x_min, y_max)])

        # creating array of distances along radial sampling line
        # TODO: Allow user specification of amount of samples
        # Possibly by clipping line first, then creating this linspace.
        rr = np.linspace(0, r, 100)

        # calculating coordinates at each point along sampling line
        coords = []
        for r in rr:
            coords.append((self.vent.x + r * np.cos(angle * (np.pi / 180)),
                           self.vent.y + r * np.sin(angle * (np.pi / 180))))

        # Creating shapely LinString object to allow clipping
        line = LineString(coords)
        xl, yl = line.xy

        # Creating rectangular representation of grid boundaries
        grid_poly = Polygon([(x_max + self.vent.x, y_max + self.vent.y),
                             (x_max + self.vent.x, y_min + self.vent.y),
                             (x_min + self.vent.x, y_min + self.vent.y),
                             (x_min + self.vent.x, y_max + self.vent.y)])

        xp, yp = grid_poly.exterior.xy

        # Intersecting radial line and grid to clip line to grid
        grid_line = line.intersection(grid_poly)
        xi, yi = grid_line.xy

        # Iterating along line, and finding nearest data point to each point
        near_points = []
        for point in grid_line.coords:
            near = shapely.ops.nearest_points(
                Point(point), MultiPoint(self.df.geometry.values))
            near_points.append(near[1])

        # Selecting each datapoint from dataframe
        sample = self.df.query("geometry in @near_points")

        sample_points = LineString(near_points)

        xg, yg = sample_points.xy

        # plot_grid(self.df, self.vent)
        plt.plot(xp, yp, 'b', label='Grid Boundaries')
        plt.plot(xg, yg, 'g', marker='o', ms=3, label='Radial Samples')
        plt.plot(xl, yl, 'r--', label='Sampling line')
        plt.legend(loc='upper left')
        plt.title(r'Radial Sampling Algorithm with Azimuth $\alpha$ = %d' % azimuth)
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
        plt.tight_layout()
        plt.savefig("rad_samp.eps", dpi=200, format='eps')
        return sample


def plot_GS(df, phi_labels, point):
    data = df[df.geometry == grd.nearest_grid_point(df, point)]
    phi_vals = data[phi_labels].transpose()
    phi_vals = phi_vals / 100
    phi_vals.plot(kind='bar', rot=0, legend=False)


def plot_grid(df, vent=None, labels=None, save=None):
    xx = df['Easting'].values
    yy = df['Northing'].values

    fig, ax = plt.subplots()
    ax.axis('equal')
    ax.plot(xx, yy, 'k.', ms=3)
    if vent is not None:
        ax.plot(vent.coords[0][0], vent.coords[0][1], 'r^', ms=6)
    if labels is None:
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
    else:
        plt.xlabel(labels[0])
        plt.ylabel(labels[0])
    if save is not None:
        plt.savefig(save, dpi=200, format='eps')
    return fig, ax


def plot_contour(xx, yy, zz, vent, title, cbar_label, cmap='', save=None):
    fig, ax = plt.subplots(1, 1)
    ax.axis('equal')
    plt.ylabel("Northing")
    plt.xlabel("Easting")
    plt.title(title)
    zz = zz + 0.000000000000000000001
    cont = ax.tricontourf(xx, yy, zz,
                          norm=LogNorm(vmin=zz.min(), vmax=zz.max()))
    cbar = fig.colorbar(cont, ax=ax)
    cbar.ax.set_ylabel(cbar_label, rotation=270)
    plt.plot(vent.x, vent.y, 'r^', ms=10)
    plt.tight_layout()
    plt.gca().set_xlim(right=max(xx), left=min(xx))
    plt.gca().set_ylim(bottom=min(yy), top=max(yy))
    if save is not None:
        plt.savefig(save, dpi=200, format='eps')
    return fig, ax


def get_voronoi(df):
    points = [
        row.geometry.coords[0]
        for i, row in df.iterrows()]

    vor = Voronoi(points)

    lines = [
        shapely.geometry.LineString(vor.vertices[line])
        for line in vor.ridge_vertices
        if -1 not in line]
    polys = shapely.ops.polygonize(lines)
    polys = [poly for poly in polys]
    df = pd.DataFrame()
    crs = {'init': 'epsg:4326'}
    vor_df = GeoDataFrame(df, crs=crs, geometry=polys)

    return vor_df


def random_sample(df, frac=0.1):
    # Takes a uniform random sample from a dataset
    # frac is the proportion of the total number of data points to be sampled.
    return df.sample(frac=0.05, replace=False)


def regular_sample(df, n=5):
    # Resamples a dataframe by taking every 5 samples in either direction
    # Meant to be used on a regular grid, but should work with irregular data
    resample = get_index_grid(df).iloc[::n, ::n]
    return df.iloc[resample.values.flatten(), :]


def get_index_grid(df):
    # Creates dataframe that is the same shape as the spatial grid
    # Makes spatial iteration easier on regular grids
    return df.reset_index().pivot(
        index='Northing', columns='Easting')['index']


def nearest_grid_point(df, point):
    # Finds existing data point in df that is closest to given point
    # TODO: rewrite using GeoPandas and shapely.geometry.distance()
    return df.ix[df['geometry'].distance(point).argsort()[0]]['geometry']


def nearest_data_point(df, point):
    return df.loc[df.index[df['geometry'] == nearest_grid_point(df, point)]]


def construct_grid(vent, north, east, south, west, elevation, spacing):
    x_east = np.linspace(
        vent.coords[0][0], vent.coords[0][0] + (east * spacing), east + 1)
    x_west = np.linspace(
        vent.coords[0][0] - west * spacing, vent.coords[0][0], west + 1)
    x = np.concatenate((x_east[:-1], x_west))
    y_south = np.linspace(
        vent.coords[0][1] - south * spacing, vent.coords[0][1], south + 1)
    y_north = np.linspace(
        vent.coords[0][1], vent.coords[0][1] + north * spacing, north + 1)
    y = np.concatenate((y_south[:-1], y_north))

    xx, yy = np.meshgrid(x, y)

    cols = ['Easting', 'Northing', 'Elev.']

    df = pd.DataFrame(columns=cols)
    for i, c in enumerate(xx):
        for j, p in enumerate(c):

            df = df.append({
                'Easting': int(p),
                'Northing': int(yy[i][j]),
                'Elev.': int(elevation)
            },
                ignore_index=True)

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    df = GeoDataFrame(df, crs=crs, geometry=geometry)
    return df


def write_grid_file(df, filename):
    df.to_csv(filename, sep=' ', columns=['Easting', 'Northing', 'Elev.'],
              index=False, header=False)


def write_config_file(config, filename):
    f = open(filename, "w")
    for key in config:
        print(key.upper() + " " + str(config[key]))
        f.write(key.upper() + " " + str(config[key]) + '\n')
    f.close()


def create_trial_configs(parameter, trial_no, values, config=None, bash=True):
    if config is None:
        config = {
            "vent_easting": 532290,
            "vent_northing": 1382690,
            "vent_elevation": 614,
            "plume_height": 25000,
            "alpha": 1.04487,
            "beta": 1.46425,
            "eruption_mass": 1.882315e+12,
            "max_grainsize": -4,
            "min_grainsize": 4,
            "median_grainsize": 0,
            "std_grainsize": 1.8,
            "eddy_const": 0.04,
            "diffusion_coefficient": 0.5,
            "fall_time_threshold": 288,
            "lithic_density": 2700.0,
            "pumice_density": 1000.0,
            "col_steps": 75,
            "part_steps": 8,
            "plume_model": 0
        }
    for val in values:
        config[parameter] = val
        fn = "./data/%s_trial_%d/%s_trial_%d.conf" % (
            parameter, trial, parameter, val)
        write_config_file(config, fn)
    if bash:
        fn2 = "./data/%s_trial_%d/%s_trial_%d.sh" % (
            parameter, trial, parameter, trial)
        f2 = open(fn2, "w")
        f2.write("#!/bin/bash\n")
        f2.write("X=%d\n" % trial)
        f2.write("for v in " + ' '.join([str(v) for v in values]) + '\n')
        f2.write("do\n")
        f2.write("./tephra2-2012 inputs/%s_trial_$X/%s_trial_$v.conf " %
                 (parameter, parameter)
                 + "inputs/%s_trial_$X/%s_trial_grid.utm " %
                 (parameter, parameter)
                 + "inputs/noWind > "
                 + "inputs/%s_trial_$X/%s_trial_$v.txt " %
                 (parameter, parameter) + '\n')
        f2.write("done\n")


if __name__ == "__main__":

    trial = 6

    cols = [5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000]

    create_trial_configs("plume_height", trial, cols)

    # test = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    # plot_grid(test.df, test.vent)
    # test.df
    #
    # plot_contour(test.df['Easting'], test.df['Northing'],
    #              test.df['MassArea'], test.vent, 'Mass/Area', 'kg/m$^2$')
    #
    # sample = test.radial_sample(60)
    #
    # len(sample)
    # plt.plot(sample['radius'], sample['MassArea'], 'r.')
    #
    # plt.plot(test.df['radius'], test.df['MassArea'], 'r.')
