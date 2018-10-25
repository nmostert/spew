import tephra2io
import gridutils as grd
import os
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from scipy.spatial import Voronoi  # , voronoi_plot_2d
import shapely.geometry
import shapely.ops
from matplotlib.colors import LogNorm
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


class Eruption:
    """Describes a set of data for a single eruption."""

    def __init__(self, data=None, vent=None, test=False):
        """Construct Eruption Object.

        This aims to be the constructor of the new Eruption class.

        """
        if os.path.isfile(data):
            print("FILE")
            self.df = tephra2io.tephra2_to_df(data)
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

        self.phi_labels = ['[-4,-3)',
                           '[-3,-2)',
                           '[-2,-1)',
                           '[-1,0)',
                           '[0,1)',
                           '[1,2)',
                           '[2,3)',
                           '[3,4)']

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


if __name__ == "__main__":
    filename = './data/ceroNegro_GITHUB.txt'

    pulu = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    plot_grid(pulu.df, pulu.vent)
    pulu.df

    plot_contour(pulu.df['Easting'], pulu.df['Northing'],
                 pulu.df['MassArea'], pulu.vent, 'Mass/Area', 'kg/m$^2$')
