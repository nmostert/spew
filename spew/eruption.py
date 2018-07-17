import tephra2io
import gridutils as grd
import os
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
import visualisation as vis
import matplotlib.pyplot as plt
from shapely.geometry import Point
from scipy.spatial import Voronoi, voronoi_plot_2d
import shapely.geometry
import shapely.ops


class Eruption:
    """Describes a set of data for a single eruption."""

    def __init__(self, data=None, vent=None, test=False):
        """Construct Eruption Object.

        This aims to be the constructor of the new Eruption class.

        """
        if os.path.isfile(data):
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

    def sample(self, n, weight=None, alpha=1, filters=None):
        """Sample the dataframe.

        Each point is assigned a uniform random number, which is weighted as
        by a factor of weight^alpha. The n points with the highest assigned
        scores are returned.

        Keyword arguments:
        weight  -- the name of the column by which the sampling probability
                   should be weighted (default None)
        alpha   -- the sampling factor (default 1.0)
        filters -- dictionary of the form {key:tuple} where key is the column
                   to be filtered by, and the tuple is the lower and upper
                   bounds of the accepted range. (default None)
                   Eg: {"radius": (1000, 500000)} will exclude all points with
                   radii not in the (1000,500000) range.

        """
        temp_df = self.df.copy()

        if filters is not None:
            for col in filters:
                print(filters[col][0])
                temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                                  (temp_df[col] < filters[col][1])]

        def prob(w, alpha):
            return np.random.random_sample() * (w**alpha)

        probs = temp_df.apply(lambda row: prob(row[weight], alpha), axis=1)
        probs = pd.Series(probs, temp_df.index)

        probs.sort_values(inplace=True, ascending=False)

        sample = temp_df.loc[probs[:n].index, :]

        return sample


def plot_GS(df, phi_labels, point):
    data = df[df.geometry == grd.nearest_grid_point(df, point)]
    phi_vals = data[phi_labels].transpose()
    phi_vals = phi_vals / 100
    phi_vals.plot(kind='bar', rot=0, legend=False)


def plot_grid(df, vent=None, labels=None):
    xx = df['Easting'].values
    yy = df['Northing'].values

    fig, ax = plt.subplots()
    ax.axis('equal')
    ax.plot(xx, yy, 'k.', ms=1)
    if vent is not None:
        ax.plot(vent.coords[0][0], vent.coords[0][1], 'r^', ms=3)
    if labels is None:
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
    else:
        plt.xlabel(labels[0])
        plt.ylabel(labels[0])
    return fig, ax


def plot_contour(xx, yy, zz, vent, title, cbar_label, cmap='magma_r', save=None):
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    ax.axis('equal')
    plt.ylabel("Northing")
    plt.xlabel("Easting")
    plt.title(title)
    plt.tricontourf(xx, yy, zz, 20, cmap=cmap)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(cbar_label, rotation=270)
    plt.plot(xx, yy, 'ko', ms=1)
    plt.plot(vent.x, vent.y, 'r^', ms=10)
    if save is not None:
        plt.savefig(save, dpi=200, format='eps')


if __name__ == "__main__":
    filename = './data/cerroNegro_fineGrid_realWind_ONLINE.txt'

    cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    filters = {
        'radius': (3000, np.inf)}
    sample = cn.sample(50, weight='Mass/Area', alpha=0.5, filters=filters)

    fig, ax = vis.plot_grid(sample, vent=cn.vent)
    fig.savefig('samp_1.eps', dpi=200, format='eps')

    plt.show()

    points = [
        row.geometry.coords[0]
        for i, row in sample.iterrows()]

    vor = Voronoi(points)

    # fig, ax1 = plt.subplots(1, 1)

    # vplt = voronoi_plot_2d(vor, ax1)
    # ax1.plot(cn.vent.coords[0][0], cn.vent.coords[0][1], 'r^', ms=10)
    # ax1.set_xlabel("Easting")
    # ax1.set_ylabel("Northing")
    # vplt.savefig('vor_50_05.eps', dpi=200, format='eps')

    lines = [
        shapely.geometry.LineString(vor.vertices[line])
        for line in vor.ridge_vertices
        if -1 not in line]

    len(lines)
    polys = shapely.ops.polygonize(lines)

    polys = [poly for poly in polys]
    print(len(polys))

    print(len(sample))

    df = pd.DataFrame()

    crs = {'init': 'epsg:4326'}
    vor_df = GeoDataFrame(df, crs=crs, geometry=polys)

    vis.plot_grid(sample, vent=cn.vent)
    plt.savefig('smp_50_5.eps', dpi=200, format='eps')
