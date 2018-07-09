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

    def sample(self, weight=None, alpha=1, filters=None):
        """Sample data.

        Rejection sampling, optionally weighted by Mass/Area

        """
        temp_df = self.df.copy()

        if filters is not None:
            for col in filters:
                print(filters[col][0])
                temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                                  (temp_df[col] < filters[col][1])]

        prob_arr = temp_df.apply(
            lambda row: np.random.random_sample() * (row[weight]**alpha), axis=1)

        prob_arr = prob_arr / prob_arr.max()

        sample = temp_df.loc[prob_arr > np.median(prob_arr), :]

        return sample


if __name__ == "__main__":
    filename = './data/cerroNegro_radGrid_noWind_ONLINE.txt'

    cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    filters = {
        'Mass/Area': (0.001, 10000),
        'radius': (5000, np.inf)}
    sample = cn.sample(weight='Mass/Area', alpha=0.1, filters=filters)
    print(len(sample))
    vis.plot_grid(sample, vent=cn.vent)
    plt.show()

    points = [
        row.geometry.coords[0]
        for i, row in sample.iterrows()]

    vor = Voronoi(points)

    voronoi_plot_2d(vor)

    lines = [
        shapely.geometry.LineString(vor.vertices[line])
        for line in vor.ridge_vertices]

    len(lines)
    polys = shapely.ops.polygonize(lines)

    polys = [poly for poly in polys]
    print(len(polys))

    print(len(sample))

    df = pd.DataFrame()

    crs = {'init': 'epsg:4326'}
    vor_df = GeoDataFrame(df, crs=crs, geometry=polys)

    vor_df.plot()
