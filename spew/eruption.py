import tephra2io
import gridutils as grd
import os
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
import visualisation as vis
import matplotlib.pyplot as plt
from shapely.geometry import Point


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

        if test:
            self.df = self.df.head()

        if isinstance(vent, Point):
            self.vent = vent
            self.grid_vent = grd.nearest_grid_point(
                self.df, vent)
            self.vent_data = grd.nearest_data_point(self.df, vent)
        self.phi_labels = ['[-4,-3)',
                           '[-3,-2)',
                           '[-2,-1)',
                           '[-1,0)',
                           '[0,1)',
                           '[1,2)',
                           '[2,3)',
                           '[3,4)']

    def sample(self, weight="Mass/Area", weight_thres=None, alpha=1, exclude=1000):
        """Sample data.

        Rejection sampling, optionally weighted by Mass/Area

        """
        temp_df = self.df.copy()

        # Pre-filtering
        if weight is not None:
            if weight_thres is None:
                bot = temp_df[weight].values.min()
                top = temp_df[weight].values.max()
                weight_thres = (bot, top)
            temp_df = temp_df[(temp_df[weight] > weight_thres[0]) &
                              (temp_df[weight] < weight_thres[1])]

        if exclude is not None:
            temp_df = temp_df[temp_df.geometry.distance(cn.vent) > exclude]

        prob_arr = temp_df.apply(
            lambda row: np.random.random_sample() * (row[weight]**alpha), axis=1)

        prob_arr = prob_arr / prob_arr.max()
        # print(np.median(prob_arr))
        sample = temp_df.loc[prob_arr > np.median(prob_arr), :]
        # print(sample)
        return sample


if __name__ == "__main__":
    filename = './data/cerroNegro_radGrid_noWind_ONLINE.txt'

    cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    sample = cn.sample(weight_thres=(0.000001, 10000), alpha=0.1)
    print(len(sample))
    vis.plot_grid(sample, vent=cn.vent)
    plt.show()
