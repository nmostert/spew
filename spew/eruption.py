import tephra2io
import os
from geopandas import GeoDataFrame
import pandas as pd
from shapely.geometry import Point


class Eruption:
    """Describes a set of data for a single eruption."""

    def __init__(self, data, vent=None):
        """Construct Eruption Object.

        This aims to be the constructor of the new Eruption class.

        """
        if os.path.isfile(data):
            self.df = tephra2io.tephra2_to_df(data)
            self.gdf = tephra2io.to_gdf(self.df)
        elif isinstance(data, GeoDataFrame):
            self.gdf = data
        elif isinstance(data, pd.DataFrame):
            self.df = data

        if isinstance(vent, Point):
            self.vent = vent


if __name__ == "__main__":
    filename = './data/cerroNegro_radGrid_noWind_ONLINE.txt'

    cn = Eruption(filename)
