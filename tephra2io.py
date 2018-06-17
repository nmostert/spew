import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point


def tephra2_to_df(filename):
    df = pd.read_csv(filename, sep=" ", header=0)
    df = df.dropna(axis=1, how='all')
    df = df.fillna(0)
    df = df.rename(columns={'#Easting': 'Easting'})

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    df = df.drop(['Easting', 'Northing'], axis=1)
    crs = {'init': 'epsg:4326'}
    gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
    return gdf
