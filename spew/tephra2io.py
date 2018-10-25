import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point
from eruption import Eruption, plot_contour, plot_grid
import re


def tephra2_to_df(filename):
    # Reads input from Tephra2 into a GeoPandas GeoDataFrame
    # Geopandas is used to enable easy spatial calculations
    # (area, distance, etc.)

    # Extract phi-classes from header and construct col names
    headers = pd.read_csv(filename, sep=" ", header=None, nrows=1)
    phi_names = []
    for name in headers[headers.columns[4:-1]].values[0]:
        m1 = re.search('\[[-+]?[0-9]*\.?[0-9]+', name)
        m2 = re.search('(?<=->)[-+]?[0-9]*\.?[0-9]+\)', name)
        phi_names.append(m1.group(0) + "," + m2.group(0))
    col_names = ["Easting", "Northing", "Elevation",
                 "MassArea"] + phi_names + ["Percent"]

    df = pd.read_csv(filename, sep=" ", header=None,
                     names=col_names, skiprows=1)
    df = df.dropna(axis=1, how='all')
    df = df.fillna(0)

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    df = GeoDataFrame(df.copy(), crs=crs, geometry=geometry)
    return df


def read_CN(filename):
    df = pd.read_csv(filename, header=0)
    df[["Northing", "Easting"]] = df[[
        "Northing", "Easting"]].fillna(method='ffill')

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    df = GeoDataFrame(df.copy(), crs=crs, geometry=geometry)
    return df


if __name__ == "__main__":
    filename1 = './data/ceroNegro_GITHUB.txt'
    cn = Eruption(data=filename1, vent=Point(532290, 1382690), test=False)

    plot_grid(cn.df, cn.vent)
    cn.df
    plot_contour(cn.df['Easting'], cn.df['Northing'], cn.df['MassArea'], cn.vent,
                 'Mass/Area', 'kg/m$^2$', save='./cerroNegro_contours.eps')

    filename2 = './data/pululagua_GITHUB.txt'
    pulu = Eruption(data=filename2, vent=Point(498000, 9630000), test=False)

    plot_grid(pulu.df, pulu.vent)
    pulu.df

    plot_contour(pulu.df['Easting'], pulu.df['Northing'], pulu.df['MassArea'], pulu.vent,
                 'Mass/Area', '', save='./pulu_contours.eps')
