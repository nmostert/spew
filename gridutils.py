from tephra2io import tephra2_to_gdf
from shapely.geometry import Point
import matplotlib.pyplot as plt
import visualisation as vis
import numpy as np
import pandas as pd
from geopandas import GeoDataFrame


def random_sample(df, frac=0.1):
    # Takes a uniform random sample from a dataset
    # frac is the proportion of the total number of data points to be sampled.
    return df.sample(frac=0.05, replace=False)


def regular_sample(gdf, n=5):
    # Resamples a dataframe by taking every 5 samples in either direction
    # Meant to be used on a regular grid, but should work with irregular data
    resample = get_index_grid(gdf).iloc[::n, ::n]
    return gdf.iloc[resample.values.flatten(), :]


def get_index_grid(df):
    # Creates dataframe that is the same shape as the spatial grid
    # Makes spatial iteration easier
    return df.reset_index().pivot(
        index='Northing', columns='Easting')['index']


def find_nearest_grid_point(df, point):
    # Finds existing data point in df that is closest to given point
    # TODO: rewrite using GeoPandas and shapely.geometry.distance()
    x_p = df.ix[(df['Easting'] - point[0]).abs().argsort()[0]]['Easting']
    y_p = df.ix[(df['Northing'] - point[1]).abs().argsort()[0]]['Northing']
    return (x_p, y_p)


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

    print(x_east)
    print(x_west)
    print(y_north)
    print(y_south)
    xx, yy = np.meshgrid(x, y)

    cols = ['Easting', 'Northing', 'Elev.']

    df = pd.DataFrame(columns=cols)
    for i, c in enumerate(xx):
        for j, p in enumerate(c):
            df = df.append(
                {'Easting': int(p), 'Northing': int(yy[i][j]), 'Elev.': int(elevation)}, ignore_index=True)

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    gdf = GeoDataFrame(df, crs=crs, geometry=geometry)
    return gdf


def write_grid_file(df, filename):
    df.to_csv(filename, sep=' ', columns=['Easting', 'Northing', 'Elev.'],
              index=False, header=False)


if __name__ == "__main__":
    vent = Point(532290, 1382690)
    grid = construct_grid(vent, 50, 50, 50, 50, 50, 1500)

    filename = 'cerroNegro_radGrid_noWind_ONLINE.txt'
    gdf = tephra2_to_gdf(filename)

    vis.plot_grid(grid, vent)
    vis.plot_grid(gdf, vent)
    plt.gca().get_lines()[0].set_color("blue")
    plt.gca().get_lines()[0].set_color("red")

    plt.show()

    write_grid_file(grid, "text.pos")
