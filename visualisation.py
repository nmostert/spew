import matplotlib.pyplot as plt
from tephra2io import tephra2_to_gdf
import sampling as smp
from shapely.geometry import Point
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


def plot_GS(df, easting, northing):
    df1 = df.loc[lambda df: df.Easting == easting]
    point = df1.loc[lambda df1: df1.Northing == northing]
    phi_vals = point[point.columns[4:]].transpose()
    phi_vals = phi_vals / 100
    phi_vals.plot(kind='bar', rot=0, legend=False)


def plot_grid(df):
    xx = df['Easting'].values
    yy = df['Northing'].values
    plt.plot(xx, yy, 'ko', ms=1)
    plt.xlabel("Easting")
    plt.ylabel("Northing")


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    gdf = tephra2_to_gdf(filename)
    gdf.head()

    x_0, y_0 = 532290, 1382690  # Actual vent location

    # Find grid coordinates closest to the vent
    x_p, y_p = smp.find_nearest_grid_point(gdf, (x_0, y_0))

    # Locate data for that grid point
    gdf.loc[gdf.index[gdf['geometry'] == Point(x_p, y_p)]]
