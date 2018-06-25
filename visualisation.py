import matplotlib.pyplot as plt
from tephra2io import tephra2_to_gdf
import gridutils as grd
from shapely.geometry import Point
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


def plot_GS(df, easting, northing):
    df1 = df.loc[lambda df: df.Easting == easting]
    point = df1.loc[lambda df1: df1.Northing == northing]
    phi_vals = point[point.columns[4:]].transpose()
    phi_vals = phi_vals / 100
    phi_vals.plot(kind='bar', rot=0, legend=False)


def plot_grid(df, vent=None):
    xx = df['Easting'].values
    yy = df['Northing'].values
    plt.axis('equal')
    plt.plot(xx, yy, 'k.', ms=1)
    if vent is not None:
        plt.plot(vent.coords[0][0], vent.coords[0][1], 'r*', ms=3)
    plt.xlabel("Easting")
    plt.ylabel("Northing")


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    vent = Point(532290, 1382690)

    gdf = tephra2_to_gdf(filename)
    gdf.head()

    plot_grid(gdf, vent)
