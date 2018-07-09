import matplotlib.pyplot as plt
from eruption import Eruption
from shapely.geometry import Point
from matplotlib import rcParams
import gridutils as grd
rcParams.update({'figure.autolayout': True})


def plot_GS(df, phi_labels, point):
    data = df[df.geometry == grd.nearest_grid_point(df, point)]
    print(data.head())
    phi_vals = data[phi_labels].transpose()
    print(phi_vals.head())
    phi_vals = phi_vals / 100
    phi_vals.plot(kind='bar', rot=0, legend=False)


def plot_grid(df, vent=None):
    xx = df['Easting'].values
    yy = df['Northing'].values
    plt.axis('equal')
    ax = plt.plot(xx, yy, 'k.', ms=1)
    if vent is not None:
        ax.plot(vent.coords[0][0], vent.coords[0][1], 'r^', ms=3)
    ax.xlabel("Easting")
    ax.ylabel("Northing")
    return ax


if __name__ == "__main__":
    filename = './data/cerroNegro_regGrid_noWind_SOURCE.txt'

    vent = Point(532290, 1382690)

    cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)
    cn.df.head()

    plot_grid(cn.df, vent)

    plot_GS(cn.df, cn.phi_labels, vent)
