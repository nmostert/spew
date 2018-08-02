import tephra2io
import gridutils as grd
import os
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from scipy.spatial import Voronoi  # , voronoi_plot_2d
import shapely.geometry
import shapely.ops
from matplotlib.colors import LogNorm
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


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

    def sample(self, n, weight="Mass/Area", alpha=0.5, beta=0.1, filters=None,
               scorefun=None):
        """Sample the dataframe.

        Each point is assigned a uniform random number, which is weighted
        according to the score function. The n points with the highest assigned
        scores are returned.

        Keyword arguments:
        weight   -- the name of the column by which the sampling probability
                    should be weighted (default None)
        alpha    -- the sampling factor (default 1.0)
        filters  -- dictionary of the form {key:tuple} where key is the column
                    to be filtered by, and the tuple is the lower and upper
                    bounds of the accepted range. (default None)
                    Eg: {"radius": (1000, 500000)} will exclude all points with
                    radii not in the (1000,500000) range.
        scorefun -- user defined score function. Should be of the form
                    f(w, a, b), where w is the weight and a and b are
                    parameters.
                    Default w^(b*log(1+alpha)).

        """
        temp_df = self.df.copy()

        if filters is not None:
            for col in filters:
                temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                                  (temp_df[col] < filters[col][1])]
        if scorefun is None:
            def scorefun(w, alpha, beta):
                r = np.random.random_sample()
                return r * (w ** (beta * np.log(1 + alpha)))
        else:
            def scorefun(w, alpha, beta):
                r = np.random.random_sample()
                return r * f(w, alpha, beta)

        temp_df['norm_weight'] = temp_df[weight].values / \
            np.max(cn.df[weight].values)

        probs = temp_df.apply(lambda row: scorefun(
            row['norm_weight'], alpha, beta), axis=1)

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
    ax.plot(xx, yy, 'k.', ms=3)
    if vent is not None:
        ax.plot(vent.coords[0][0], vent.coords[0][1], 'r^', ms=6)
    if labels is None:
        plt.xlabel("Easting (m)")
        plt.ylabel("Northing (m)")
    else:
        plt.xlabel(labels[0])
        plt.ylabel(labels[0])
    return fig, ax


def plot_contour(xx, yy, zz, vent, title, cbar_label, cmap='', save=None):
    fig, ax = plt.subplots(1, 1)
    ax.axis('equal')
    plt.ylabel("Northing")
    plt.xlabel("Easting")
    plt.title(title)
    zz = zz + 0.000000000000000000001
    V = [0.0001, 0.000316, 0.001, 0.00316, 0.01, 0.0316,
         0.1, 0.316, 1, 3.16, 10, 31.6, 100, 316, 1000, 3162]
    plt.tricontour(xx, yy, zz, V,
                   norm=LogNorm(vmin=zz.min(), vmax=zz.max()))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(cbar_label, rotation=270)
    plt.plot(vent.x, vent.y, 'r^', ms=10)
    plt.gca().set_xlim(right=535000)
    plt.gca().set_ylim(bottom=1370000)
    if save is not None:
        plt.savefig(save, dpi=200, format='eps')


def get_voronoi(df):
    points = [
        row.geometry.coords[0]
        for i, row in df.iterrows()]

    vor = Voronoi(points)

    lines = [
        shapely.geometry.LineString(vor.vertices[line])
        for line in vor.ridge_vertices
        if -1 not in line]
    polys = shapely.ops.polygonize(lines)
    polys = [poly for poly in polys]
    df = pd.DataFrame()
    crs = {'init': 'epsg:4326'}
    vor_df = GeoDataFrame(df, crs=crs, geometry=polys)

    return vor_df


if __name__ == "__main__":
    filename = './data/cerroNegro_fineGrid_realWind_ONLINE.txt'

    cn = Eruption(data=filename, vent=Point(532290, 1382690), test=False)

    filters = {
        'radius': (5000, np.inf)}
    sample = cn.sample(50, weight='Mass/Area', alpha=0.1,
                       beta=0.1, filters=filters)

    fig, ax = plot_grid(sample, vent=cn.vent)
    fig, ax = plot_grid(cn.df, vent=cn.vent)
    temp_df = cn.df.copy()
    for col in filters:
        temp_df = temp_df[(temp_df[col] > filters[col][0]) &
                          (temp_df[col] < filters[col][1])]
    norm_mass = temp_df["Mass/Area"].values / \
        np.max(temp_df["Mass/Area"].values)

    def f1(w, a, b):
        return w ** a

    def f2(w, a, b):
        return np.exp(a * np.sqrt(w))

    def f3(w, a, b):
        return w ** (b * np.log(1 + a))

    def f4(w, a, b):
        return w ** (b * (a ** 2))

    xx = np.linspace(0, 1, num=200)

    a1s = np.linspace(0, 1, num=5)
    a2s = np.linspace(0, 1, num=5)
    a3s = np.linspace(0, 1, num=5)
    a4s = np.linspace(0, 1, num=5)

    yf1 = [f1(xx, a, 0) for a in a1s]
    yf2 = [f2(xx, a, 0) for a in a2s]
    yf3 = [f3(xx, a, 0.1) for a in a3s]
    yf4 = [f4(xx, a, 0.1) for a in a4s]

    for y1, a1 in zip(yf1, a1s):
        plt.plot(xx, y1)
        plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a1s])
        plt.title(r'$f_1(w) = w^{\alpha}$', fontsize=16)
        plt.xlabel(r'$w$ (Normalised)')
        plt.ylabel(r'$f_1(w)$')
        plt.savefig('f1.eps', dpi=200, format='eps')

    for y2, a2 in zip(yf2, a2s):
        plt.plot(xx, y2)
        plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a2s])
        plt.title(r'$f_2(w) = e^{\alpha \sqrt{w}}$', fontsize=16)
        plt.xlabel(r'$w$ (Normalised)')
        plt.ylabel(r'$f_2(w)$')
        plt.savefig('f2.eps', dpi=200, format='eps')

    for y3, a3 in zip(yf3, a3s):
        plt.plot(xx, y3)
        plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a2s])
        plt.title(r'$f_3(w) = w^{\beta \log(1+\alpha)}, \quad \beta=0.1$',
                  fontsize=16)
        plt.xlabel(r'$w$ (Normalised)')
        plt.ylabel(r'$f_3(w)$')
        plt.savefig('f3.eps', dpi=200, format='eps')

    for y4, a4 in zip(yf4, a4s):
        plt.plot(xx, y4)
        plt.legend([r'$\alpha = ' + str(a) + r'$' for a in a4s])
        plt.title(
            r'$f_4(w) = w^{\beta \alpha^2}, \quad \beta=0.1$', fontsize=16)
        plt.xlabel(r'$w$ (Normalised)')
        plt.ylabel(r'$f_4(w)$')
        plt.savefig('f4.eps', dpi=200, format='eps')

    tits = [
        r'$f_1(w) = w^{\alpha}$',
        r'$f_2(w) = e^{\alpha \sqrt{w}}$',
        r'$f_3(w) = w^{\beta \log(1+\alpha)}, \quad \beta=0.1$',
        r'$f_4(w) = w^{\beta \alpha^2}, \quad \beta=0.1$'
    ]
    ylabs = [
        r'$f_1(w)$',
        r'$f_2(w)$',
        r'$f_3(w)$',
        r'$f_4(w)$'
    ]
    count = 0
    count2 = 0
    for f, tit, yl in zip([f1, f2, f3, f4], tits, ylabs):
        count += 1
        count2 = 0
        for a, b in zip([0, 0.1, 0.5, 0.9, 1], [0.1] * 5):
            count2 += 1
            y = f(norm_mass, a, b)
            yr = [yy * np.random.random_sample() for yy in y]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3))
            sc = ax1.scatter(temp_df["Easting"], temp_df["Northing"],
                             marker=',', c=y, vmin=0, vmax=1)
            cbar = fig.colorbar(sc, ax=ax1)

            ax1.plot(cn.vent.x, cn.vent.y, 'r^', ms=8)
            samp = cn.sample(50, weight='Mass/Area', alpha=a,
                             beta=b, filters=filters, scorefun=f)
            ax1.scatter(samp["Easting"], samp["Northing"],
                        marker='.', color='r', s=5)
            ax2.hist(norm_mass, bins=100, log=True)
            ax3 = ax2.twinx()
            ax3.plot(xx, f(xx, a, b), 'r-')
            ax1.set(aspect='equal')
            ax1.set_title(tit, fontsize=14)
            ax1.set_ylabel(r'Northing')
            ax1.set_xlabel(r'Easting')
            ax2.set_title(r'$\alpha = $ ' + str(a), fontsize=14)
            ax2.set_xlabel(r'$w$ (Normalised)')
            ax2.set_ylabel(r'log(Frequency)')
            ax3.set_ylabel(yl)
            ax3.set_ylim([0, 1.1])
            plt.tight_layout()
            # plt.savefig('samp_f' + str(count) + '_a' + str(count2) +
            # '.eps', dpi=200, format='eps')
