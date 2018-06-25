from tephra2io import tephra2_to_gdf
from shapely.geometry import Point
import gridutils as grd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    gdf = tephra2_to_gdf(filename)
    gdf.head()

    x_0, y_0 = 532290, 1382690

    # Find grid point closest to the vent (pseudo origin)
    x_p, y_p = grd.find_nearest_grid_point(gdf, (x_0, y_0))

    p_origin = gdf.loc[gdf.index[gdf['geometry'] == Point(x_p, y_p)]]

    p_origin

    # Construct grid  with indices of gdf to allow spatial iteration.
    ix_df = gdf.reset_index().pivot(
        index='Northing', columns='Easting')['index']

    # Find grid index of pseudo origin.
    y_idx = ix_df.index.get_loc(y_p)
    x_idx = ix_df.columns.get_loc(x_p)

    # construct a radial line from the pseudo origin to the edge (East)
    line = gdf.iloc[ix_df.iloc[y_idx, x_idx:len(ix_df.columns)]]

    def dist(row):
        return p_origin.geometry.distance(row['geometry'])

    radius = line.apply(lambda row: dist(row), axis=1)
    line = line.assign(r=radius)

    phi_labels = ['[-4,-3)',
                  '[-3,-2)',
                  '[-2,-1)',
                  '[-1,0)',
                  '[0,1)',
                  '[1,2)',
                  '[2,3)',
                  '[3,4)']

    log_S = []
    r_sqr = line['r'].values ** 2

    plt.xlabel("r^2")
    plt.ylabel("ln(S)")
    for phi in phi_labels:
        S = (line['Mass/Area'].values * (line[phi]))
        plt.plot(r_sqr, np.log(S), '.')
        plt.legend(phi_labels)
    plt.show()

    plt.xlabel("r^2")
    plt.ylabel("S")
    for phi in phi_labels:
        S = (line['Mass/Area'].values * (line[phi]))
        plt.plot(r_sqr, S, '.')
        plt.legend(phi_labels)
    plt.show()

    def func(x, m, c):
        return m * x + c

    for phi in phi_labels:
        log_S = np.log(line['Mass/Area'].values * (line[phi]))
        log_S = log_S.replace([np.inf, -np.inf], np.nan)
        nan_idx = [not na for na in log_S.isna()]
        clean_log_S = log_S.dropna()
        xx = r_sqr[nan_idx]
        params = curve_fit(func, xx, clean_log_S)
        print(params)
        plt.xlabel("r^2")
        plt.ylabel("ln(S)")
        plt.plot(xx, clean_log_S, '.')
        plt.plot(xx, func(
            xx, params[0][0], params[0][1]))
        lab = 'y = ' + '{:0.3e}'.format(params[0][0]) + \
            'x + ' + str(round(params[0][1], 4))
        plt.legend([phi, lab])
        plt.show()
