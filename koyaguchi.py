from tephra2io import tephra2_to_gdf
from shapely.geometry import Point
import sampling as smp
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    gdf = tephra2_to_gdf(filename)
    gdf.head()

    x_0, y_0 = 532290, 1382690

    # Find grid point closest to the vent (pseudo origin)
    x_p, y_p = smp.find_nearest_grid_point(gdf, (x_0, y_0))

    p_origin = gdf.loc[gdf.index[gdf['geometry'] == Point(x_p, y_p)]]

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
