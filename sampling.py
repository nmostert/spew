from tephra2io import tephra2_to_gdf
# from shapely.geometry import Point
import matplotlib.pyplot as plt
import visualisation as vis


def random_sample(df, frac=0.1):
    # Takes a uniform random sample from a dataset
    # frac is the proportion of the total number of data points to be sampled.
    return df.sample(frac=0.05, replace=False)


def regular_sample(df, n=5):
    # Resamples a dataframe by taking every 5 samples in either direction
    # Meant to be used on a regular grid, but should work with irregular data
    resample = get_index_grid(gdf).iloc[::n, ::n]
    return df.iloc[resample.values.flatten(), :]


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


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    gdf = tephra2_to_gdf(filename)
    gdf.head()

    gdf1 = regular_sample(gdf, 5)
    gdf2 = random_sample(gdf, 0.1)

    vis.plot_grid(gdf)
    plt.show()

    plt.figure()
    vis.plot_grid(gdf1)
    plt.show()

    plt.figure()
    vis.plot_grid(gdf2)
    plt.show()
