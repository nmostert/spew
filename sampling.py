
def random_sample(df, frac=0.1):
    return df.sample(frac=0.05, replace=False)


def resample(df, n=5):
    return df.iloc[::n, :]


def find_nearest_grid_point(df, point):
    x_p = df.ix[(df['Easting'] - point[0]).abs().argsort()[0]]['Easting']
    y_p = df.ix[(df['Northing'] - point[1]).abs().argsort()[0]]['Northing']
    return (x_p, y_p)
