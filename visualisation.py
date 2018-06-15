import matplotlib.pyplot as plt


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
