from tephra2io import tephra2_to_df
# import visualisation as vis
import sampling as smp
# import matplotlib.pyplot as plt


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'
    df = tephra2_to_df(filename)

    x_0, y_0 = 532290, 1382690
    x_p, y_p = smp.find_nearest_grid_point(df, (x_0, y_0))

    xx = list(set(df['Easting']))
    xx.sort()

    r_x = [i for i in xx if i >= x_p]
    r_y = [y_p] * len(r_x)

    df1 = df.loc[lambda df: df.Easting in r_x]
    print(df1)
