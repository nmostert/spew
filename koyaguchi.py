from tephra2io import tephra2_to_df
import visualisation as vis
import sampling as smp
import matplotlib.pyplot as plt


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'
    df = tephra2_to_df(filename)

    x_0, y_0 = 532290, 1382690
    x_p, y_p = smp.find_nearest_grid_point(df, (x_0, y_0))

    vis.plot_grid(df)
    plt.plot(x_0, y_0, 'r*')
    plt.plot(x_p, y_p, 'b*')
    plt.show()
