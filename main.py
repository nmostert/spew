from tephra2io import tephra2_to_df
import visualisation as vis
import matplotlib.pyplot as plt
import sampling as smp


def grid_sampling(df):
    plt.figure()
    vis.plot_grid(smp.random_sample(df, frac=0.05))
    plt.show()

    plt.figure()
    vis.plot_grid(smp.resample(df, n=10))
    plt.show()


def k23a(phi):
    pass


def koyaguchi_test(df):
    pass


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'
    df = tephra2_to_df(filename)
    grid_sampling(df)
