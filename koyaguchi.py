from tephra2io import tephra2_to_df
# import visualisation as vis
# import sampling as smp
# import matplotlib.pyplot as plt


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    gdf = tephra2_to_df(filename)
    gdf.head()
