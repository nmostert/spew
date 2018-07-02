import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point


def tephra2_to_df(filename):
    # Reads input from Tephra2 into a GeoPandas GeoDataFrame
    # Geopandas is used to enable easy spatial calculations
    # (area, distance, etc.)
    df = pd.read_csv(filename, sep=" ", header=0)
    df = df.dropna(axis=1, how='all')
    df = df.fillna(0)
    df = df.rename(columns={'#Easting': 'Easting',
                            '[-4->-3)': '[-4,-3)',
                            '[-3->-2)': '[-3,-2)',
                            '[-2->-1)': '[-2,-1)',
                            '[-1->0)': '[-1,0)',
                            '[0->1)': '[0,1)',
                            '[1->2)': '[1,2)',
                            '[2->3)': '[2,3)',
                            '[3->4)': '[3,4)'})

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    df = GeoDataFrame(df.copy(), crs=crs, geometry=geometry)
    return df


if __name__ == "__main__":
    filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

    df = tephra2_to_df(filename)
    df.head()

    df['Mass/Area'] = df['Mass/Area'] * (1500 * 1500)
    df['Mass/Area'].sum()