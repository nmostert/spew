import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import Point
from tephra2io import tephra2_to_df
from eruption import plot_grid, plot_contour


def read_CN(filename):
    df = pd.read_csv(filename, header=0)
    df[["Northing", "Easting"]] = df[[
        "Northing", "Easting"]].fillna(method='ffill')

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    crs = {'init': 'epsg:4326'}
    df = GeoDataFrame(df.copy(), crs=crs, geometry=geometry)
    return df


filename = './data/cerroNegro_regGrid_noWind_SOURCE.txt'

df = tephra2_to_df(filename)
df.head()

filename = './data/mass_sanitized.csv'

df2 = read_CN(filename)

df2.head()
vent = Point(532290, 1382690)

plot_contour(df2['Easting'], df2['Northing'], df2['Thickness'], vent,
             'Thickness', "m", save='./cn_thickness.eps')
plot_contour(df2['Easting'], df2['Northing'], df2['Density'], vent,
             'Density', "kg/m^3", save='./cn_density.eps')

plot_contour(df2['Easting'], df2['Northing'], df2['Isomass'], vent,
             'Isomass', "kg/m^2", save='./cn_isomass.eps')

dfA = df2[df2["Layer"] == "A"]
dfB = df2[df2["Layer"] == "B"]
len(dfA)
len(dfB)

plot_contour(dfA['Easting'], dfA['Northing'], dfA['Thickness'], vent,
             'Layer A Thickness', "m", save='./cnA_thickness.eps')

plot_contour(dfA['Easting'], dfA['Northing'], dfA['Density'],
             vent, 'Layer A Density', "kg/m^3", save='./cnA_density.eps')

plot_contour(dfA['Easting'], dfA['Northing'], dfA['Isomass'],
             vent, 'Layer A Isomass', "kg/m^2", save='./cnA_isomass.eps')

plot_contour(dfB['Easting'], dfB['Northing'], dfB['Thickness'],
             vent, 'Layer B Thickness', "m", save='./cnB_thickness.eps')

plot_contour(dfB['Easting'], dfB['Northing'], dfB['Density'],
             vent, 'Layer B Density', "kg/m^3", save='./cnB_density.eps')

plot_contour(dfB['Easting'], dfB['Northing'], dfB['Isomass'],
             vent, 'Layer B Isomass', "kg/m^2", save='./cnB_isomass.eps')

plot_grid(df, vent=vent)
plot_grid(df2, vent=vent)