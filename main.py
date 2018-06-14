from tephra2io import tephra2_to_df


filename = 'cerroNegro_regGrid_noWind_SOURCE.txt'

df = tephra2_to_df(filename)

print(df)
