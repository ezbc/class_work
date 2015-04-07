import os
import numpy as np
import pandas as pd
import pickle

dir_files = os.listdir('./')

zip_files = []

ids = ['08839', '08840', '08841', '08842']
Zs = [0.0001, 0.001, 0.01, 0.03]

df = pd.DataFrame(columns=[str(Z) for Z in Zs])
df = {}
for i, Z in enumerate(Zs):
    df[str(Z)] = {}
    df[str(Z)]['id'] = ids[i]

if 0:
    df['Mass [Msun]'] = np.zeros(20)
    df['Luminosity [Lsun]'] = np.zeros(20)
    df['Teff [K]'] = np.zeros(20)

count = 0
masses = []

#for dir_file in dir_files:
for i, ezid in enumerate(ids):
    os.system('rm -rf summary.txt')

    os.system('unzip ezweb_' + ezid + '.zip')

    data = np.loadtxt('summary.txt')

    if 0:
        df['mass [msun]'][count] = data[2]
        df['luminosity [lsun]'][count]= 10**data[3]
        df['teff [k]'][count] = 10**data[5]

        count += 1

    #nested_df = pd.DataFrame()
    nested_df = {}
    nested_df['Time [years]'] = data[:, 1]
    nested_df['Mass [Msun]'] = data[:, 2]
    nested_df['Luminosity [Lsun]'] = 10**data[:, 3]
    nested_df['Radius [Rsun]'] = 10**data[:, 4]
    nested_df['Teff [K]'] = 10**data[:, 5]
    nested_df['Tc [K]'] = 10**data[:, 6]
    nested_df['Degeneracy'] = data[:, 9]
    nested_df['Central Hydrogen Mass Fraction'] = data[:, 10]
    nested_df['L_PP [Lsun]'] = data[:, 18]
    nested_df['L_CNO [Lsun]'] = data[:, 19]
    nested_df['L_3a [Lsun]'] = data[:, 20]
    nested_df['L_Z [Lsun]'] = data[:, 21]
    nested_df['L_V [Lsun]'] = data[:, 22]

    if 0:
        df[str(Zs[i])] = ('Mass [Msun]', 'Luminosity [Lsun]', 'Teff [K]',)

    df[str(Zs[i])] = nested_df

    #masses.append(data[2])

#df = df.sort(df.columns[0])

#df.to_csv('../evol_data.csv', index=False)
#df.save('../evol_data.pkl')

with open('../evol_data.pkl', 'wb') as handle:
    pickle.dump(df, handle)



'''
0      1               2               3               4               5
index  t               mass            Lsun            Radius          Ts
00000  0.00000000E+00  1.00000000E+00 -1.54655885E-01 -5.27369747E-02  3.74939401E+00  7.12628079E+00  1.89376232E+00  1.71504221E+01 -1.72462659E+00  6.98287816E-01  2.81392784E-01  2.32069350E-03  2.48802121E-03  1.00995897E-02  1.87353434E+03  1.27401361E+07  1.48570470E+10  6.51226079E-01  5.22548254E-02  2.63533706E-34  0.00000000E+00  1.62294544E-08  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
'''


