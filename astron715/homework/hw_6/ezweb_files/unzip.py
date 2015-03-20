import os
import numpy as np
import pandas as pd

dir_files = os.listdir('./')

zip_files = []

df = pd.DataFrame(columns=('Mass [Msun]', 'Luminosity [Lsun]', 'Teff [K]',))
df['Mass [Msun]'] = np.zeros(20)
df['Luminosity [Lsun]'] = np.zeros(20)
df['Teff [K]'] = np.zeros(20)

count = 0
masses = []

for dir_file in dir_files:
    os.system('rm -rf summary.txt')
    if '.zip' in dir_file:
        os.system('unzip ' + dir_file)

        data = np.loadtxt('summary.txt')[0]


        if data[2] not in masses:
            df['Mass [Msun]'][count] = data[2]
            df['Luminosity [Lsun]'][count]= 10**data[3]
            df['Teff [K]'][count] = 10**data[5]

            count += 1

        masses.append(data[2])

df = df.sort(df.columns[0])

df.to_csv('../cluster_data.csv', index=False)




'''
0      1               2               3               4               5
index  t               mass            Lsun            Radius          Ts
00000  0.00000000E+00  1.00000000E+00 -1.54655885E-01 -5.27369747E-02  3.74939401E+00  7.12628079E+00  1.89376232E+00  1.71504221E+01 -1.72462659E+00  6.98287816E-01  2.81392784E-01  2.32069350E-03  2.48802121E-03  1.00995897E-02  1.87353434E+03  1.27401361E+07  1.48570470E+10  6.51226079E-01  5.22548254E-02  2.63533706E-34  0.00000000E+00  1.62294544E-08  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
'''


