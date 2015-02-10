#!/usr/bin/python

def prob1():

    import numpy as np
    from pandas import DataFrame
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import kurtosis

    data_files = ('ccd_response.tsv', 'doi10_filter_response.tsv',
                  'i_band_johnson_response.tsv')

    for data_file in data_files:
        df = pd.read_csv(data_file, sep='\t')

        if data_file == 'doi10_filter_response.tsv':
            df['i'] = df['i']**(df['air_mass'])

        i_trans_sub = df['i'][df['i'] > 0.01 * df['i']]
        #i_trans_sub = df['i']

        print('\nStats on ' + data_file)
        print('Mean = {0:.2f}'.format(np.mean(i_trans_sub)))
        print('Std = {0:.2f}'.format(np.std(i_trans_sub)))
        print('Kurtosis = {0:.2f}'.format(kurtosis(i_trans_sub)))

        #plt.plot(df['wavelength'], df['i'])
        #plt.show()

def prob2():

    import numpy as np
    from pandas import DataFrame
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import kurtosis

    print('\n\nProblem 2\n')

    data_files = ('ccd_response.tsv', 'doi10_filter_response.tsv',
                  'i_band_johnson_response.tsv')
    names = ('Interference filter', 'SDSS i-band', 'Johnson I')

    df_atmo = pd.read_csv('kpno_winter.tsv', sep='\t')
    df_atmo['i'][np.isnan(df_atmo['i'])] = 0.0

    airmasses = (1.0, 2.0)

    for airmass in airmasses:

        print('\nAirmass of ' + str(airmass) + '\n')

        for i, data_file in enumerate(data_files):
            df = pd.read_csv(data_file, sep='\t')

            if data_file == 'doi10_filter_response.tsv':
                df['i'] = df['i']**(df['air_mass'])

            df_atmo_interp = np.interp(df['wavelength'],
                                       df_atmo['wavelength'],
                                       df_atmo['i'])
            df['i'] *= df_atmo_interp

            # Correct for airmass
            df['i'] **= 1.0 / airmass

            i_trans_sub = df['i'][df['i'] > 0.01 * df['i']]
            #i_trans_sub = df['i']

            print('\nStats on ' + names[i])
            print('Mean = {0:.2f}'.format(np.mean(i_trans_sub)))
            print('Std = {0:.2f}'.format(np.std(i_trans_sub)))
            print('Kurtosis = {0:.2f}'.format(kurtosis(i_trans_sub)))

            if airmass == 1.0:
                plt.plot(df['wavelength'], df['i'],
                        label=names[i])
                plt.xlim([6500, 12000])
                plt.ylim([0, 1.1])
                plt.xlabel('Wavelength [Angstroms]')
                plt.ylabel('Transmission')
                plt.legend(loc='best')

        if airmass == 1.0:
            plt.savefig('fig_prob2.png')
            #plt.show()

def prob3():

    import numpy as np

    effic = 0.55 * 0.9 * 0.995

    m_V = 18

    flux = 10**((m_V - 0.03) / -2.5)

    area = 8 * 0.0729 # m^2
    neper_Vband = 0.16
    flux_photon = flux * 15.1*10**6 * neper_Vband * 2.9 * 1.9
    solid_angle = 2.9 * 1.9 * 3600.0**2
    S_gamma = flux_photon * area * solid_angle

    print('flux', flux)
    print('S_gamma', S_gamma)
    print('photon flux', flux_photon)

    npix = 6.4 * 10**7
    I_sky = S_gamma / area / npix
    RN = 10.0 # e / s
    t = (3.0 * RN)**2 / (I_sky * effic)

    print('I_sky =', I_sky)
    print('efficiency =', effic)
    print(t / 60.0, 'min')

    print('\n\n3b\n\n')

    t = 3600.0 # s
    bandpass = (3 * RN)**2 * npix / (effic * t * flux_photon * area * solid_angle)

    print('bandpass = ', bandpass)

def main():
    prob1()
    prob2()
    prob3()

if __name__ == '__main__':
    main()
