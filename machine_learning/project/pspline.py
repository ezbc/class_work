#!/usr/bin/python




def main():

    import pickle

    test_data = pickle.load(open('HT2003_data_test100.pickle'))

    for key in test_data: print key


if __name__ == '__main__':
    main()



