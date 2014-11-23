
def main():

    import numpy as np
    import matplotlib.pyplot as plt

    noise_scalings = (0.001, 0.06325, 0.1255, 0.18775, 0.25)
    resid_frees = (8.4, 105.3, 134.3, 144.9, 150.7)
    lams = (1.3e-5, 0.00079, 0.0033, 0.0075, 0.014)


    plt.plot(lams, noise_scalings)
    plt.show()

if __name__ == '__main__':
    main()
