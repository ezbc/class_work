#!/usr/bin/python

def main():

    import numpy.f2py as f2py

    with open('Pspline.f', 'r') as f:
        code = f.read()

    f2py.compile(code, modulename='pspline')

    import pspline

    help(pspline.pspline)


if __name__ == '__main__':
    main()
