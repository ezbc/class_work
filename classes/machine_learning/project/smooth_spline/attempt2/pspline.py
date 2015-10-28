#!/usr/bin/python

def main():

    from scipy.weave import inline
    import numpy.f2py as f2py

    with open('Pspline.f', 'r') as f:
        code = f.read()

    #inline(code)
    f2py.compile(code, modulename='pspline')

if __name__ == '__main__':
    main()
