#!/usr/bin/python

def main():

    from scipy.weave import inline

    with open('pspline.c', 'r') as f:
        code = f.read()

    inline(code)

if __name__ == '__main__':
    main()
