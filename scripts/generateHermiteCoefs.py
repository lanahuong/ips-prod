#!/usr/bin/env python3
import sys

import numpy


def format_matrix(poly):
    return 'arma::mat{{' + ','.join('{' + ','.join(map(str, p)) + '}' for p in poly) + '}}'


def get_coefs(i):
    return numpy.polynomial.hermite.hermgauss(i) if i != 0 else [[0], [0]]


pre_dir1 = '#ifndef SOLVER_HERM_COEFS'
pre_dir2 = '#define SOLVER_HERM_COEFS'
pre_dir3 = '#define HERM_QUADRA_N_MAX '
pre_dir4 = '#define HERMITE_QUADRA_COEFS '
pre_dir5 = '#endif'

if __name__ == "__main__":
    if len(sys.argv) == 2:
        n = int(sys.argv[1])
        data = '{' + ','.join(format_matrix(get_coefs(i)) for i in range(0, n + 1)) + '}'
        print(pre_dir1)
        print(pre_dir2)
        print(pre_dir3 + str(n))
        print(pre_dir4 + data)
        print(pre_dir5)
