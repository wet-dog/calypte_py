from __future__ import division, print_function, absolute_import

import math
import tmm
from typing import NamedTuple

class Parameters(NamedTuple):
    f_k: int
    n_k: int
    f_m: int
    n_m: int
    k_m: int
    f_a: int
    n_a: int

def cauchy(A, B, wavelength):
    return A + B*wavelength**-2

def melanin_imaginary(wavelength):
    a_m = 0.56
    b_m = 270
    return a_m * math.exp(-wavelength / b_m)

def keratin_index(wavelength):
    A_k = 1.532
    B_k = 5890
    return cauchy(A_k, B_k, wavelength)

def melanin_index(wavelength):
    A_m = 1.648
    B_m = 23700
    return cauchy(A_m, B_m, wavelength)

def calc_top_melanin_layer(wavelength):
    f_k = 0.05
    n_k = keratin_index(wavelength)
    f_m = 0.4
    n_m = melanin_index(wavelength)
    k_m = melanin_imaginary(wavelength)
    f_a = 0.55
    n_a = 1
    return Parameters(f_k, n_k, f_m, n_m, k_m, f_a, n_a)

def calc_keratin_layer(wavelength):
    f_k = 1
    n_k = keratin_index(wavelength)
    f_m = 0
    n_m = 0
    k_m = 0
    f_a = 0
    n_a = 0
    return Parameters(f_k, n_k, f_m, n_m, k_m, f_a, n_a)

def calc__layer(wavelength):
    f_k = 1
    n_k = keratin_index(wavelength)
    f_m = 0
    n_m = 0
    k_m = 0
    f_a = 0
    n_a = 0
    return Parameters(f_k, n_k, f_m, n_m, k_m, f_a, n_a)

def calc_refractive_index(params: Parameters):
    n_eff = (params.f_k*params.n_k
            + params.f_m*(params.n_m - complex(params.k_m))
            + params.f_a*params.n_a)

    return n_eff

if __name__ == "__main__":
    f_m = 0
    wavelength = 500
    keratin_layer_index = calc_refractive_index(calc_keratin_layer(wavelength))
    # Need to calculate melanosome refractive index and have a layer separate
    # to it from the internal melanosome layer
    print("keratin refractive index", keratin_layer_index)
