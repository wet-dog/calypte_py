from __future__ import division, print_function, absolute_import

import tmm
import math
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from typing import NamedTuple

class Parameters(NamedTuple):
    f_k: int
    n_k: int
    f_m: int
    n_m: int
    k_m: int
    f_a: int
    n_a: int

class Melanosome(NamedTuple):
    top_membrane: complex
    internal_layer: complex
    bottom_membrane: complex

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

def calc_melanosome_internal_layer_index(wavelength):
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

def calc_melanin_layer(wavelength):
    f_k = 0
    n_k = 0
    f_m = 1
    n_m = melanin_index(wavelength)
    k_m = melanin_imaginary(wavelength)
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
    melanosome_membrane_index = calc_refractive_index(calc_melanin_layer(wavelength))
    melanosome_internal_index = calc_refractive_index(calc_melanosome_internal_layer_index(wavelength))

    # Thicknesses in nm
    keratin_cortex_thickness = 5
    keratin_separation_thickness = 50
    melanosome_membrane_thickness = 30
    melanosome_internal_thickness = 50
    top_melanosome_internal_thickness = 100

    melanosome_layers = 12

    melanosome = Melanosome(melanosome_membrane_index, melanosome_internal_index, melanosome_membrane_index)

    print("keratin refractive index", keratin_layer_index)
    print("melanosome membrane refractive index", melanosome_membrane_index)
    print("melanosome internal layer index", melanosome_internal_index)

    #index of refraction of my material: wavelength in nm versus index.
    material_nk_data = np.array([[200, 2.1+0.1j],
                            [300, 2.4+0.3j],
                            [400, 2.3+0.4j],
                            [500, 2.2+0.4j],
                            [750, 2.2+0.5j]])
    material_nk_fn = interpolate.interp1d(material_nk_data[:,0].real,
                            material_nk_data[:,1], kind='quadratic')
    d_list = [np.inf,300,np.inf] #in nm
    lambda_list = np.linspace(200,750,400) #in nm
    T_list = []
    for lambda_vac in lambda_list:
        n_list = [1, material_nk_fn(lambda_vac), 1]
        T_list.append(tmm.coh_tmm('s',n_list,d_list,0,lambda_vac)['T'])
    plt.figure()
    plt.plot(lambda_list,T_list)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fraction of power transmitted')
    plt.title('Transmission at normal incidence')
