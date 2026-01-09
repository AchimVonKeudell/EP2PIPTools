import numpy as np
from . import weights

h = 6.62607004e-34  # Planck's constant
c = 299792458  # speed of light in m s-1
N_A = 6.02214129e23  # avogadro's number in mol-1
k_B = 1.38064852e-23  # Boltzmann constant in m2 kg s-2 K-1
c2 = h * 100 * c / k_B  # second radiative constant in cm K
line_strength_constant = 8 * np.pi * 100 * c  # in cm/s

CO_dunham_matrices = [
    np.array([
        [0, 1.9312808724, -6.121468e-6, 5.8272e-12],
        [2169.8135802, -1.75044121e-2, 1.1526e-9, 1.7375e-13],
        [-13.2883076, 5.487e-7, -1.8050e-10, 0],
        [1.051127e-2, 2.541e-8, 0, 0],
        [5.7440e-5, 0, 0, 0],
        [9.831e-7, 0, 0, 0],
        [-3.1660e-8, 0, 0, 0]
    ]),
    np.array([
        [0, 1.846151573, -5.593e-6, 5.0131e-12],
        [2121.439929, -1.635955118e-2, 0.9637e-9, -1.4379e-13],
        [-12.70252823, 4.8339e-7, 1.5103e-10, 0],
        [0.983714e-2, 2.1862e-8, 0, 0],
        [5.1647e-5, 0, 0, 0],
        [8.995e-7, 0, 0, 0],
        [-2.7838e-8, 0, 0, 0]
    ]),
    np.array([
        [0, 1.839113772, -5.5505e-6, 4.956e-12],
        [2117.399036, -1.626607176e-2, 0.9546e-9, -1.4188e-13],
        [-12.65420846, 4.797e-7, -1.493e-10, 0],
        [0.978093e-2, 2.165e-8, 0, 0],
        [5.1254e-5, 0, 0, 0],
        [8.909e-7, 0, 0, 0],
        [-2.752e-8, 0, 0, 0]
    ])
]
CO_statistical_weights = [
    weights.FixedRotationalStatisticalWeight(1, 1),
    weights.FixedRotationalStatisticalWeight(2, 1),
    weights.FixedRotationalStatisticalWeight(1, 1),
    weights.FixedRotationalStatisticalWeight(6, 1),
    weights.FixedRotationalStatisticalWeight(2, 1),
    weights.FixedRotationalStatisticalWeight(12, 1),
]

# vibrational constants for isotopologues of CO2
CO2_626_vibrational_constants = [
    np.zeros((3)),  # first order vibrational constants
    np.zeros((3, 3)),  # second order vibrational constants
    np.zeros((3, 3, 3)),  # third order vibrational constants
    0  # vibrational constant for angular momentum
]
CO2_636_vibrational_constants = [
    np.zeros((3)),  # first order vibrational constants
    np.zeros((3, 3)),  # second order vibrational constants
    np.zeros((3, 3, 3)),  # third order vibrational constants
    0  # vibrational constant for angular momentum
]
CO2_628_vibrational_constants = [
    np.zeros((3)),  # first order vibrational constants
    np.zeros((3, 3)),  # second order vibrational constants
    np.zeros((3, 3, 3)),  # third order vibrational constants
    0  # vibrational constant for angular momentum
]
# Source: Suzuki, I. (1968). General anharmonic force constants of carbon dioxide. Journal of Molecular Spectroscopy, 25(4), 479–500. https://doi.org/10.1016/S0022-2852(68)80018-9
CO2_626_vibrational_constants[0][0] = 1354.31
CO2_626_vibrational_constants[0][1] = 672.85
CO2_626_vibrational_constants[0][2] = 2396.32
CO2_626_vibrational_constants[1][0, 0] = -2.93
CO2_626_vibrational_constants[1][0, 1] = -4.61
CO2_626_vibrational_constants[1][0, 2] = -19.82
CO2_626_vibrational_constants[1][1, 1] = 1.35
CO2_626_vibrational_constants[1][1, 2] = -12.31
CO2_626_vibrational_constants[1][2, 2] = -12.47
CO2_626_vibrational_constants[2] = -0.97

# Source: Suzuki, I. (1968). General anharmonic force constants of carbon dioxide. Journal of Molecular Spectroscopy, 25(4), 479–500. https://doi.org/10.1016/S0022-2852(68)80018-9
CO2_636_vibrational_constants[0][0] = 1354.31
CO2_636_vibrational_constants[0][1] = 653.70
CO2_636_vibrational_constants[0][2] = 2328.12
CO2_636_vibrational_constants[1][0, 0] = -2.93
CO2_636_vibrational_constants[1][0, 1] = -4.46
CO2_636_vibrational_constants[1][0, 2] = -19.33
CO2_636_vibrational_constants[1][1, 1] = 1.28
CO2_636_vibrational_constants[1][1, 2] = -11.54
CO2_636_vibrational_constants[1][2, 2] = -11.71
CO2_636_vibrational_constants[2] = -0.91

# Source: Suzuki, I. (1968). General anharmonic force constants of carbon dioxide. Journal of Molecular Spectroscopy, 25(4), 479–500. https://doi.org/10.1016/S0022-2852(68)80018-9
CO2_628_vibrational_constants[0][0] = 1315.21
CO2_628_vibrational_constants[0][1] = 667.72
CO2_628_vibrational_constants[0][2] = 2378.53
CO2_628_vibrational_constants[1][0, 0] = -2.76
CO2_628_vibrational_constants[1][0, 1] = -4.45
CO2_628_vibrational_constants[1][0, 2] = -19.04
CO2_628_vibrational_constants[1][1, 1] = 1.34
CO2_628_vibrational_constants[1][1, 2] = -12.18
CO2_628_vibrational_constants[1][2, 2] = -12.34
CO2_628_vibrational_constants[2] = -0.95

# the resonance constants for the isotopologues of CO2
# Source: Stull, V. R., Wyatt, P. J., & Plass, G. N. (1962). Vibrational Energies of the CO 2 Molecule. The Journal of Chemical Physics, 37(7), 1442–1445. https://doi.org/10.1063/1.1733302
CO2_626_resonance_constants = [51.31, 0.15, 0.41, 0.78, 0]
CO2_636_resonance_constants = [46.52, 0.10, 0.37, 0.74, 0]
CO2_628_resonance_constants = [52.5, 0.10, 0.40, 0.75, 0]

# rotational constants for isotopologues of CO2
CO2_626_rotational_constants = [
    0,  # rotational B constant
    np.zeros((3)),  # first order vibrational correction factors
    np.zeros((3, 3)),  # second order vibrational correction factors
    0,  # centrifugal constant
    0  # rotational H constant
]
CO2_636_rotational_constants = [
    0,  # rotational B constant
    np.zeros((3)),  # first order vibrational correction factors
    np.zeros((3, 3)),  # second order vibrational correction factors
    0,  # centrifugal constant
    0  # rotational H constant
]
CO2_628_rotational_constants = [
    0,  # rotational B constant
    np.zeros((3)),  # first order vibrational correction factors
    np.zeros((3, 3)),  # second order vibrational correction factors
    0,  # centrifugal constant
    0  # rotational H constant
]

# Source: Courtoy, C. P. (1956). Spectres de vibration-rotation de molecules simples diatomiques ou polyatomiques avec long parcours d’absorption.
CO2_626_rotational_constants[0] = 0.39022
CO2_626_rotational_constants[1][0] = 12.40e-4
CO2_626_rotational_constants[1][1] = -7.4e-4
CO2_626_rotational_constants[1][2] = 30.60e-4
CO2_626_rotational_constants[2] = 0.1333e-6
CO2_626_rotational_constants[3] = 0.009e-12

# Source: Klarenaar, B. L. M., Engeln, R., Bekerom, D. C. M. Van Den, Sanden, M. C. M. Van De, & Guaitella, O. (2017). Time evolution of vibrational temperatures in a CO 2 glow discharge measured with infrared absorption spectroscopy Experimental methods, 1–12. https://doi.org/10.1088/1361-6595/aa902e
# Change of rotational spacing by vibrational quanta:
# Suzuki, I. (1968). General anharmonic force constants of carbon dioxide. Journal of Molecular Spectroscopy, 25(4), 479–500. https://doi.org/10.1016/S0022-2852(68)80018-9
CO2_636_rotational_constants[0] = 0.39024
CO2_636_rotational_constants[1][0] = 12.41e-4
CO2_636_rotational_constants[1][1] = -7.01e-4
CO2_636_rotational_constants[1][2] = 29.48e-4
CO2_636_rotational_constants[2] = 0.1332e-6
CO2_636_rotational_constants[3] = 0.009e-12

# Source: Klarenaar, B. L. M., Engeln, R., Bekerom, D. C. M. Van Den, Sanden, M. C. M. Van De, & Guaitella, O. (2017). Time evolution of vibrational temperatures in a CO 2 glow discharge measured with infrared absorption spectroscopy Experimental methods, 1–12. https://doi.org/10.1088/1361-6595/aa902e
# Change of rotational spacing by vibrational quanta:
# Suzuki, I. (1968). General anharmonic force constants of carbon dioxide. Journal of Molecular Spectroscopy, 25(4), 479–500. https://doi.org/10.1016/S0022-2852(68)80018-9
CO2_628_rotational_constants[0] = 0.36819
CO2_628_rotational_constants[1][0] = 11.38e-4
CO2_628_rotational_constants[1][1] = -7.03e-4
CO2_628_rotational_constants[1][2] = 28.79e-4
CO2_628_rotational_constants[2] = 0.1187e-6
CO2_628_rotational_constants[3] = 0.0075e-12 

CO2_statistical_weights = [
    weights.CO2RotationalStatisticalWeight(1),
    weights.CO2RotationalStatisticalWeight(2),
    weights.FixedRotationalStatisticalWeight(1, 1),
    weights.FixedRotationalStatisticalWeight(6, 1),
    weights.FixedRotationalStatisticalWeight(2, 1),
    weights.FixedRotationalStatisticalWeight(12, 1),
    weights.CO2RotationalStatisticalWeight(1),
    weights.FixedRotationalStatisticalWeight(6, 1),
]

# some constants as defined in McDowell, R. S. (1988). Rotational partition functions for linear molecules. The Journal of Chemical Physics, 88(1), 356–361. https://doi.org/10.1063/1.454608
# to calculate the rotational partition sum
CO2_molecular_properties = [
    [
        2,  # classical symmetry number sigma
        1,  # kappa
        1  # nuclear spin degeneracy I
    ],
    [
        1,  # classical symmetry number
        1,  # kappa
        1  # nuclear spin degeneracy
    ],
    [
        1,  # classical symmetry number
        1,  # kappa
        1  # nuclear spin degeneracy
    ]
]

CO2_vibrational_constants = [
    CO2_626_vibrational_constants,
    CO2_636_vibrational_constants,
    CO2_628_vibrational_constants
]

CO2_resonance_constants = [
    CO2_626_resonance_constants,
    CO2_636_resonance_constants,
    CO2_628_resonance_constants
]

CO2_rotational_constants = [
    CO2_626_rotational_constants,
    CO2_636_rotational_constants,
    CO2_628_rotational_constants
]


NH3_statistical_weights = [
    weights.NH3RotationalStatisticalWeight(6, 1),
    weights.FixedRotationalStatisticalWeight(4, 1)
]


# fit the ExoMol data of the vibrational ground level (sigma=6), constants are in accordance with:
#   Dowling, Journal of Molecular Spectrscopy 27, 527-538 (1968), 10.1016/0022-2852(68)90058-1
# B = 9.94 cm^-1
NH3_rotational_constants = [
    np.array([[0.00000000e+00, -3.67019269e+00, -9.20868377e-04],
              [9.92368290e+00, 1.23156886e-03, 6.33613450e-07],
              [-7.00411058e-04, -1.14646013e-07, -3.05081580e-10]])
]

NH3_vibrational_constants = [
    [[3423.5, -67.002], {'s': [0, 758.9, 45.608], 'a': [0.793403, 857.780, 44.350]}, [3272.5, -36.113], 1608.7]
]