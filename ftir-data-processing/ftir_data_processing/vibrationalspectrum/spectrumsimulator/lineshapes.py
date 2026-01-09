import numpy as np
from scipy.special import wofz

s2pi = np.sqrt(2*np.pi)
s2 = np.sqrt(2.0)


def gaussian(x, center, sigma, amplitude=1.0):
    """
    Return a 1-dimensional Gaussian function.
    gaussian(x, amplitude, center, sigma) =
        (amplitude/(s2pi*sigma)) * exp(-(1.0*x-center)**2 / (2*sigma**2))
    """
    return (amplitude/(s2pi*sigma)) * np.exp(-(1.0*x-center)**2 / (2*sigma**2))


def lorentzian(x, center, gamma, amplitude=1.0):
    """
    Return a 1-dimensional Lorentzian function.
    lorentzian(x, amplitude, center, sigma) =
        (amplitude/(1 + ((1.0*x-center)/sigma)**2)) / (pi*sigma)
    """
    return (amplitude/(1 + ((1.0*x-center)/gamma)**2)) / (np.pi*gamma)


def voigt(x, center, sigma, gamma, amplitude=1.0):
    """
    Return a 1-dimensional Voigt function.
    voigt(x, amplitude, center, sigma, gamma) =
        amplitude*wofz(z).real / (sigma*s2pi)
    see http://en.wikipedia.org/wiki/Voigt_profile
    """
    z = (x-center + 1j*gamma) / (sigma*s2)
    return amplitude * wofz(z).real / (sigma*s2pi)
