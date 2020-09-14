"""Correlated speed-dependent Nelkin-Ghatak profile.

The definitions of all the profiles available in the module use the
unitless formulas from [1]_. This means that to obtain unity, the
profiles should be integrated over variable :math:`u=\nu/\omega_D$`,
where :math:`\omega_D` is the Gaussian Doppler width of the profile.
Functions which after integration over frequencies expressed in Hz or in
cm-1 return unity, can be obtained by dividing the functions here by the
Doppler width expressed in appropriate units.  The frequency arguments
of the functions should be given in Hz.

.. [1] Ciuryło, Roman, Shapes of pressure- and Doppler-broadened
       spectral lines in the core and near wings, Physical Review A,
       58(2), 1029–1039 (1998).
       http://dx.doi.org/10.1103/PhysRevA.58.1029
"""
import numpy as np
import scipy.constants as C
from scipy.integrate import simps
from scipy.special import wofz

amu2kg = C.value('atomic mass constant') #: atomic mass unit to kilogram


def vp(x, x0, gam, dop):
    """The Voigt profile.

    Parameters
    ----------
    x : ndarray
        Frequency axis.
    x0 : float
        Line position with shift.
    gam : float
        Lorentzian half-width.
    dop : float
        Doppler half-width.
    """
    sigma = dop/np.sqrt(2*np.log(2))

    return wofz(((x - x0) + 1j*gam)/sigma/np.sqrt(2))/sigma/np.sqrt(2*np.pi)


def ngp(x, x0, gam, dop, nur, nui):
    """The Nelkin-Ghatak profile.

    Parameters
    ----------
    x : ndarray
        Frequency axis.
    x0 : float
        Line position with shift.
    gam : float
        Lorentzian half-width.
    dop : float
        Doppler half-width.
    nur : float
        Real part of frequency of optical collisions.
    nui : float
        Imaginary part of frequency of optical collisions.
    """
    nuopt = nur+1.0j*nui
    voigt = vp(x, x0, gam+nuopt, dop)

    return voigt/(1-np.pi*nuopt*voigt)


def sdvp(x, x0, v, vm, gamma, delta):
    """The speed-dependent Voigt profile.

    All arguments (except speeds) are in Hz.  Speeds are in m/s.

    Args:
    - x - calculated positions,
    - x0 - line position,
    - v - speeds at which other params are available,
    - vm - the most probable speed of the absorber,
    - gamma - speed-dependent pressure width HWHM,
    - delta - speed-dependent pressure shift,

    """
    norm = 2/np.sqrt(np.pi)**3
    vred = mirror(v/vm)
    vred[:v.size] *= -1
    sigma = x0/C.c*vm

    boltz = vred*np.exp(-vred**2)
    arg = (x[:, np.newaxis] - x0 - mirror(delta) + sigma*vred)/mirror(gamma)
    expr = boltz*(np.arctan(arg) + 1.0j/2*np.log(1 + arg**2))

    return norm*simps(expr, vred)


def csdngp_dop(x, x0, vreds, sigma, gamma, delta, nur, nui):
    """Calculate the CSDNGP profile.

    All arguments should be in the same frequency units.

    Parameters
    ----------
    x : ndarray
        x-axis
    x0 : float
        transition frequency
    vreds : ndarray
        reduced speeds (v/v_mps) at which parameters are available
    sigma : float
        doppler width of the line, not HWHM - Gaussian width (x0/C.c*vm).
    gamma : ndarray
        lorentzian halfwidth
    delta : ndarray
        pressure shift
    nur : ndarray
        real part of nu_opt
    nui : ndarray
        imaginary part of nu_opt

    Returns
    -------
    ndarray of complex
        see the module-level note for the interpretation of the result
    """
    norm = 2/np.sqrt(np.pi)**3
    vred = mirror(vreds)
    vred[:vreds.size] *= -1

    boltz = vred*np.exp(-vred**2)
    arg = (x[:, np.newaxis] - x0 - mirror(delta) - mirror(nui) + sigma*vred)/mirror(gamma + nur)
    Gexpr = boltz*(np.arctan(arg) + 1.0j/2*np.log(1 + arg**2))
    Hexpr = mirror(nur + 1.0j*nui)/sigma*Gexpr
    G = norm*simps(Gexpr, vred)
    H = norm*simps(Hexpr, vred)

    return G/(1 - np.pi*H)


def csdngp(x, x0, v, vm, gamma, delta, nur, nui):
    """Calculate the CSDNGP profile.

    All arguments (except speeds) are in Hz.  Speeds are in m/s.

    Args:
    - x - calculated positions,
    - x0 - line position,
    - v - speeds at which other params are available,
    - vm - the most probable speed of the absorber,
    - gamma - speed-dependent pressure width HWHM,
    - delta - speed-dependent pressure shift,
    - nur, nui - speed-depedent velocity of optical collisions.
    """
    norm = 2/np.sqrt(np.pi)**3
    vred = mirror(v/vm)
    vred[:v.size] *= -1
    sigma = x0/C.c*vm

    boltz = vred*np.exp(-vred**2)
    arg = (x[:, np.newaxis] - x0 - mirror(delta) - mirror(nui) + sigma*vred)/mirror(gamma + nur)
    Gexpr = boltz*(np.arctan(arg) + 1.0j/2*np.log(1 + arg**2))
    Hexpr = mirror(nur + 1.0j*nui)/sigma*Gexpr
    G = norm*simps(Gexpr, vred)
    H = norm*simps(Hexpr, vred)

    return G/(1 - np.pi*H)


def mirror(arr):
    """Return double-size array mirrored around zero.

    mirror([1, 2]) -> [2, 1, 1, 2]
    """
    if len(arr) == 1:
        return arr

    res = np.empty(2*arr.size, dtype=arr.dtype)
    res[:arr.size] = arr[::-1]
    res[arr.size:] = arr

    return res


def doppler_hwhm(nu0, T, m):
    """Calculate Doppler HWHM.

    Parameters
    ----------
    nu0 : float
        Transition frequency (in any unit proportional to Hz)
    T : float
        Temperature in Kelvin.
    m : float
        Absorber/emitter mass in g/mol.

    Returns
    -------
    hwhm : float
        Doppler half-width at half-maximum in the same unit as `nu0`.
    """
    return nu0*np.sqrt(C.k*T/amu2kg/m)/C.c*np.sqrt(2*np.log(2))
