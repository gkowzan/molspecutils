"""Functions related to optical cavities."""
from typing import List, Set, Dict, Tuple, Optional, Callable, Union
import numpy as np
import scipy.constants as C
from shed.units import wn2lambda
np_float = Union[float, np.ndarray]

#######################################################################
# Cavity-enhanced spectroscopies                                      #
#######################################################################
def reflected_field(omg, r, fsr):
    """Reflected field coefficient. E_ref/E_in."""
    return r*(np.exp(1.0j*omg/fsr)-1)/(1-r**2*np.exp(1.0j*omg/fsr))


def pdh_error_signal(omg, r, fsr):
    """For high frequency modulation, the sine part."""
    return -np.imag(reflected_field(omg, r, fsr))


def cav_trans_full(nu: np_float, alpha: np_float, refr: np_float,
                   R: np_float, L: float, f0: float) -> np_float:
    """Cavity transmission function.
    
    Eq. 5 from A. Cygan, D. Lisak, P. Morzyński, M. Bober, M. Zawada,
    E. Pazderski, and R. Ciuryło, Opt. Express 21, 29744 (2013).

    This form does not assume a Lorentzian line shape.

    Parameters
    ---------
    nu : float or ndarray
        Optical frequency in Hz.
    alpha : float or ndarray
        Absorption coefficient in 1/cm.
    refr : float or ndarray
        Refractive index.
    R : float or ndarray
        Intensity reflection coefficient.
    L : float
        Cavity length in cm.
    f0 : float
        Offset frequency of the cavity (caused by Guoy phase shift).

    Returns
    -------
    float or ndarray
        Ratio of transmitted to incident light.
    """
    nu_fsr = C.c/(2*L*refr)
    top = (1-R)**2*np.exp(-alpha*L)
    bot = (1-R*np.exp(-alpha*L))**2 +\
        4*R*np.exp(-alpha*L)*np.sin(np.pi*(nu-f0)/nu_fsr)**2

    return top/bot
    

def cav_trans_exact(nu, alpha, nu0, gam, n, R, L):
    """Cavity transmission function for CMWS/CMDS.

    Eq. 5 from A. Cygan, D. Lisak, P. Morzyński, M. Bober, M. Zawada,
    E. Pazderski, and R. Ciuryło, Opt. Express 21, 29744 (2013).

    This form of `nu_fsr` assumes a Lorentzian absorption line shape.

    Parameters
    ----------
    nu : ndarray
        x-axis frequencies.
    alpha : ndarray
        Absorption coefficients corresponding to frequencies in `nu`.
    nu0 : float
        Center frequency of absorption line.
    gam : float
        Half-width at half-maximum of Lorentzian absorption line.
    n : float
        Refractive index without absorption.
    R : float
        Intensity reflection coefficient.
    L : float
        Cavity length (meters)

    Returns
    -------
    I : ndarray
        Ratio of transmitted to incident intensity.
    """
    nu_fsr = C.c/(2*L*n*(1-(alpha*C.c*(nu-nu0))/(4*np.pi*n*nu*gam)))
    print(nu_fsr)

    top = (1-R)**2*np.exp(-alpha*L)
    bot = (1-R*np.exp(-alpha*L))**2 +\
        4*R*np.exp(-alpha*L)*np.sin(np.pi*nu/nu_fsr)**2

    return top/bot


def cav_trans(R, a, L, disp=None):
    """This is exactly the same as in Maslowski MATLAB function, which is
    different from supplemental information of 'Quantum-noise-limited
    optical frequency comb spectroscopy'.
    """
    if disp is not None:
        delta = np.real(a)*L
        phi = np.imag(a)*L
        T = 1 - R
        x1 = T*np.exp(1j*disp/2 - delta/2 - 1j*phi/2)
        x2 = 1 - R*np.exp(1j*disp - delta - 1j*phi)
        y1 = T*np.exp(1j*disp/2)
        y2 = 1 - R*np.exp(1j*disp)
        x = np.abs(x1/x2)
        x = x**2
        y = np.abs(y1/y2)
        y = y**2

        return x/y
    else:
        return (1-R)**2*np.exp(-a*L)/((1-R*np.exp(-a*L))**2)


#######################################################################
# Short and simple functions                                          #
#######################################################################
def refraction2loss(wn, refraction):
    return refraction*4*np.pi/100/wn2lambda(wn)


def loss2refraction(wn, loss):
    """Convert effective loss per unit length to refractive index.

    Parameters
    ----------
    wn : float
        Wavenumbers in cm-1.
    loss : float
        Effective loss in 1/cm.

    Returns
    -------
    refraction : float
        Refractive index.
    """
    return loss*100*wn2lambda(wn)/4/np.pi


def loss2gam(loss):
    return 100*C.c*loss/4/np.pi


def gam2loss(gam):
    r"""Return effective loss per unit length in cm-1.

    Parameters
    ----------
    gam : float
        HWHM of cavity mode.

    Returns
    -------
    float
        Effective loss per unit length in cm-1.

    Notes
    -----
    The formula is as follows:

    .. math:: \alpha = \frac{4\pi n(\nu)}{c 100}\gamma

    So the value returned by this function is appropriate for vacuum.
    If background gases are present, then the result has to be
    multiplied by the refractive index.
    """
    return 4*np.pi*gam/C.c/100.0


def freq2phase(freq, frep):
    """Convert frequency offset to phase offset for PDH distortions."""
    return freq*np.pi/frep


def phase2freq(phase, frep):
    """Convert phase offset to frequency offset for PDH distortions."""
    return phase/np.pi*frep


def R_tau(tau, L):
    """Intensity reflection coeff. from ringdown time."""
    L1 = -np.sqrt((2*C.c*L*tau)**2 + L**4) + 2*(C.c*tau)**2 + L**2
    L2 = 2*(C.c*tau)**2

    return L1/L2


def Fin_R(R):
    """Finesse from intensity reflection coefficient."""
    return np.pi*np.sqrt(R)/(1-R)


def R_Fin(Fin):
    """Intensity reflection coefficient from finesse."""
    L1 = 2*Fin**2 - np.pi*np.sqrt((2*Fin)**2 + np.pi**2) + np.pi**2
    L2 = 2*Fin**2

    return L1/L2


def R_F(F):
    """Intensity reflection coefficient from coefficient of finesse."""
    return (np.sqrt(F**2 + 4) - 2)/F


def F_Fin(Fin):
    """Coefficient of finesse from finesse."""
    return 1/np.sin(np.pi/2/Fin)**2


# def cav_trans(R, a, L):
#     return (1-R)**2*np.exp(-a*L)/((1-R*np.exp(-a*L))**2)


def leff_trans(a, L):
    return np.exp(-a*L)


def gam_fsr2finesse(gam, fsr):
    """Calculate Finesse from FWHM and FSR.

    .. math:: \\mathcal{F} = \\frac{\\mathrm{FSR}}{\\mathrm{FWHM}} = \\pi \\frac{\\sqrt R}{1-R}

    Parameters
    ----------
    gam : float
        FWHM of the cavity mode in Hz,
    fsr : float
        Free spectral range of the cavity in Hz.

    Returns
    -------
    float
        Finesse of the cavity.
    """
    return fsr/gam

def circulating(r1, rm):
    """Fraction of circulating power.

    Parameters
    ----------
    r1 : float
        Field reflection from the input mirror.
    rm : float
        Field reflection from the rest.
    """
    t1 = np.sqrt(1-r1**2)

    return t1**2/(1-r1*rm)


def reflected(r1, rm):
    """Fraction of reflected power.

    Parameters
    ----------
    r1 : float
        Field reflection from the input mirror.
    rm : float
        Field reflection from the rest.
    """
    t1 = np.sqrt(1-r1**2)

    return (r1-rm)**2/(1-r1*rm)**2
