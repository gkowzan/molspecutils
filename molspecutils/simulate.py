"""Simulate frequency comb pulses."""
# * Notes
# Functions which should logically accept frequencies or wavelengths accept frequencies in Hz. Functions returning amplitudes or powers return appropriately normalized spectral densities.  Time is expressed in seconds.
# * Imports
import numpy as np
import scipy.constants as C
from scipy.stats import linregress
from scipy.signal import detrend
from scipy.misc import derivative
from knickknacks.experiment import find_index
import pyfftw.interfaces.scipy_fftpack as fftp


# * Functions
# ** Analyze pulses
def carrier_freq(interf: np.ndarray, debug: bool=False):
    """Get carrier frequency from phase ramp."""
    ih = hilbert(interf)
    iph = np.unwrap(np.angle(ih))
    ix = np.arange(iph.size)
    regres = linregress(ix, iph)

    if debug:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(ix, regres.slope*ix+regres.intercept)
        ax.plot(iph, 'o', ms=1)

    return regres.slope


def hilbert(x, N=None, axis=-1):
    """
    Compute the analytic signal, using the Hilbert transform.
    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.
    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`
    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = fftp.fft(x, N, axis=axis, planner_effort='FFTW_ESTIMATE')
    h = np.zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if x.ndim > 1:
        ind = [np.newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    x = fftp.ifft(Xf*h, axis=axis, planner_effort='FFTW_ESTIMATE')

    return x

# ** Generate pulses
def p_gauss(t, tg, E_0=1.0, a=None, carrier=None):
    """Gaussian pulse with linear chirp normalized to :math:`E_0`."""
    norm = 1/np.sqrt(np.pi)/tg

    pulse = norm*E_0*np.exp(-(t/tg)**2)

    if carrier is not None:
        if a is not None:
            return pulse*np.cos(2*np.pi*carrier*t-a*(t/tg)**2)
        return pulse*np.cos(2*np.pi*carrier*t)
    else:
        return pulse

def carrier_align(carrier, tr):
    """Align carrier to tr to get zero-offset comb."""
    Tc = 1/carrier

    return 1/(tr/np.round(tr/Tc))


def t_comb_adjust(t, tr):
    """Adjust time axis for pulse train.

    Modifies `t` in place.
    """
    to_shift = np.where(t>tr/2)[0]
    if to_shift.size == 0:
        to_shift = np.where(t<-tr/2)[0]
        if to_shift.size == 0:
            return t
        else:
            t[to_shift] += tr
            return t_comb_adjust(t, tr)
    else:
        t[to_shift] -= tr
        return t_comb_adjust(t, tr)
    

def p_gauss_comb(t, tg, tr, carrier, E_0=1.0, a=None):
    carrier = carrier_align(carrier, tr)
    # t = t_comb_adjust(t, tr)

    return p_gauss(t_comb_adjust(t.copy(), tr), tg, E_0, a, carrier)


def p_gauss_autocorrelation(t, tg, E_0=1.0, carrier=None):
    """Gaussian pulse autocorrelation."""
    norm = 1/np.sqrt(2*np.pi)/tg
    pulse = norm*E_0**2*np.exp(-(t**2/2/tg**2))

    if carrier is not None:
        return pulse*np.exp(1.0j*2*np.pi*carrier*t)
    else:
        return pulse


def comb_autocorrelation(t, tg, tr, carrier):
    """Frequency comb autocorrelation."""
    t = np.modf(t/tr)[0]*tr - tr/2

    return p_gauss_autocorrelation(t, tg, carrier=carrier)

# ** Modify pulses
def add_phase_shift(fs, ft, ps):
    """Add spectral phase shift to FT of a pulse.

    Parameters
    ----------
    fs : ndarray
        Frequency axis in Hz.
    ft : ndarray
        Fourier transform of a pulse.
    ps : callable
        Returns phase shift for given frequencies.

    Returns
    -------
    rft : ndarray
        Fourier transform with applied phase shift.
    """
    return ft*np.exp(1.0j*ps(fs))


# ** Materials
def n_bk7(f):
    """Index of refraction for N-BK7 glass.

    Parameters
    ----------
    f : float
        Frequency in Hz.
    """
    x = C.c/f*1e6

    n = (1+1.03961212/(1-0.00600069867/x**2)+0.231792344/(1-0.0200179144/x**2)+1.01046945/(1-103.560653/x**2))**.5

    return n


def beta_coeffs(n, omega0, order=2):
    """Return beta coeffs based on refractive index function.

    Parameters
    ----------
    n : callable
        Refractive index function, frequency -> refractive index.
    omega0 : float
        Central angular frequency.
    order : int
        Highest order of beta parameter.

    Returns
    -------
    betas : ndarray
        Array of beta coefficients.
    """
    def beta(omega):
        return omega*n(omega/2/np.pi)/C.c

    return np.array([derivative(beta, omega0, dx=2.0*np.pi*100e6, n=cur_order, order=101)
                     for cur_order in range(order+1)])


def spectral_ps(fmin, fmax, d, n, nonlinear=True):
    """Return callable returning spectral phase shift.

    Parameters
    ----------
    fmin : float
        Minimum frequency.
    fmax : float
        Maximum frequency.
    d : float
        Thickness of the material.
    n : callable
        Callable return refractive index of the material for a frequency
        in Hz.
    nonlinear : bool
        Remove linear part.

    Returns
    -------
    ps : callable
        Function accepting frequencies and returning phase shift between
        `fmin` and `fmax` for a material of thickness `d`.  Returns 0.0
        for frequencies beyond the defined frequencies range.
    """
    def ps(f):
        fs = np.atleast_1d(f)
        rps = np.zeros(fs.size)
        imin, imax = sorted((find_index(fmin, fs), find_index(fmax, fs)))
        rps[imin:imax] = detrend(-2*np.pi*n(fs[imin:imax])*fs[imin:imax]*d/C.c)
        rps[imin:imax] = np.remainder(rps[imin:imax], -2*np.pi)
        imin, imax = sorted((find_index(-fmin, fs), find_index(-fmax, fs)))
        rps[imin:imax] = detrend(-2*np.pi*n(fs[imin:imax])*fs[imin:imax]*d/C.c)
        rps[imin:imax] = np.remainder(rps[imin:imax], -2*np.pi)

        return rps

    return ps


# ** Noise
def fftnoise(f):
    f = np.array(f, dtype='complex')
    Np = (len(f) - 1) // 2
    phases = np.random.rand(Np) * 2 * np.pi
    phases = np.cos(phases) + 1j * np.sin(phases)
    f[1:Np+1] *= phases
    f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
    return np.fft.ifft(f).real


def low_pass(data, max_freq, samplerate):
    data_fft = fftp.fft(data)
    freqs = np.abs(fftp.fftfreq(data.size, 1/samplerate))
    idx = np.where(freqs<=max_freq)[0]
    data_fft[idx] = 0.0

    return fftp.ifft(data_fft)


def band_limited_noise(min_freq, max_freq, samples=1024, samplerate=1):
    freqs = np.abs(np.fft.fftfreq(samples, 1/samplerate))
    f = np.zeros(samples)
    idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
    f[idx] = 1
    return fftnoise(f)


def oneoverf(f, knee, alpha):
    desc = np.ones_like(f)
    desc[f<knee] = np.abs((f[f<knee]/knee)**(-alpha))
    desc[0] = 1

    return desc


def oneoverf_noise(sigma, knee, alpha, samples, samplerate):
    # https://gist.github.com/zonca/979729
    wn = np.random.normal(0.0, sigma, samples)
    s = fftp.rfft(wn)
    f = fftp.rfftfreq(samples, d=1/samplerate)[:len(s)]
    f[-1] = np.abs(f[-1])
    fft_sim = s*oneoverf(f, knee, alpha)

    return fftp.irfft(fft_sim)


# ** Utilities
def lambda_span(lam, nu_span):
    """Span of wavelengths around `lam` corresponding to `nu_span`."""
    return lam**2/C.c*nu_span
