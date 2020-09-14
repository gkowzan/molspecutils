"""Auxilliary functions extending the HAPI library.
"""
import numpy as np
import spectroscopy.foreign.hapi3 as h3
from spectroscopy.data.pytips import (Tdat, TIPS_ISO_HASH, TIPS_GSI_HASH,
                                TIPS_NPT)


#: Boltzmann constant in erg/K
cBolts = 1.380648813E-16


def alpha2sw(alpha, conc):
    """Line strength from absoption coefficient and concentration.

    Parameters
    ----------
    alpha : float
        Absorption coefficient in cm-1.
    conc : float
        Volume concentration in cm-3.

    Returns
    -------
    float
        Line strength in cm-1/(molecule*cm-2)
    """
    return alpha/conc


def AtoB(aa, A, B, npt):
    """Lagrange 3- and 4-point interpolation.

    Parameters
    ----------
    aa : float
        Point at which `B` will be evaluated.
    A : ndarray
        x data.
    B : ndarray
        y data.
    npt : int
        Size of `A` and `B`

    Returns
    -------
    bb : float
        Interpolated value of `B` a `aa`.

    Notes
    -----
    Looks like a very literal translation from some Fortran code.
    """
    for I in range(2,npt+1):
        if A[I-1] >= aa:
            if I < 3 or I == npt:
                J = I
                if I < 3: J = 3
                if I == npt: J = npt
                J = J-1   # zero index correction
                A0D1=A[J-2]-A[J-1]
                if A0D1 == 0.0: A0D1=0.0001
                A0D2=A[J-2]-A[J]
                if A0D2 == 0.0: A0D2=0.0000
                A1D1=A[J-1]-A[J-2]
                if A1D1 == 0.0: A1D1=0.0001
                A1D2=A[J-1]-A[J]
                if A1D2 == 0.0: A1D2=0.0001
                A2D1=A[J]-A[J-2]
                if A2D1 == 0.0: A2D1=0.0001
                A2D2=A[J]-A[J-1]
                if A2D2 == 0.0: A2D2=0.0001

                A0=(aa-A[J-1])*(aa-A[J])/(A0D1*A0D2)
                A1=(aa-A[J-2])*(aa-A[J])/(A1D1*A1D2)
                A2=(aa-A[J-2])*(aa-A[J-1])/(A2D1*A2D2)

                bb = A0*B[J-2] + A1*B[J-1] + A2*B[J]

            else:
                J = I
                J = J-1   # zero index correction
                A0D1=A[J-2]-A[J-1]
                if A0D1 == 0.0: A0D1=0.0001
                A0D2=A[J-2]-A[J]
                if A0D2 == 0.0: A0D2=0.0001
                A0D3 = (A[J-2]-A[J+1])
                if A0D3 == 0.0: A0D3=0.0001
                A1D1=A[J-1]-A[J-2]
                if A1D1 == 0.0: A1D1=0.0001
                A1D2=A[J-1]-A[J]
                if A1D2 == 0.0: A1D2=0.0001
                A1D3 = A[J-1]-A[J+1]
                if A1D3 == 0.0: A1D3=0.0001

                A2D1=A[J]-A[J-2]
                if A2D1 == 0.0: A2D1=0.0001
                A2D2=A[J]-A[J-1]
                if A2D2 == 0.0: A2D2=0.0001
                A2D3 = A[J]-A[J+1]
                if A2D3 == 0.0: A2D3=0.0001

                A3D1 = A[J+1]-A[J-2]
                if A3D1 == 0.0: A3D1=0.0001
                A3D2 = A[J+1]-A[J-1]
                if A3D2 == 0.0: A3D2=0.0001
                A3D3 = A[J+1]-A[J]
                if A3D3 == 0.0: A3D3=0.0001

                A0=(aa-A[J-1])*(aa-A[J])*(aa-A[J+1])
                A0=A0/(A0D1*A0D2*A0D3)
                A1=(aa-A[J-2])*(aa-A[J])*(aa-A[J+1])
                A1=A1/(A1D1*A1D2*A1D3)
                A2=(aa-A[J-2])*(aa-A[J-1])*(aa-A[J+1])
                A2=A2/(A2D1*A2D2*A2D3)
                A3=(aa-A[J-2])*(aa-A[J-1])*(aa-A[J])
                A3=A3/(A3D1*A3D2*A3D3)

                bb = A0*B[J-2] + A1*B[J-1] + A2*B[J] + A3*B[J+1]

            break
    return bb


def BD_TIPS_2011_PYTHON(M, I, T):
    # out of temperature range
    if T<70. or T>3000.:
        raise Exception('TIPS: T must be between 70K and 3000K.')

    try:
        # get statistical weight for specified isotopologue
        gi = TIPS_GSI_HASH[(M,I)]
        # interpolate partition sum for specified isotopologue
        Qt = AtoB(T, Tdat, TIPS_ISO_HASH[(M, I)], TIPS_NPT)
    except KeyError:
        raise Exception('TIPS: no data for M,I = %d,%d.' % (M,I))

    return gi, Qt


def partitionSum(M, I, T, step=None):
    """Calculate range of partition sums at different temperatures.

    Output depends on a structure of input parameter T so that:

    - If T is a scalar/list and step IS NOT provided, then calculate
      partition sums over each value of T.
    - If T is a list and step parameter IS provided, then calculate
      partition sums between T[0] and T[1] with a given step.

    Parameters
    ----------
    M : int
        HITRAN molecule number
    I : int
        HITRAN isotopologue number
    T : float
        Temperature in Kelvin
    step: float, optional
        Step to calculate temperatures

    Returns
    -------
    TT : array_like, optional
        List of temperatures (present only if T is a list)
    PartSum : array_like
        Partition sums calculated on a list of temperatures

    References
    ----------
    [1] A. L. Laraia, R. R. Gamache, J. Lamouroux, I. E. Gordon,
    L. S. Rothman.  Total internal partition sums to support planetary
    remote sensing.  Icarus, Volume 215, Issue 1, September 2011, Pages
    391–400 http://dx.doi.org/10.1016/j.icarus.2011.06.004
    """
    # partitionSum
    if not step:
        if type(T) not in set([list, tuple]):
            return BD_TIPS_2011_PYTHON(M, I, T)[1]
        else:
            return [BD_TIPS_2011_PYTHON(M, I, temp)[1] for temp in T]
    else:
        TT = np.arange(T[0], T[1], step)
        return TT, np.array([BD_TIPS_2011_PYTHON(M, I, temp)[1] for temp in TT])


def PYTIPS(M, I, T):
    """Return total internal partition sum.

    Parameters
    ----------
    M : int
        Molecule id.
    I : int
        Isotopologue id.
    T : float
        Temperature (Kelvin)

    Returns
    -------
    TIPS : float
        Total internal partition sum.
    """
    return BD_TIPS_2011_PYTHON(M, I, T)[1]


def EnvironmentDependency_Intensity(
        LineIntensityRef, T, Tref, SigmaT, SigmaTref,
        LowerStateEnergy, LineCenter):
    """Return line intensity at given temperature `T`.

    Parameters
    ----------
    LineIntensityRef : float
        HITRAN line intensity at `Tref`.
    T : float
        Temperature for returned line intensity (Kelvin).
    Tref : float
        Reference temperature (Kelvin).
    SigmaT : float
        TIPS at temperature `T`.
    SigmaTref : float
        TIPS at temperature `Tref`.
    LowerStateEnergy: float
        Lower state energy (cm-1).
    LineCenter : float
        Transition energy (cm-1).

    Returns
    -------
    LineIntensity : float
        Line intensity at temperature `T`.
    """
    const = np.float64(1.4388028496642257)
    ch = np.exp(-const*LowerStateEnergy/T)*(1-np.exp(-const*LineCenter/T))
    zn = np.exp(-const*LowerStateEnergy/Tref)*(1-np.exp(-const*LineCenter/Tref))
    LineIntensity = LineIntensityRef*SigmaTref/SigmaT*ch/zn

    return LineIntensity


def EnvironmentDependency_Gamma0(Gamma0_ref, T, Tref, p, pref,
                                 TempRatioPower):
    return Gamma0_ref*p/pref*(Tref/T)**TempRatioPower


def EnvironmentDependency_Delta0(Delta0_ref, p, pref):
    return Delta0_ref*p/pref


def volumeConcentration(p, T):
    """Return concentration (per unit volume).

    Parameters
    ----------
    p : float
        Pressure (atm)
    T : float
        Temperature (Kelvin)
    """
    return (p/9.869233e-7)/(cBolts*T) # CGS


def cross_sections(table, T, p):
    """Calculate cross sections for lines in `table` at given T and p.

    Parameters
    ----------
    table : ndarray
        HAPI table with line parameters.
    T : float
        Temperature in Kelvin.
    p : float
        Pressure in atms.

    Returns
    -------
    lines : ndarray
        2D array with columns: (`nu`, `xs`), where `nu` is the
        transition frequency in cm-1 and `xs` is the cross section in
        cm**2.  The number of rows is the same as in `table`, i.e. each
        row corresponds to a different transition.
    """
    Tref = 296.0
    nline = table['header']['number_of_rows']

    lines = np.empty((nline, 2))
    for row_num in range(nline):
        LineCenterDB = table['data']['nu'][row_num]
        LineIntensityDB = table['data']['sw'][row_num]
        LowerStateEnergyDB = table['data']['elower'][row_num]
        MoleculeNumberDB = table['data']['molec_id'][row_num]
        IsoNumberDB = table['data']['local_iso_id'][row_num]

        SigmaT = PYTIPS(MoleculeNumberDB, IsoNumberDB, T)
        SigmaTref = PYTIPS(MoleculeNumberDB, IsoNumberDB, Tref)
        LineIntensity = EnvironmentDependency_Intensity(
            LineIntensityDB, T, Tref, SigmaT, SigmaTref,
            LowerStateEnergyDB, LineCenterDB)

        lines[row_num] = (LineCenterDB, LineIntensity)

    lines[:, 1] *= volumeConcentration(p, T)

    return lines


def init_params_p(table_name, p, T):
    """Prepare initial parameters for fitting spectra.

    Gammas and deltas are pressure-independent, cros sections are given
    for specific pressure.

    Parameters
    ----------
    table_name : str
        HAPI table name with molecular data.
    T : float
        Temperature in Kelvin.
    p : float
        Pressure in atm.

    Returns
    -------
    lines : list
        Names of transitions.
    params : ndarray
        ndarray with columns :math:`\\sigma`, :math:`\\nu_0`,
        :math:`\\gamma`.
    """
    Tref = 296.0

    nus = np.array(h3.getColumn(table_name, 'nu'))
    deltas = np.array(h3.getColumn(table_name, 'delta_air'))
    gammas = np.array(h3.getColumn(table_name, 'gamma_air'))
    n_air = np.array(h3.getColumn(table_name, 'n_air'))

    params = np.empty((nus.size, 4))
    params[:, 2] = EnvironmentDependency_Gamma0(
        gammas, T, Tref, p, 1.0, n_air
    )                           # halfwidth
    params[:, 0] = cross_sections(
        h3.LOCAL_TABLE_CACHE[table_name], T, p
    )[:, 1]                     # cross section
    params[:, 1] = nus
    params[:, 3] = EnvironmentDependency_Delta0(
        deltas, 1.0, 1.0
    )

    lines = h3.LOCAL_TABLE_CACHE[table_name]['data']['local_lower_quanta']
    lines = [l.strip().replace(' ', '') for l in lines]

    return lines, params


def init_params(table_name, p, T):
    """Prepare initial parameters for fitting spectra.

    Parameters
    ----------
    table_name : str
        HAPI table name with molecular data.
    p : float
        Pressure in atm.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    lines : list
        Names of transitions.
    params : ndarray
        ndarray with columns :math:`\\sigma`, :math:`\\nu_0`,
        :math:`\\gamma`.
    """
    Tref = 296.0

    nus = np.array(h3.getColumn(table_name, 'nu'))
    deltas = np.array(h3.getColumn(table_name, 'delta_air'))
    gammas = np.array(h3.getColumn(table_name, 'gamma_air'))
    n_air = np.array(h3.getColumn(table_name, 'n_air'))

    params = np.empty((nus.size, 3))
    params[:, 2] = EnvironmentDependency_Gamma0(
        gammas, T, Tref, p, 1.0, n_air
    )                           # halfwidth
    params[:, 0] = cross_sections(
        h3.LOCAL_TABLE_CACHE[table_name], T, p
    )[:, 1]                     # cross section
    params[:, 1] = nus + EnvironmentDependency_Delta0(
        deltas, p, 1.0
    )

    lines = h3.LOCAL_TABLE_CACHE[table_name]['data']['local_lower_quanta']
    lines = [l.strip().replace(' ', '') for l in lines]

    return lines, params
