"""Auxilliary functions extending the HAPI library. """
# * Imports and constants
from pathlib import Path
import sys
from collections import namedtuple
import importlib.util
import re
import appdirs
import numpy as np
import scipy.constants as C
from molspecutils.data.pytips import (Tdat, TIPS_ISO_HASH, TIPS_GSI_HASH,
                                      TIPS_NPT)
import molspecutils.utils as u

def lazy(fullname):
  try:
    return sys.modules[fullname]
  except KeyError:
    spec = importlib.util.find_spec(fullname)
    module = importlib.util.module_from_spec(spec)
    loader = importlib.util.LazyLoader(spec.loader)
    # Make module with proper locking and get it inserted into sys.modules.
    loader.exec_module(module)
    return module

#: Lazy import of hapi3
h3 = lazy('molspecutils.foreign.hapi3')

dirs = appdirs.AppDirs('happier', 'gkowzan')
hitran_cache = str(Path(dirs.user_cache_dir) / 'db')
#: Boltzmann constant in erg/K
cBolts = 1.380648813E-16
Tref = 296.0

# * Retrieve molecular data from HITRAN
HITRANRow = namedtuple("HITRANRow", "llq, luq, glq, guq, nu, elower, sw, a, gair, gself, dair, nair, gp, gpp")


def molname2molid(molname: str) -> int:
    for v in h3.ISO_ID.values():
        if v[5] == molname:
            return v[0]

digits = [str(x) for x in range(1, 10)]

def CH3Cl_nu_to_dict(snu):
    if snu[0] in digits:
        return {'nu'+snu[2]: int(snu[0])}
    else:
        return {'nu'+snu[1]: 1}


def CH3Cl_gq_to_dict(gq):
    """Convert global quanta to dict."""
    gdict = dict(nu1=0, nu2=0, nu3=0, nu4=0, nu5=0, nu6=0)
    gq = gq.strip()
    if not gq == 'GROUND':
        if '+' in gq:
            nuf, nul = gq.split('+')
            gdict.update(CH3Cl_nu_to_dict(nuf))
            gdict.update(CH3Cl_nu_to_dict(nul))
        else:
            gdict.update(CH3Cl_nu_to_dict(gq))

    return gdict


def CH3Cl_lq_to_dict(lq):
    """Convert local quanta to dict."""
    fstring = lq[11:].strip()
    if fstring:
        f = float(fstring)
    else:
        f = 0.0
    return dict(J=int(lq[:3].strip()), K=int(lq[3:6].strip()),
                l=int(lq[6:8].strip()), rovib=lq[8:10].strip(),
                f=f)
            

def CO_llq_to_pair(llq):
    """Convert `local_lower_quanta` string to (Jpp, Jp)."""
    match = re.fullmatch('\s*([PR])\s*([0-9]+)\s*', llq)
    branch, j = match.group(1), int(match.group(2))
    if branch=='R':
        jp = j+1
    elif branch=='P':
        jp = j-1
    else:
        raise ValueError('Branch is neither "R" nor "P".')

    return (j, jp)


def energy_levels(nus: list, M: int, I: int, db: str=hitran_cache):
    """Retrieve energy levels from HITRAN for.

    Parameters
    ----------
    nu: list of int
        Vibrational numbers.
    M: int
        HITRAN molecule id.
    I: int
        HITRAN isotopologue id for molecule `M`.
    db: str, optional
        Path to the database folder.

    Returns
    -------
    levels: dict
        keys are tuples (nu, j), values are energy levels in cm-1
    """
    table_name = "{:d}_{:d}".format(M, I)
    h3.db_begin(db)
    if table_name not in h3.getTableList():
        h3.fetch(table_name, M, I, 0, 99999)

    levels = {}
    for nu in nus:
        conds = ('AND', ('=', 'molec_id', M), ('=', 'local_iso_id', I),
                 ('MATCH', ('STR', r'\s+{:d}\s*'.format(nu)), 'global_lower_quanta'),
                 ('MATCH', ('STR', r'\s+{:d}\s*'.format(nu)), 'global_upper_quanta'))
        dest_table = '{:s}_nu{:d}'.format(table_name, nu)
        h3.select(table_name, Conditions=conds, DestinationTableName=dest_table, Output=False)
        llqs, elowers = h3.getColumns(dest_table, ParameterNames=('local_lower_quanta', 'elower'))
        for llq, elower in zip(llqs, elowers):
            levels[(nu, CO_llq_to_pair(llq)[0])] = elower

    return levels


def equilibrium_pops(levels: dict, T: float, M: int, I: int):
    """Calculate equilibrium populations at temperature T.

    This is population summed over degenerate magnetic levels.
    """
    kt = u.joule2wn(C.k*T)
    # return {k: (2*k[1]+1)*np.exp(-v/kt)/PYTIPS(M, I, T) for k, v in levels.items()}
    # adding square root factor, don't remember why, have to check Vaccaro
    return {k: np.sqrt(2*k[1]+1)*np.exp(-v/kt)/PYTIPS(M, I, T) for k, v in levels.items()}


def line_params(bands: list, M: int, I: int, db: str=hitran_cache):
    r"""Collect parameters of rovib transitions.

    `mu` parameter is :math:`\sqrt{S_{J',J''}}=\sqrt{3\frac{\epsilon_0\hbar\pi c^3}{\omega^3}A_{J'\to J''}(2J'+1)}`.

    Parameters
    ----------
    bands: list of tuple
        List of (nupp, nup) tuples.
    M: int
        HITRAN molecule id.
    I: int
        HITRAN isotopologue id.
    db: str
        Path to HITRAN db.

    Returns
    -------
    params: dict of dict
        Nested dict: ((nupp, jpp), (nup, jp)) -> parameter -> value, where
        parameter is one of: 'mu' , 'gam', 'del', 'nu'.
    """
    table_name = "{:d}_{:d}".format(M, I)
    h3.db_begin(db)
    if table_name not in h3.getTableList():
        h3.fetch(table_name, M, I, 0, 99999)

    params = {}
    for nupp, nup in bands:
        conds = ('AND', ('=', 'molec_id', M), ('=', 'local_iso_id', I),
                 ('MATCH', ('STR', r'\s+{:d}\s*'.format(nupp)), 'global_lower_quanta'),
                 ('MATCH', ('STR', r'\s+{:d}\s*'.format(nup)), 'global_upper_quanta'))
        dest_name = '{:s}_band_{:d}_{:d}'.format(table_name, nupp, nup)
        h3.select(table_name, DestinationTableName=dest_name, Conditions=conds, Output=False)
        param_list = ('local_lower_quanta', 'a', 'gamma_air', 'delta_air', 'nu', 'n_air')
        llqs, As, gams, delts, nus, n_airs = h3.getColumns(dest_name, ParameterNames=param_list)
        for llq, a, gam, delt, nu, n_air in zip(llqs, As, gams, delts, nus, n_airs):
            jpp, jp = CO_llq_to_pair(llq)
            # mu = np.sqrt(a*(2*jp+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*u.wn2nu(nu))**3)
            # See Eqs. (5.4) and (3.33) in rotsim2d_roadmap.pdf. `mu` here
            # includes the Hönl-London factor, so it is actually sqrt(S_j2,j1/(2j_1+1))
            # and not mu_j2,j1
            # mu = np.sqrt(a*(2*jp+1)/(2*jpp+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*u.wn2nu(nu))**3)
            mu = np.sqrt(a*(2*jp+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*u.wn2nu(nu))**3)
            params[((nupp, jpp), (nup, jp))] = {
                'mu': mu, 'gam': gam, 'del': delt, 'n_air': n_air ,'nu': nu}
    return params

def generate_line_params(bands: list, elevels: dict, line_params: dict):
    """Generate line parameters based on HITRAN data.

    Parameters
    ----------
    bands: list of tuple
        (nupp, nup, dj) tuples.
    elevels: dict
        Result of :func:`energy_levels`.
    line_params: dict
        Result of :func:`line_params`.

    Returns
    -------
    params: dict of dict
        Nested dict: ((nupp, jpp), (nup, jp)) -> parameter -> value, where
        parameter is one of: 'mu' , 'gam', 'del', 'nu'.
    """
    params = {}

    # generate line positions
    for nupp, nup, dj in bands:
        if nupp == nup:
            djs = (dj, )
        else:
            djs = (dj, -dj)
        for nu, j in elevels:
            if nu == nupp:
                for dj in djs:
                    if (nup, j+dj) in elevels:
                        params.setdefault(((nupp, j), (nup, j+dj)), {})['nu'] = elevels[(nup, j+dj)]-elevels[(nupp, j)]

    # generate line parameters, remove lines for which no parameters can be found
    for line in list(params.keys()):
        found = False
        for old_line in line_params:
            if line[0] == old_line[0] and line[0][0] == old_line[0][0] and line[1][0] == old_line[1][0]:
                params[line]['gam'] = line_params[old_line]['gam']
                params[line]['n_air'] = line_params[old_line]['n_air']
                params[line]['del'] = 0.0
                # this is almost completely wrong
                params[line]['mu'] = line_params[old_line]['mu']
                found = True
                break
        if not found:
            del params[line]

    return params

def calc_relaxation():
    """Calculate rotational relaxation rate for diagonal elements."""
    return NotImplementedError

def apply_env(line_params: dict, env: dict, to_nu: bool=False):
    """Apply environment (p, T) dependency to line parameters."""
    ret = {}
    for k, v in line_params.items():
        gam = EnvironmentDependency_Gamma0(v['gam'], env['T'], Tref, env['p'], 1.0, v['n_air'])
        nu = v['nu'] + EnvironmentDependency_Delta0(v['del'], env['p'], 1.0)
        if to_nu:
            gam = u.wn2nu(gam)
            nu = u.wn2nu(nu)
        ret[k] = {'gam': gam, 'nu': nu}
        if 'mu' in v:
            ret[k]['mu'] = v['mu']

    return ret

# * For line-shape calculations
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
