"""Classes representing rovibrational molecular states and vibrational modes."""
import abc
import logging
from typing import Mapping, Tuple, Union

import numpy as np
import scipy.constants as C
from attrs import define, field
from sqlalchemy import select
from sqlalchemy.engine import Engine

import molspecutils.alchemy.C2H2_nu1nu3 as C2H2_nu1nu3
import molspecutils.alchemy.CH3Cl_nu3 as CH3Cl_nu3
import molspecutils.alchemy.CO as CO
import molspecutils.happier as hap
import molspecutils.utils as u
from molspecutils.alchemy.convert import get

log = logging.getLogger('__name__')

class MissingStateError(Exception):
    pass

class MissingLineError(Exception):
    pass


@define(frozen=True, cache_hash=True, order=True)
class DiatomState:
    nu: int
    j: int
    name: str = field(init=False, eq=False, repr=False)

    def __attrs_post_init__(self):
        object.__setattr__(self, "name",
                           "{:d},{:d}".format(self.nu, self.j))

    @classmethod
    def from_symtop(cls, symtop: "SymTopState"):
        return cls(symtop.nu, symtop.j)


@define(frozen=True, cache_hash=True, order=True)
class SymTopState:
    nu: int
    j: int
    k: int
    name: str = field(init=False, eq=False, repr=False)

    def __attrs_post_init__(self):
        object.__setattr__(self, "name",
                           "{:d},{:d},{:d}".format(self.nu, self.j, self.k))

    @classmethod
    def from_diatom(cls, diatom: DiatomState, k: int):
        return cls(diatom.nu, diatom.j, k)


RotState = Union[DiatomState, SymTopState]


class VibrationalMode(abc.ABC):
    """Interface of a class representing a molecular vibrational mode."""
    @abc.abstractmethod
    def gamma(self, pair: Tuple[RotState, RotState]) -> float:
        """Pressure broadening coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def delta(self, pair: Tuple[RotState, RotState]) -> float:
        """Pressure shift coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def mu(self, pair: Tuple[RotState, RotState]) -> float:
        """Reduced matrix element for dipole transition between `pair` states."""

    @abc.abstractmethod
    def tips(self, T: float) -> float:
        """Total internal partition function."""

    @abc.abstractmethod
    def energy(self, state: RotState) -> float:
        """Return energy of `state`."""
        pass

    @abc.abstractmethod
    def degeneracy(self, state: RotState) -> float:
        """Return quantum state degeneracy."""

    def nu(self, pair: Tuple[RotState, RotState]) -> float:
        try:
            return u.wn2nu(self.energy(pair[1])-self.energy(pair[0]))
        except KeyError as e:
            raise MissingStateError("Energy for state `{!s}` is not available".format(
                e.args[0]))

    def difference_pop(self, pair: Tuple[RotState, RotState], T: float) -> float:
        r"""Population difference factor between nondegenerate states.

        .. math::

            \frac{N_{\nu_1,j_1,d_1}-N_{\nu_2,j_2,d_2}}{N} = \frac{1}{Z} e^{-E_{\nu_1,j_1/kT}} \left( 1 - e^{-h\nu/kT} \right)

        where :math:`d_i` are non-rotational degenerate states labels.
        """
        e1, e2 = self.energy(pair[0]), self.energy(pair[1])
        ediff = e2-e1
        kT = C.k*T/C.h/100.0/C.c

        return np.exp(-e1/kT)*(1 - np.exp(-ediff/kT))/self.tips(T)


    def equilibrium_pop(self, state: RotState, T: float) -> float:
        e = self.energy(state)
        kt = u.joule2wn(C.k*T)

        return self.degeneracy(state)*np.exp(-e/kt)/self.tips(T)


@define
class LineParams:
    A: float
    gamma: float
    delta: float
    sw: float


class AlchemyModeMixin(VibrationalMode):
    """Molecular vibrational mode backed by SQLAlchemy sqlite database.

    Subclasses need to define and fill the following attributes:

    - ``lines``, dict from 2-tuples of :class:`RotState` to :class:`LineParams`
    - ``elevels``, dict from :class:`RotState` to energy
    - ``degeneracies``, dict from :class:`RotState` to quantum state degeneracy
    """
    def line_params(self, pair: Tuple[RotState, RotState]) -> Tuple[LineParams, int]:
        result = self.lines.get(pair)
        if result:
            return (result, 1)
        result = self.lines.get(pair[::-1])
        if result is None:
            log.debug("Missing line parameters for: %r", pair)
            result = LineParams(0, 0, 0, 0)
        return (result, -1)

    def gamma(self, pair: Tuple[RotState, RotState]) -> float:
        params, _ = self.line_params(pair)
        if params is None or params.gamma == 0:
            gam = self._fake_gamma(pair)
        else:
            gam = params.gamma

        return u.wn2nu(gam)

    def sw(self, pair: Tuple[RotState, RotState]) -> float:
        "HITRAN line strength in cm-1 cm**2."
        return self.line_params(pair)[0].sw

    def delta(self, pair: Tuple[RotState, RotState]) -> float:
        params, _ = self.line_params(pair)
        delt = params.delta

        return u.wn2nu(delt)

    def mu(self, pair: Tuple[RotState, RotState]) -> float:
        r"""Reduced matrix element for `pair[0]` to `pair[1]` transition.

        Obtained from HITRAN's Einsten A-coefficient:

        .. math::

            \sqrt{S(J', J'')} |\langle \nu' | \mu | \nu'' \rangle| =
            |\langle \nu';J'\|\mu^{(1)}\|\nu'';J''\rangle|^{2} =
            A_{\nu'J'\to\nu''J''}\frac{\epsilon_{0}hc^{3}d'}{16\pi^{3}\nu^{3}_{\nu'J',\nu''J''}}

        where :math:`d'` is total degeneracy of the upper state.
        """
        return rme(pair, self)

    def energy(self, state: RotState) -> float:
        """Return energy of ``state``."""
        return self.elevels[state]

    def degeneracy(self, state: RotState) -> float:
        """Return quantum state degeneracy."""
        return self.degeneracies[state]

    def _fake_gamma(self, pair: Tuple[RotState, RotState]) -> float:
        for kpp, kp in self.lines.keys():
            if pair[0]==kp or pair[0]==kpp or pair[1]==kp or pair[1]==kpp:
                return self.lines[(kpp, kp)].gamma
        raise ValueError("Can't generate fake pressure width for coherence '{:s}".format(pair))


class CH3ClAlchemyMode(AlchemyModeMixin):
    molecule = 'CH3Cl'
    def __init__(self, engine_or_path=None, iso=1):
        """Provide transition and energy level data of nu3 mode.

        If `engine_or_path` is a directory, then it will be searched for sqlite3
        database with appropriate structure (see
        :mod:`molspecutils.alchemy.CH3Cl`). If not present, it will search for
        HAPI db file `CH3Cl_nu3_<iso>.data` and attempt to extract the
        data and put it in sqlite3 db. If neither is present, it will fetch the
        data from HITRAN first. If `engine_or_path` is None, it will do the
        above for default user cache directory. If `engine_or_path` is a
        :class:`sqlalchemy.Engine` instance, it will be assumed to contain the
        required molecular data.

        Parameters
        ----------
        engine_or_path : str
            Path to directory with HAPI or sqlite3 database or
            :class:`sqlachemy.Engine` instance of opened sqlite3 database.
        iso : int
            Isotopologue number, 1 for 35Cl and 2 for 37Cl. Required for
            fetching data and calculating appropriate total partition function.
        """
        if isinstance(engine_or_path, Engine):
            engine = engine_or_path
        else:
            engine = get(engine_or_path, 'CH3Cl_nu3', iso)

        self._iso = iso
        self._generate_elevels(engine)
        self._generate_lines(engine)

    def _generate_elevels(self, engine: Engine):
        self.elevels = {}
        self.degeneracies = {}
        with engine.begin() as conn:
            st = CH3Cl_nu3.states
            result = conn.execute(
                select(st.c.energy, st.c.g, st.c.nu3, st.c.J, st.c.K))
            for row in result:
                self.elevels[SymTopState(*row[2:])] = row[0]
                self.degeneracies[SymTopState(*row[2:])] = row[1]

    def _generate_lines(self, engine: Engine):
        self.lines = {}
        with engine.begin() as conn:
            lp = CH3Cl_nu3.line_parameters
            statepp = CH3Cl_nu3.states.alias()
            statep = CH3Cl_nu3.states.alias()
            result = conn.execute(
                select(lp.c.a, lp.c.gair, lp.c.dair, lp.c.sw,
                       statepp.c.nu3, statepp.c.J, statepp.c.K,
                       statep.c.nu3, statep.c.J, statep.c.K).\
                join_from(lp, statepp, lp.c.statepp==statepp.c.id).\
                join_from(lp, statep, lp.c.statep==statep.c.id))

            for row in result:
                spp = SymTopState(*row[4:7])
                sp = SymTopState(*row[7:10])
                self.lines[(spp, sp)] = LineParams(*row[:4])

    def _custom_degeneracy(self, state: RotState) -> float:
        gnuc = 16
        k3 = 2 if state.k > 0 and state.k % 3 == 0 else 1

        return gnuc*k3*(2*state.j+1)

    def equilibrium_pop(self, state: RotState, T: float) -> float:
        return super().equilibrium_pop(state, T)*0.5

    def tips(self, T: float) -> float:
        """Total internal partition function."""
        return hap.PYTIPS(24, self._iso, T)


class COAlchemyMode(AlchemyModeMixin):
    molecule = "CO"
    def __init__(self, engine_or_path=None, iso=1):
        """Provide transition and energy level data for CO.

        Parameters
        ----------
        engine_or_path
            Path to directory with HAPI or sqlite3 database or
            :class:`sqlachemy.Engine` instance of opened sqlite3 database.
        iso
            Isotopologue number. Required for fetching data and calculating
            appropriate total partition function.
        """
        if isinstance(engine_or_path, Engine):
            engine = engine_or_path
        else:
            engine = get(engine_or_path, 'CO', iso)

        self._iso = iso
        self._generate_lines(engine)
        self._generate_elevels(engine)

    def _generate_lines(self, engine: Engine):
        self.lines = {}
        with engine.begin() as conn:
            lp = CO.line_parameters
            statepp = CO.states.alias()
            statep = CO.states.alias()
            result = conn.execute(
                select(lp.c.a, lp.c.gair, lp.c.dair, lp.c.sw,
                       statepp.c.nu, statepp.c.J,
                       statep.c.nu, statep.c.J).\
                join_from(lp, statepp, lp.c.statepp==statepp.c.id).\
                join_from(lp, statep, lp.c.statep==statep.c.id))

            for row in result:
                spp = DiatomState(*row[4:6])
                sp = DiatomState(*row[6:8])
                self.lines[(spp, sp)] = LineParams(*row[:4])

    def _generate_elevels(self, engine: Engine):
        self.elevels = {}
        self.degeneracies = {}
        with engine.begin() as conn:
            st = CO.states
            result = conn.execute(
                select(st.c.energy, st.c.g, st.c.nu, st.c.J))
            for row in result:
                self.elevels[DiatomState(*row[2:])] = row[0]
                self.degeneracies[DiatomState(*row[2:])] = row[1]

    def _custom_degeneracy(self, state: RotState) -> float:
        return 2*state.j+1

    def tips(self, T: float) -> float:
        """Total internal partition function."""
        return hap.PYTIPS(5, self._iso, T)


class C2H2AlchemyMode(AlchemyModeMixin):
    molecule = "C2H2"
    def __init__(self, engine_or_path=None, iso=1):
        """Provide transition and energy level data for C2H2 nu1+nu3 mode.

        Parameters
        ----------
        engine_or_path
            Path to directory with HAPI or sqlite3 database or
            :class:`sqlachemy.Engine` instance of opened sqlite3 database.
        iso
            Isotopologue number. Required for fetching data and calculating
            appropriate total partition function.
        """
        if isinstance(engine_or_path, Engine):
            engine = engine_or_path
        else:
            engine = get(engine_or_path, 'C2H2_nu1nu3', iso)

        self._iso = iso
        self._generate_lines(engine)
        self._generate_elevels(engine)

    def _generate_lines(self, engine: Engine):
        self.lines = {}
        with engine.begin() as conn:
            lp = C2H2_nu1nu3.line_parameters
            statepp = C2H2_nu1nu3.states.alias()
            statep = C2H2_nu1nu3.states.alias()
            result = conn.execute(
                select(lp.c.a, lp.c.gair, lp.c.dair, lp.c.sw,
                       statepp.c.nu1, statepp.c.J,
                       statep.c.nu1, statep.c.J).\
                join_from(lp, statepp, lp.c.statepp==statepp.c.id).\
                join_from(lp, statep, lp.c.statep==statep.c.id))

            for row in result:
                spp = DiatomState(*row[4:6])
                sp = DiatomState(*row[6:8])
                self.lines[(spp, sp)] = LineParams(*row[:4])

    def _generate_elevels(self, engine: Engine):
        self.elevels = {}
        self.degeneracies = {}
        with engine.begin() as conn:
            st = C2H2_nu1nu3.states
            result = conn.execute(
                select(st.c.energy, st.c.g, st.c.nu1, st.c.J))
            for row in result:
                self.elevels[DiatomState(*row[2:])] = row[0]
                self.degeneracies[DiatomState(*row[2:])] = row[1]

    def _custom_degeneracy(self, state: RotState) -> int:
        """Return quantum state degeneracy.

        For ground state, even is 1 and odd is 3. For excited state, even is 3 and odd is 1."""
        jdeg = 2*state.j+1
        if state.nu == 0:
            if state.j % 2 == 0:
                nucdeg = 1
            else:
                nucdeg = 3
        elif state.nu == 1:
            if state.j % 2 == 0:
                nucdeg = 3
            else:
                nucdeg = 1
        else:
            raise ValueError("Can't calculate degeneracy for '{:s}'".format(state))

        return jdeg*nucdeg

    def tips(self, T: float) -> float:
        """Total internal partition function."""
        return hap.PYTIPS(26, self._iso, T)


class LinearManifold:
    def __init__(self, origin: float, B: float, D: float=0.0, HJ: float=0.0):
        self.origin = origin
        self.B = B
        self.D = D
        self.HJ = HJ

    def energy(self, J: int) -> float:
        JJ1 = J*(J+1)
        return self.origin + self.B*JJ1 - self.D*JJ1**2 + self.HJ*JJ1**3


class SymTopManifold:
    def __init__(self, origin: float, B: float, C: float,
                 D: float=0.0, DJK: float=0.0, DK:float =0.0,
                 HJ: float=0.0, HJK: float=0.0, HKJ: float=0.0,
                 HK: float=0.0):
        self.origin = origin
        self.B = B
        self.C = C              # A or C
        self.D = D
        self.DJK = DJK
        self.DK = DK
        self.HJ = HJ
        self.HJK = HJK
        self.HKJ = HKJ
        self.HK = HK

    def energy(self, J: int, K: int) -> float:
        if K>J:
            raise ValueError("K cannot be larger than J!")
        JJ1 = J*(J+1)
        return self.origin + self.B*JJ1 + (self.C-self.B)*K**2\
            - self.D*JJ1**2 - self.DK*K**4 - self.DJK*JJ1*K**2\
            + self.HJ*JJ1**3 + self.HK*K**6 + self.HJK*JJ1**2*K**2\
            + self.HKJ*JJ1*K**4


class COEffectiveMode(COAlchemyMode):
    def __init__(self, manifolds: Mapping[int, LinearManifold],
                 engine_or_path=None, iso=1):
        COAlchemyMode.__init__(self, engine_or_path, iso)
        self.manifolds = manifolds

    @classmethod
    def from_constants(cls, we: float, Be: float, alphae: float, maxnu: int,
                       wexe: float=0.0, gammae: float=0.0, betae: float=0.0,
                       De: float=0.0) -> "COEffectiveMode":
        nus = list(range(maxnu+1))
        manifolds = []
        for v in nus:
            v12 = v+0.5
            manifolds.append(
                LinearManifold(
                    origin=we*v12 - wexe*v12**2,
                    B=Be - alphae*v12 + gammae*v12**2,
                    D=De + betae*v12
                ))

        return COEffectiveMode(dict(zip(nus, manifolds)))

    def nu(self, pair: Tuple[RotState, RotState]):
        Ep = self.manifolds[pair[1].nu].energy(pair[1].j)
        Epp = self.manifolds[pair[0].nu].energy(pair[0].j)

        return Ep-Epp


class CH3ClEffectiveMode(CH3ClAlchemyMode):
    def __init__(self, manifolds: Mapping[int, SymTopManifold],
                 engine_or_path=None, iso=1):
        CH3ClAlchemyMode.__init__(self, engine_or_path, iso)
        self.manifolds = manifolds

    def nu(self, pair: Tuple[SymTopState, SymTopState]) -> float:
        Ep = self.manifolds[pair[1].nu].energy(pair[1].j, pair[1].k)
        Epp = self.manifolds[pair[0].nu].energy(pair[0].j, pair[0].k)

        return Ep-Epp


# Functions below are mostly for sanity checking my (re)calculations
def tme_sum(pair: Tuple[RotState, RotState],
            mode: AlchemyModeMixin) -> float:
    r"""Unweighted sum of transition matrix elements.

    The sum goes over all degenerate states corresponding to the energy levels
    in ``pair``. The degeneracy may include spatial degeneracy (magnetic
    states), hyperfine degeneracy or other kinds. It is related to the Einstein
    A coefficient for spontaneous emission by:

    .. math::

        \sum_{m_i,m_f,d_if} |\langle \alpha_f J_f m_f d_f|\vec{mu}| J_i m_i d_i \rangle |^2 = \frac{3 A_{fi} \epsilon_0 h c^3 g_f}{16\pi^3\nu^3}

    where :math:`d_x` are labels for degenerate states other than spatial
    states, :math:`g_f=(2J_f+1)D_f` is the total upper state degeneracy
    (statistical weight). For rovibrational transitions, this can be expressed as:

    .. math:

        D_i S(J_f, J_i) \langle \nu_f | \mu | \nu_i \rangle

    See Eq. (32) in M. Simečková, M., Jacquemart, D., Rothman, L. S., Gamache,
    R. R. & Goldman, A. Einstein a-coefficients and statistical weights for
    molecular absorption transitions in the hitran database. Jqsrt 98, 130–155
    (2006)
    and
    Gamache, R. R. & Rothman, L. S. Extension of the HITRAN database to non-LTE
    applications. Journal of quantitative spectroscopy and radiative transfer
    48, 519–525 (1992).
    """
    params, _ = mode.line_params(pair)
    if params is None:
        print(pair)
    A = params.A
    nu = abs(mode.nu(pair))
    if _ == 1:
        statep = pair[1]
    else:
        statep = pair[0]
    deg = mode.degeneracy(statep)

    sum = np.abs(3*A*C.epsilon_0*C.h*C.c**3*deg/(16*np.pi**3*nu**3))

    return sum


def tme_sum_from_line_strength(pair: Tuple[RotState, RotState],
                               mode: AlchemyModeMixin) -> float:
    """The magnitude is not correct."""
    params, _ = mode.line_params(pair)
    if params is None:
        print(pair)
    sw = params.sw
    molid = hap.molname2molid(mode.molecule)
    abundance = hap.h3.ISO[(molid, 1)][2]
    pop_factor = mode.difference_pop(pair, 296.0)
    sum = sw/abundance/pop_factor

    return sum


def rme(pair: Tuple[RotState, RotState],
        mode: AlchemyModeMixin) -> float:
    r"""Calculate reduced matrix element from Einstein A coefficient.

    .. math::

        \sqrt{S(J', J'')} |\langle \nu' | \mu | \nu'' \rangle|
    """
    tme = tme_sum(pair, mode)
    nonrot_deg1 = mode.degeneracy(pair[0])/(2*pair[0].j+1)
    # nonrot_deg2 = mode.degeneracy(pair[1])/(2*pair[1].j+1)
    # print("Degeneracies:", nonrot_deg1, nonrot_deg2)
    rme = np.sqrt(tme/nonrot_deg1)
    if pair[0].j < pair[1].j:
        rme = -rme

    return rme


def hitran_line_strength_rme(pair: Tuple[RotState, RotState],
                             mode: AlchemyModeMixin,
                             unit: str='SI') -> float:
    molid = hap.molname2molid(mode.molecule)
    abundance = hap.h3.ISO[(molid, 1)][2]
    omega = mode.nu(pair)*2*np.pi
    RME = rme(pair, mode)
    nucdeg = mode.degeneracy(pair[0])/(2*pair[0].j+1)
    tme = RME**2*nucdeg
    pop_factor = mode.difference_pop(pair, 296.0)
    fac = omega/(C.epsilon_0*C.c*C.hbar*6.0)
    ls = abundance*tme*pop_factor*fac
    if unit == 'HITRAN':
        ls = ls/C.c/100*1e4

    return ls


def hitran_line_strength(pair: Tuple[RotState, RotState],
                         mode: AlchemyModeMixin,
                         unit: str='SI') -> float:
    """Return HITRAN line strength for `pair`.

    Parameters
    ----------
    pair
        Transition.
    mode
        Molecular vibrational mode.
    unit
        'SI' (default) for m**2 Hz, 'HITRAN' for cm-1 cm**2.
    """
    molid = hap.molname2molid(mode.molecule)
    abundance = hap.h3.ISO[(molid, 1)][2]
    omega = mode.nu(pair)*2*np.pi
    tme = tme_sum(pair, mode)
    pop_factor = mode.difference_pop(pair, 296.0)
    fac = omega/(C.epsilon_0*C.c*C.hbar*6.0)
    ls = abundance*tme*pop_factor*fac
    if unit == 'HITRAN':
        ls = ls/C.c/100*1e4

    return ls


def einstein_from_sw(pair: Tuple[RotState, RotState],
                     mode: AlchemyModeMixin):
    """Calculate Einstein coefficient A21 from HITRAN line strength."""
    molid = hap.molname2molid(mode.molecule)
    nu0 = u.nu2wn(mode.nu(pair))
    abundance = hap.h3.ISO[(molid, 1)][2]
    deg = mode.degeneracy(pair[1])
    c2 = C.h*C.c*100/C.k/296.0

    A = 8*np.pi*100.0*C.c*nu0**2*mode.tips(296.0)*mode.sw(pair)/\
        np.exp(-c2*mode.energy(pair[0]))/(1 - np.exp(-c2*nu0))/\
        abundance/deg

    return A
