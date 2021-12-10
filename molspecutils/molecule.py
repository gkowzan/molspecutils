"""Molecular effective Hamiltonians for calculations of energy levels."""
import logging
from typing import Tuple, Mapping
import abc
from collections import namedtuple
import numpy as np
import scipy.constants as C
from sqlalchemy import select
from sqlalchemy.engine import Engine
import molspecutils.utils as u
import molspecutils.happier as hap
import molspecutils.alchemy.CO as CO
import molspecutils.alchemy.CH3Cl_nu3 as CH3Cl_nu3
from molspecutils.alchemy.convert import get
from typing import Union

class MissingStateError(Exception):
    pass

class MissingLineError(Exception):
    pass

log = logging.getLogger('__name__')


class DiatomState(namedtuple("DiatomState", ["nu", "j"])):
    """Named tuple representing diatom rovib state with `nu`, `j` quantum numbers."""
    __slots__ = ()

    @classmethod
    def from_symtop(cls, symtop):
        """Drop `k` from :class:`SymTopState`."""
        return cls(nu=symtop.nu, j=symtop.j)

    @property
    def name(self):
        return "{:d},{:d}".format(self.nu, self.j)

    def __repr__(self):
        return "DiatomState(nu={:d}, j={:d})".format(self.nu, self.j)

    def __eq__(self, o):
        if not isinstance(o, DiatomState):
            return NotImplemented
        return self.nu == o.nu and self.j == o.j

    def __hash__(self):
        return hash((self.nu, self.j))


class SymTopState(namedtuple("SymTopState", ["nu", "j", "k"])):
    """Named tuple representing symmetric top rovib state with `nu`, `j`, `k`
    quantum numbers."""
    __slots__ = ()

    @property
    def name(self):
        return "{:d},{:d},{:d}".format(self.nu, self.j, self.k)

    def __repr__(self):
        return "SymTopState(nu={:d}, j={:d}, k={:d})".format(self.nu, self.j, self.k)

    def __eq__(self, o):
        if not isinstance(o, SymTopState):
            return NotImplemented
        return self.nu == o.nu and self.k == o.k and self.j == o.j

    def __hash__(self):
        return hash((self.nu, self.j, self.k))


RotState = Union[DiatomState, SymTopState]


class VibrationalMode(abc.ABC):
    """Interface of a class representing a molecular vibrational mode."""
    @abc.abstractmethod
    def gamma(self, pair: Tuple[RotState, RotState]):
        """Pressure broadening coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def delta(self, pair: Tuple[RotState, RotState]):
        """Pressure shift coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def mu(self, pair: Tuple[RotState, RotState]):
        """Reduced matrix element for dipole transition between `pair` states."""

    @abc.abstractmethod
    def nu(self, pair: Tuple[RotState, RotState]):
        """Frequency of molecular coherence `pair`."""

    @abc.abstractmethod
    def equilibrium_pop(self, state: RotState, T: float):
        """Fractional population of `state` at temperature `T` in thermal
        equilibrium."""


class AlchemyModeMixin:
    """Molecular vibrational mode backed by SQLAlchemy sqlite database.

    This is a mix-in."""
    def line_params(self, pair: Tuple[RotState, RotState]):
        result = self._line_params(pair)
        if result:
            return (result, 1)
        result = self._line_params(pair[::-1])
        if result is None:
            log.debug("Missing line parameters for: %r", pair)
            result = dict(zip(['A', 'gamma', 'delta', 'sw'], [0]*3))
        return (result, -1)

    def gamma(self, pair: Tuple[RotState, RotState]):
        params, _ = self.line_params(pair)
        if params is None or params['gamma'] == 0:
            gam = self._fake_gamma(pair)
        else:
            gam = params['gamma']

        return u.wn2nu(gam)

    def sw(self, pair: Tuple[RotState, RotState]):
        params, _ = self.line_params(pair)
        sw = params['sw']

        return sw/1e4           # m**2/molecule

    def delta(self, pair: Tuple[RotState, RotState]):
        params, _ = self.line_params(pair)
        delt = params['delta']

        return u.wn2nu(delt)

    def nu(self, pair: Tuple[RotState, RotState]):
        try:
            return u.wn2nu(self.elevels[pair[1]]-self.elevels[pair[0]])
        except KeyError as e:
            raise MissingStateError("Energy for state `{!s}` is not available".format(
                e.args[0]
            ))

    def _fake_gamma(self, pair: Tuple[RotState, RotState]):
        for kpp, kp in self.lines.keys():
            if pair[0]==kp or pair[0]==kpp or pair[1]==kp or pair[1]==kpp:
                return self.lines[(kpp, kp)]['gamma']

    def _line_params(self, pair: Tuple[RotState, RotState]):
        return self.lines.get(pair)


class CH3ClAlchemyMode(AlchemyModeMixin, VibrationalMode):
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
        engine_or_path
            Path to directory with HAPI or sqlite3 database or
            :class:`sqlachemy.Engine` instance of opened sqlite3 database.
        iso
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
        with engine.begin() as conn:
            st = CH3Cl_nu3.states
            result = conn.execute(
                select(st.c.energy, st.c.nu3, st.c.J, st.c.K))
            for row in result:
                self.elevels[SymTopState(*row[1:])] = row[0]

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
                self.lines[(spp, sp)] = dict(
                    zip(['A', 'gamma', 'delta', 'sw'], row[:4]))

    def mu(self, pair: Tuple[RotState, RotState]):
        r"""Reduced matrix element for `pair[0]` to `pair[1]` transition.

        Obtained from HITRAN's Einsten A-coefficient:

        .. math::

            |\langle \nu';J'\|\mu^{(1)}\|\nu'';J''\rangle|^{2} = A_{\nu'J'\to\nu''J''}\frac{\epsilon_{0}hc^{3}(2J'+1)}{16\pi^{3}\nu^{3}_{\nu'J',\nu''J''}}
        """
        params, _ = self.line_params(pair)
        if params is None:
            print(pair)
        A = params['A']
        nu = abs(self.nu(pair))
        if _ == 1:
            j = pair[1].j
        elif _ == -1:
            j = pair[0].j
        # rmu = np.sqrt(A*(2*pair[0].j+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*nu)**3)
        # Eq. (5.9) from rotsim2d_roadmap
        rmu = np.sqrt(3*A*C.epsilon_0*C.h*C.c**3*(2*j+1)/(16*np.pi**3*nu**3))
        if pair[0].j < pair[1].j:
            rmu = -rmu

        return rmu

    def equilibrium_pop(self, state: RotState, T: float):
        kt = u.joule2wn(C.k*T)
        e = self.elevels[state]
        gnuc = 16
        k3 = 2 if state.k > 0 and state.k % 3 == 0 else 1
        fudge_factor = 0.5

        # return np.sqrt(kfac*(2*state.j+1))*np.exp(-e/kt)/hap.PYTIPS(24, 1, T)
        return gnuc*k3*(2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(24, 1, T)*fudge_factor


class COAlchemyMode(AlchemyModeMixin, VibrationalMode):
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
                self.lines[(spp, sp)] = dict(
                    zip(['A', 'gamma', 'delta', 'sw'], row[:4]))

    def _generate_elevels(self, engine: Engine):
        self.elevels = {}
        with engine.begin() as conn:
            st = CO.states
            result = conn.execute(
                select(st.c.energy, st.c.nu, st.c.J))
            for row in result:
                self.elevels[DiatomState(*row[1:])] = row[0]

    def mu(self, pair: Tuple[RotState, RotState]):
        r"""Reduced matrix element for the `pair` transitions.

        Obtained from HITRAN's Einsten A-coefficient:

        .. math::

            |\langle \nu';J'\|\mu^{(1)}\|\nu'';J''\rangle|^{2} = A_{\nu'J'\to\nu''J''}\frac{\epsilon_{0}hc^{3}(2J'+1)}{16\pi^{3}\nu^{3}_{\nu'J',\nu''J''}}
        """
        params, _ = self.line_params(pair)
        if params is None:
            print(pair)
        A = params['A']
        nu = abs(self.nu(pair))
        if _ == 1:
            j = pair[1].j
        elif _ == -1:
            j = pair[0].j
        # rmu = np.sqrt(A*(2*pair[0].j+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*nu)**3)
        # Eq. (5.9) from rotsim2d_roadmap
        rmu = np.sqrt(3*A*C.epsilon_0*C.h*C.c**3*(2*j+1)/(16*np.pi**3*nu**3))

        return rmu

    def equilibrium_pop(self, state: RotState, T: float):
        r"""Fractional population of spatial sublevel of `state` in thermal equilibrium,
        :math:`e^{-E_{\nu,j}/kT}/Q`.
        """
        kt = u.joule2wn(C.k*T)
        e = self.elevels[state]

        # return np.sqrt(2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(5, 1, T)
        return (2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(5, self._iso, T)


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

    def nu(self, pair: Tuple[RotState]):
        Ep = self.manifolds[pair[1].nu].energy(pair[1].j)
        Epp = self.manifolds[pair[0].nu].energy(pair[0].j)

        return Ep-Epp
