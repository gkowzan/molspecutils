"""Molecular effective Hamiltonians for calculations of energy levels."""
from typing import Tuple
import abc
from collections import namedtuple
import numpy as np
import scipy.constants as C
from sqlalchemy.orm import aliased, selectinload, Session
from sqlalchemy import select
import shed.units as u
import spectroscopy.happier as hap
import spectroscopy.alchemy.CO as CO
import spectroscopy.alchemy.CH3Cl as CH3Cl

class RotState(abc.ABC):
    pass


@RotState.register
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


@RotState.register
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


class VibrationalMode(abc.ABC):
    """Interface of a class representing a molecular vibrational mode."""
    @abc.abstractmethod
    def gamma(self, pair: Tuple[RotState]):
        """Pressure broadening coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def delta(self, pair: Tuple[RotState]):
        """Pressure shift coefficient for `pair` molecular coherence."""

    @abc.abstractmethod
    def mu(self, pair: Tuple[RotState]):
        """Reduced matrix element for dipole transition between `pair` states."""

    @abc.abstractmethod
    def nu(self, pair: Tuple[RotState]):
        """Frequency of molecular coherence `pair`."""

    @abc.abstractmethod
    def equilibrium_pop(self, state: RotState, T: float):
        """Fractional population of `state` at temperature `T` in thermal
        equilibrium."""


class AlchemyModeMixin:
    """Molecular vibrational mode backed by SQLAlchemy sqlite database.

    This is a mix-in."""
    def line_params(self, pair: Tuple[RotState]):
        result = self._line_params(pair)
        if result:
            return (result, 1)
        return (self._line_params(pair[::-1]), -1)

    def gamma(self, pair: Tuple[RotState]):
        params, _ = self.line_params(pair)
        if params is None:
            gam = self._fake_gamma(pair)
        else:
            gam = params['gamma']

        return u.wn2nu(gam)

    def delta(self, pair: Tuple[RotState]):
        params, _ = self.line_params(pair)
        delt = params['delta']

        return u.wn2nu(delt)

    def nu(self, pair: Tuple[RotState]):
        return u.wn2nu(self.elevels[pair[1]]-self.elevels[pair[0]])

    def _fake_gamma(self, pair: Tuple[RotState]):
        for kpp, kp in self.lines.keys():
            if pair[0]==kp or pair[0]==kpp or pair[1]==kp or pair[1]==kpp:
                return self.lines[(kpp, kp)]['gamma']

    def _line_params(self, pair: Tuple[RotState]):
        return self.lines.get(pair)


class CH3ClAlchemyMode(AlchemyModeMixin, VibrationalMode):
    def __init__(self, engine):
        self._generate_dicts(engine)

    def _generate_dicts(self, engine):
        self.elevels, self.lines = {}, {}
        rovpp, rovp = aliased(CH3Cl.RovibState), aliased(CH3Cl.RovibState)
        nupp, nup = aliased(CH3Cl.VibState), aliased(CH3Cl.VibState)

        with Session(bind=engine) as session:
            input_result = session.execute(
                select(CH3Cl.LineParameters).options(selectinload('*'))\
                .join(CH3Cl.LineParameters.transition)\
                .join(rovpp, rovpp.id==CH3Cl.TransitionPair.statepp_id)\
                .join(rovp, rovp.id==CH3Cl.TransitionPair.statep_id)\
                .join(nupp, nupp.id==rovpp.nu_id).join(nup, nup.id==rovp.nu_id).
                where(nupp.nu1==0, nupp.nu2==0, nupp.nu4==0, nupp.nu5==0, nupp.nu6==0,
                      nup.nu1==0, nup.nu2==0, nup.nu4==0, nup.nu5==0, nup.nu6==0)).scalars()
            for lp in input_result:
                alch_statepp = lp.transition.statepp
                alch_statep = lp.transition.statep
                if alch_statepp.j.f != 0.0 or alch_statep.j.f != 0.0:
                    continue
                statepp = SymTopState(nu=alch_statepp.nu.nu3, j=alch_statepp.j.j,
                                      k=alch_statepp.j.k)
                statep = SymTopState(nu=alch_statep.nu.nu3, j=alch_statep.j.j,
                                     k=alch_statep.j.k)

                self.lines[(statepp, statep)] = dict(A=lp.A, gamma=lp.gamma_air, delta=lp.delta_air)
                self.elevels[statepp] = alch_statepp.energy
                self.elevels[statep] = alch_statep.energy

    def mu(self, pair: Tuple[RotState]):
        params, _ = self.line_params(pair)
        if params is None:
            print(pair)
        A = params['A']
        nu = abs(self.nu(pair))
        rmu = np.sqrt(A*(2*pair[0].j+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*nu)**3)

        return rmu

    def equilibrium_pop(self, state: RotState, T: float):
        kt = u.joule2wn(C.k*T)
        e = self.elevels[state]
        kfac = 1 if state.j==0 else 2

        # return np.sqrt(kfac*(2*state.j+1))*np.exp(-e/kt)/hap.PYTIPS(24, 1, T)
        return kfac*(2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(24, 1, T)


class COAlchemyMode(AlchemyModeMixin, VibrationalMode):
    def __init__(self, engine):
        self._generate_dicts(engine)

    def _generate_dicts(self, engine):
        self.elevels, self.lines = {}, {}

        with Session(bind=engine) as session:
            input_result = session.execute(
                select(CO.LineParameters).options(selectinload('*'))).scalars()
            for lp in input_result:
                alch_statepp = lp.transition.statepp
                alch_statep = lp.transition.statep
                statepp = DiatomState(nu=alch_statepp.nu.nu, j=alch_statepp.j.j)
                statep = DiatomState(nu=alch_statep.nu.nu, j=alch_statep.j.j)

                self.lines[(statepp, statep)] = dict(A=lp.A, gamma=lp.gamma_air, delta=lp.delta_air)
                self.elevels[statepp] = alch_statepp.energy
                self.elevels[statep] = alch_statep.energy

    def mu(self, pair: Tuple[RotState]):
        params, _ = self.line_params(pair)
        if params is None:
            print(pair)
        A = params['A']
        nu = abs(self.nu(pair))
        # rmu = np.sqrt(A*(2*pair[0].j+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*nu)**3)
        # Eq. (5.9) from rotsim2d_roadmap
        rmu = np.sqrt(A*C.epsilon_0*C.h*C.c**3*(2*pair[1].j+1)/(16*np.pi**3*nu**3)) 

        return rmu

    def equilibrium_pop(self, state: RotState, T: float):
        kt = u.joule2wn(C.k*T)
        e = self.elevels[state]

        # return np.sqrt(2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(5, 1, T)
        return (2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(5, 1, T)


# class COAlchemyMode(AlchemyModeMixin, VibrationalMode):
#     def _line_params(self, pair: Tuple[RotState]):
#         """Returns a tuple of (nu, A, gamma_air, delta_air).

#         Don't return the whole LineParameters object for performance reasons."""
#         nupp, nup = aliased(CO.VibState), aliased(CO.VibState)
#         jpp, jp = aliased(CO.RotState), aliased(CO.RotState)
#         rovpp, rovp = aliased(CO.RovibState), aliased(CO.RovibState)
        
#         subjpp = select(CO.RotState.id).filter_by(j=pair[0].j).scalar_subquery()
#         subnupp = select(CO.VibState.id).filter_by(nu=pair[0].nu).scalar_subquery()
#         subjp = select(CO.RotState.id).filter_by(j=pair[1].j).scalar_subquery()
#         subnup = select(CO.VibState.id).filter_by(nu=pair[1].nu).scalar_subquery()

#         result = self.sess.execute(
#             select(CO.LineParameters.nu, CO.LineParameters.A, CO.LineParameters.gamma_air,
#                    CO.LineParameters.delta_air).options(selectinload('*')).\
#             join(CO.LineParameters.transition).\
#             join(rovpp, CO.TransitionPair.statepp).join(nupp, rovpp.nu).\
#             join(jpp, rovpp.j).\
#             join(rovp, CO.TransitionPair.statep).join(nup, rovp.nu).\
#             join(jp, rovp.j).\
#             where(jpp.id==subjpp, nupp.id==subnupp, jp.id==subjp, nup.id==subnup)
#         ).one_or_none()

#         return result
        
#     def _fake_gamma(self, pair: Tuple[RotState]):
#         return self.sess.execute(
#             select(CO.LineParameters.gamma_air).join(CO.LineParameters.transition).\
#             join(CO.TransitionPair.statepp).join(CO.RovibState.j).join(CO.RovibState.nu).\
#             where(CO.RotState.j==pair[0].j, CO.VibState.nu==pair[0].nu)
#         ).scalars().first()
    
#     def _state_energy(self, session, state):
#         return session.execute(
#             select(CO.RovibState.energy).join(CO.VibState).join(CO.RotState).\
#             where(CO.VibState.nu==state.nu, CO.RotState.j==state.j)
#         ).scalars().one()

#     # def nu(self, pair: Tuple[RotState]):
#     #     params, sign = self.line_params(pair)
#     #     nu = params.nu
#     #     if self.frequency:
#     #         nu = u.wn2nu(nu)

#     #     return nu

#     def mu(self, pair: Tuple[RotState]):
#         params, sign = self.line_params(pair)
#         nu, A = params[:2]
#         rmu = np.sqrt(A*(2*pair[0].j+1)*C.c**3*C.hbar*np.pi*C.epsilon_0*3/(2*np.pi*u.wn2nu(nu))**3)

#         return rmu

#     def equilibrium_pop(self, state: RotState, T: float):
#         kt = u.joule2wn(C.k*T)
#         e = self._state_energy(self.sess, state)

#         return np.sqrt(2*state.j+1)*np.exp(-e/kt)/hap.PYTIPS(5, 1, T)
    
# class Manifold:
#     def __init__(self, origin: float, B: float, C: float,
#                  D: float=0.0, DJK: float=0.0, DK:float =0.0,
#                  HJ: float=0.0, HJK: float=0.0, HKJ: float=0.0,
#                  HK: float=0.0):
#         self.origin = origin
#         self.B = B
#         self.C = C              # A or C
#         self.D = D
#         self.DJK = DJK
#         self.DK = DK
#         self.HJ = HJ
#         self.HJK = HJK
#         self.HKJ = HKJ
#         self.HK = HK

#     def energy(self, J: int, K: int):
#         if K>J:
#             return ValueError("K cannot be larger than J!")
#         JJ1 = J*(J+1)
#         return self.origin + self.B*JJ1 + (self.C-self.B)*K**2\
#             - self.DJ*JJ1**2 - self.DK*K**4 - self.DJK*JJ1*K**2\
#             + self.HJ*JJ1**3 + self.HK*K**6 + self.HJK*JJ**2*K**2 + self.HKJ*JJ1*K**4
