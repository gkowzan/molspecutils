"""Extract SQL table of LineParameter's into a dict.

dict will be indexed by SymTopState instances."""
# * Imports
from pathlib import Path
from sqlalchemy import create_engine, select
from sqlalchemy.orm import aliased, Session, selectinload
import molspecutils.happier as h
from molspecutils.alchemy import CH3Cl
from molspecutils.molecule import SymTopState

# * DBs
input_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
input_engine = create_engine("sqlite:///" + str(input_path))
input_session = Session(bind=input_engine)

# * Query DB
rovpp, rovp = aliased(CH3Cl.RovibState), aliased(CH3Cl.RovibState)
nupp, nup = aliased(CH3Cl.VibState), aliased(CH3Cl.VibState)
input_result = input_session.execute(
    select(CH3Cl.LineParameters).options(selectinload('*'))\
    .join(CH3Cl.LineParameters.transition)\
    .join(rovpp, rovpp.id==CH3Cl.TransitionPair.statepp_id)\
    .join(rovp, rovp.id==CH3Cl.TransitionPair.statep_id)\
    .join(nupp, nupp.id==rovpp.nu_id).join(nup, nup.id==rovp.nu_id).
    where(nupp.nu1==0, nupp.nu2==0, nupp.nu4==0, nupp.nu5==0, nupp.nu6==0,
          nup.nu1==0, nup.nu2==0, nup.nu4==0, nup.nu5==0, nup.nu6==0)).scalars()

# * Convert to dict
elevel_dict, trans_dict = {}, {}
for lp in input_result:
    alch_statepp = lp.transition.statepp
    statepp = SymTopState(nu=alch_statepp.nu.nu3, j=alch_statepp.j.j,
                          k=alch_statepp.j.k)
    alch_statep = lp.transition.statep
    statep = SymTopState(nu=alch_statep.nu.nu3, j=alch_statep.j.j,
                         k=alch_statep.j.k)

    trans_dict[(statepp, statep)] = dict(A=lp.A, gamma=lp.gamma_air, delta=lp.delta_air)
    elevel_dict[statepp] = alch_statepp.energy
    elevel_dict[statep] = alch_statep.energy

