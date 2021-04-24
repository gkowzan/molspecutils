"""Make a smaller database with only nu3 transitions."""
# * Imports
from pathlib import Path
from sqlalchemy import create_engine, select
from sqlalchemy.orm import aliased, Session, selectinload
import molspecutils.happier as h
from molspecutils.alchemy import CH3Cl

# * DBs
input_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
input_engine = create_engine("sqlite:///" + str(input_path))
input_session = Session(bind=input_engine)

output_path = Path(h.hitran_cache) / 'CH3Cl_nu3.sqlite3'
output_engine = create_engine("sqlite:///" + str(output_path))
output_session = Session(bind=output_engine)
CH3Cl.Base.metadata.create_all(output_engine)

# * Query input DB
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
for lp in input_result:
    output_session.merge(lp)
    
output_session.commit()
