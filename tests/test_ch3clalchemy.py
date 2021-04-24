from pathlib import Path
from sqlalchemy.orm import Session, selectinload, aliased
from sqlalchemy import create_engine, select
import molspecutils.happier as h
from molspecutils.molecule import CH3ClAlchemyMode, SymTopState
from molspecutils.alchemy import CH3Cl

sql_path = Path(h.hitran_cache) / 'CH3Cl.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))

ch3cl_mode = CH3ClAlchemyMode(engine)
statepp = SymTopState(nu=0, j=3, k=3)
statep = SymTopState(nu=1, j=3, k=3)

params = ch3cl_mode.line_params((statepp, statep))

with ch3cl_mode.sess() as session:
    print(ch3cl_mode._state_energy(session, SymTopState(nu=0, j=3, k=3)))

with ch3cl_mode.sess() as session:
    state = SymTopState(nu=1, j=3, k=3)
    result = session.execute(
        select(CH3Cl.RovibState).join(CH3Cl.VibState).join(CH3Cl.RotState).\
            where(CH3Cl.VibState.nu1==0, CH3Cl.VibState.nu2==0, CH3Cl.VibState.nu3==state.nu,
                  CH3Cl.VibState.nu4==0, CH3Cl.VibState.nu5==0, CH3Cl.VibState.nu6==0,
                  CH3Cl.RotState.j==state.j, CH3Cl.RotState.k==state.k, CH3Cl.RotState.l==0,
                  CH3Cl.RotState.f==0.0)
    ).scalars().first()
    print(result)

pair = (statepp, statep)
with ch3cl_mode.sess() as session:
    rovpp, rovp = aliased(CH3Cl.RovibState), aliased(CH3Cl.RovibState)
    subjpp = select(CH3Cl.RotState.id).filter_by(
        j=pair[0].j, k=pair[0].k, l=0, f=0.0
    ).subquery()
    subnupp = select(CH3Cl.VibState.id).filter_by(
        nu1=0, nu2=0, nu3=pair[0].nu, nu4=0, nu5=0, nu6=0
    ).subquery()
    subjp = select(CH3Cl.RotState.id).filter_by(
        j=pair[1].j, k=pair[1].k, l=0, f=-0.0
    ).subquery()
    subnup = select(CH3Cl.VibState.id).filter_by(
        nu1=0, nu2=0, nu3=pair[1].nu, nu4=0, nu5=0, nu6=0
    ).subquery()
    result = session.execute(
        select(CH3Cl.LineParameters).options(selectinload('*')).\
        join(CH3Cl.LineParameters.transition).\
        join(rovpp, CH3Cl.TransitionPair.statepp).join(subnupp, subnupp.c.id == rovpp.nu_id).\
        join(subjpp, subjpp.c.id == rovpp.j_id).\
        join(rovp, CH3Cl.TransitionPair.statep).join(subnup, subnup.c.id == rovp.nu_id).\
        join(subjp, subjp.c.id == rovp.j_id)).scalars().all()
