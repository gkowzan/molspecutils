# * Imports
from typing import Union
from importlib import import_module
from pathlib import Path
from sqlalchemy import create_engine, select
import sqlalchemy.exc as exc
from sqlalchemy.orm import Session, selectinload
import spectroscopy.happier as h

ParameterNames = ('local_lower_quanta', 'local_upper_quanta', 'global_lower_quanta',
                  'global_upper_quanta', 'nu', 'elower', 'sw', 'a', 'gamma_air',
                  'gamma_self', 'delta_air', 'n_air', 'gp', 'gpp')

# * Functions
def only_jnu(cls, kwargs):
    if cls.__name__ == 'RovibState':
        return dict(j=kwargs['j'], nu=kwargs['nu'])
    else:
        return kwargs

def add_or_create(session, cls, kwargs):
    obj = session.execute(select(cls).filter_by(**only_jnu(cls, kwargs))).scalars().one_or_none()
    if obj is None:
        obj = cls(**kwargs)
        session.add(obj)
    return obj

def convert(mol_name: str, cache: Union[str, Path], overwrite_hitran=False, overwrite_alchemy=False):
    """Save HITRAN molecule data into sqlite3 DB."""
    hapi_path = Path(cache) / (mol_name + '.data')
    sql_path = Path(cache) / (mol_name + '.sqlite3')

    mol_id = h.molname2molid(mol_name)
    h.h3.db_begin(cache)
    if overwrite_hitran:
        hapi_path.unlink(missing_ok=True)
        hapi_path.with_suffix('.header').unlink(missing_ok=True)
        h.h3.fetch(mol_name, mol_id, 1, 0, 20000.0, ParameterGroups=['160-char', 'Labels'])
    if overwrite_alchemy:
        sql_path.unlink(missing_ok=True)

    engine = create_engine("sqlite:///" + str(sql_path))
    molmod = import_module('spectroscopy.alchemy.'+mol_name)
    molmod.Base.metadata.create_all(engine)
    
    with Session(engine) as session:
        for row in zip(*h.h3.getColumns(mol_name, ParameterNames)):
            llq, luq, glq, guq, nu, el, sw, a, gair, gself, dair, nair, gp, gpp = row
            rotppd, rotpd = molmod.local_state_convert(llq, luq)
            vibppd, vibpd = molmod.global_state_convert(glq, guq)

            rotpp_state = add_or_create(session, molmod.RotState, rotppd)
            rotp_state = add_or_create(session, molmod.RotState, rotpd)

            vibpp_state = add_or_create(session, molmod.VibState, vibppd)
            vibp_state = add_or_create(session, molmod.VibState, vibpd)

            rovibpp_state = add_or_create(
                session, molmod.RovibState, dict(j=rotpp_state, nu=vibpp_state, energy=el, g=gpp))
            rovibp_state = add_or_create(
                session, molmod.RovibState, dict(j=rotp_state, nu=vibp_state, energy=el+nu, g=gp))

            transition_pair = add_or_create(
                session, molmod.TransitionPair, dict(statepp=rovibpp_state, statep=rovibp_state))

            line_params = molmod.LineParameters(
                nu=nu, sw=sw, A=a, gamma_air=gair, gamma_self=gself, n_air=nair, delta_air=dair,
                transition=transition_pair)
            # HITRAN database has duplicate rows, great...
            result = session.execute(
                    select(molmod.LineParameters).options(selectinload('*'))\
                    .where(molmod.LineParameters.transition==transition_pair)
            ).scalars().one_or_none()
            if result is not None:
                print('Skipping duplicate:')
                print(line_params)
                print('HITRAN row:', row)
                print('state fragment', ''.join([guq, glq, luq, llq]))
                print('of')
                print(result, '\n')
                continue
            session.add(line_params)
        session.commit()


if __name__ == '__main__':
    from sqlalchemy.orm import aliased
    from sqlalchemy import and_
    import spectroscopy.alchemy.CO as CO
    nupp, nup = aliased(CO.VibState), aliased(CO.VibState)
    jpp, jp = aliased(CO.RotState), aliased(CO.RotState)
    rovpp, rovp = aliased(CO.RovibState), aliased(CO.RovibState)
    sql_path = Path(h.hitran_cache) / 'CO.sqlite3'
    engine = create_engine("sqlite:///" + str(sql_path))
    with Session(engine) as session:
        result = session.execute(
            select(CO.LineParameters.gamma_air).join(CO.LineParameters.transition).\
            join(rovpp, CO.TransitionPair.statepp).join(nupp, rovpp.nu).\
            join(jpp, rovpp.j).\
            join(rovp, CO.TransitionPair.statep).join(nup, rovp.nu).\
            join(jp, rovp.j).filter(and_(jpp==CO.RotState(j=1), nupp.nu==0, jp.j==2, nup.nu==1))
        ).scalars()
        print(result.one())

    with Session(engine) as session:
        subjpp = select(CO.RotState.id).filter_by(j=1).scalar_subquery()
        subnupp = select(CO.VibState.id).filter_by(nu=0).scalar_subquery()
        subjp = select(CO.RotState.id).filter_by(j=2).scalar_subquery()
        subnup = select(CO.VibState.id).filter_by(nu=1).scalar_subquery()
        result = session.execute(
            select(CO.LineParameters.gamma_air).join(CO.LineParameters.transition).\
            join(rovpp, CO.TransitionPair.statepp).join(nupp, rovpp.nu).\
            join(jpp, rovpp.j).\
            join(rovp, CO.TransitionPair.statep).join(nup, rovp.nu).\
            join(jp, rovp.j).\
            where(jpp.id==subjpp, nupp.id==subnupp, jp.id==subjp, nup.id==subnupp)
        ).scalars()
        print(result.all())

    with Session(engine) as session:
        print(session.execute(
            select(CO.RovibState.energy).join(CO.VibState).join(CO.RotState).\
            where(CO.VibState.nu==0, CO.RotState.j==1)
        ).scalars().all())
