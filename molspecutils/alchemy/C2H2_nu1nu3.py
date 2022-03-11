import molspecutils.alchemy.meta as meta
from attrs import asdict
from molspecutils.happier import C2H2_llq_to_pair, C2H2GlobalState, HITRANRow
from sqlalchemy import (Column, Float, ForeignKey, Integer, MetaData, String,
                        Table)


def state_dicts(row: HITRANRow):
    """Return state dicts for insertion into db."""
    jpp, jp = C2H2_llq_to_pair(row.llq)

    statepp = asdict(C2H2GlobalState.from_str(row.glq))
    statepp.update({'g': row.gpp,
                    'energy': row.elower,
                    'J': jpp})

    statep = asdict(C2H2GlobalState.from_str(row.guq))
    statep.update({'g': row.gp,
                   'energy': row.elower+row.nu,
                   'J': jp})

    return statepp, statep


def skip_func(statepp: dict, statep: dict) -> bool:
    pp = [statepp['nu'+str(nu)]==0 for nu in (1, 2, 3, 4, 5)]
    p = [statep['nu'+str(nu)]==0 for nu in (2, 4, 5)]
    p1 = [statep['nu1']==1, statep['nu3']==1]

    return not all(pp + p + p1)

metadata = MetaData()

states = Table(
    "states", metadata,
    Column('id', Integer, primary_key=True),
    Column('energy', Float),    # general
    Column('g', Float),
    Column('J', Integer),     # rotational

    Column('nu1', Integer),   # vibrational
    Column('nu2', Integer),
    Column('nu3', Integer),
    Column('nu4', Integer),
    Column('nu5', Integer),
    Column('l4', Integer),
    Column('l5', Integer),

    Column('plus', String(1)),
    Column('S', String(1))
)

line_parameters = meta.line_parameters.to_metadata(metadata)
