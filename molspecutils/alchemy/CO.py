from sqlalchemy import Column, Integer, String, ForeignKey, Float, MetaData, Table
from molspecutils.happier import CO_llq_to_pair, HITRANRow
import molspecutils.alchemy.meta as meta


def state_dicts(row: HITRANRow):
    """Return state dicts for insertion into db."""
    Jpp, Jp = CO_llq_to_pair(row.llq)
    nupp, nup = int(row.glq.strip()), int(row.guq.strip())

    return (dict(J=Jpp, nu=nupp, g=row.gpp, energy=row.elower),
            dict(J=Jp, nu=nup, g=row.gp, energy=row.elower+row.nu))


def line_dict(row: HITRANRow, statepp: int, statep: int):
    """Return line dict for insertion into db."""
    row = row._asdict()
    line = {k: row[k] for k in meta.line_columns}
    line['statepp'] = statepp
    line['statep'] = statep

    return line


metadata = MetaData()

states = Table(
    "states", metadata,
    Column('id', Integer, primary_key=True),
    Column('energy', Float),    # general
    Column('g', Float),
    Column('J', Integer),     # rotational
    Column('nu', Integer),
)

line_parameters = meta.line_parameters.to_metadata(metadata)
