"""Store molecular data extracted from HITRAN in an sqlite3 database.

1. If HITRAN .par/.data file is not present then fetch it with HAPI.
2. If HITRAN file is present but sqlite3 db isn't, then initialize the database
   and populate it. If HITRAN file path is not provided then look for it under
   $XDG_DATA_CACHE/happier/db.
3. If sqlite3 DB is present, then load it whole into memory for faster access.

The main function is :func:`get` which does all these steps if necessary and
returns :class:`sqlalchemy.Engine.engine` instance for the sql databse.
"""
# * Imports
from typing import Union, Sequence
from types import ModuleType
from importlib import import_module
from pathlib import Path
from sqlalchemy import create_engine, select, insert
from sqlalchemy.orm import Session, selectinload
import molspecutils.happier as h
import molspecutils.mirs as mirs
from molspecutils.utils import chunked

StrPath = Union[str, Path]
ParameterNames = ('local_lower_quanta', 'local_upper_quanta', 'global_lower_quanta',
                  'global_upper_quanta', 'nu', 'elower', 'sw', 'a', 'gamma_air',
                  'gamma_self', 'delta_air', 'n_air', 'gp', 'gpp')
VALUES_LIMIT=999

# * Functions
def fetch(hapi_path, molecule_mode, iso):
    try:
        hapi_path.unlink(missing_ok=True)
    except FileNotFoundError:
        pass
    try:
        hapi_path.with_suffix('.header').unlink(missing_ok=True)
    except FileNotFoundError:
        pass

    molname = molecule_mode.split('_')[0]
    h.h3.db_begin(str(hapi_path.parent))
    h.h3.fetch(molecule_mode+'_'+str(iso), h.molname2molid(molname), iso, 0, 20000.0,
               ParameterGroups=['160-char', 'Labels'])


def _append_or_get(insert_states, state):
    try:
        idx = insert_states.index(state)
    except ValueError:
        idx = len(insert_states)
        insert_states.append(state)

    return idx


def _with_id(d, i):
    d['id'] = i

    return d


def convert(sql_path, molecule_mode, iso):
    try:
        sql_path.unlink(missing_ok=True)
    except FileNotFoundError:
        pass

    molmod = import_module('molspecutils.alchemy.' + molecule_mode)
    skip_func = getattr(molmod, 'skip_func', None)

    # prepare data for insertion into DB
    if list(h.h3.LOCAL_TABLE_CACHE.keys()) == ['sampletab']:
        # not initialized or empty directory
        h.h3.db_begin(str(sql_path.parent))
    insert_lines, insert_states = [], []
    for row in zip(*h.h3.getColumns(molecule_mode+'_'+str(iso), ParameterNames)):
        row = h.HITRANRow(*row)
        statepp, statep = molmod.state_dicts(row)
        if skip_func is not None and skip_func(statepp, statep):
            continue
        idpp = _append_or_get(insert_states, statepp)
        idp = _append_or_get(insert_states, statep)

        insert_lines.append(molmod.line_dict(row, idpp+1, idp+1))

    # insert data
    engine = create_engine("sqlite:///" + str(sql_path))
    molmod.metadata.create_all(engine)
    with engine.begin() as conn:
        conn.execute(insert(molmod.states),
                     [_with_id(d, i) for i, d in enumerate(insert_states, start=1)])
        for lines_chunk in chunked(insert_lines, VALUES_LIMIT//len(insert_lines[0])):
            conn.execute(insert(molmod.line_parameters), lines_chunk)
    engine.dispose()


def get(db_path, molecule_mode, iso):
    if db_path is None:
        db_path = h.hitran_cache
    print(f"Using database directory {db_path}")
    Path(db_path).mkdir(parents=True, exist_ok=True)
    sql_path = Path(db_path) / (molecule_mode+'_'+str(iso) + '.sqlite3')
    if not sql_path.is_file():
        print(f"sqlite3 database not found: {sql_path}")
        hapi_path = Path(db_path) / (molecule_mode+'_'+str(iso) + '.data')
        if not hapi_path.is_file():
            print(f"HITRAN data not found, fetching with HAPI")
            fetch(hapi_path, molecule_mode, iso)
        print("Converting from HITRAN to sqlite3")
        convert(sql_path, molecule_mode, iso)

    return create_engine("sqlite:///" + str(sql_path))
