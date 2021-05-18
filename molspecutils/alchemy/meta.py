"""Take data from HAPI and put it in an sqlite database with sqlalchemy."""
from pathlib import Path
from pathlib import Path
import appdirs
from sqlalchemy import Column, Integer, String, ForeignKey, Float, UniqueConstraint, MetaData, Table

dirs = appdirs.AppDirs('happier', 'gkowzan')
hitran_cache = str(Path(dirs.user_cache_dir) / 'db')
_metadata = MetaData()

line_parameters = Table(
    "line_parameters", _metadata,
    Column('id', Integer, primary_key=True),
    Column('nu', Float),
    Column('sw', Float),
    Column('a', Float),
    Column('gair', Float),
    Column('gself', Float),
    Column('dair', Float),
    Column('nair', Float),
    Column('statepp', Integer, ForeignKey('states.id'), nullable=False),
    Column('statep', Integer, ForeignKey('states.id'), nullable=False)
)
line_columns = ['nu', 'sw', 'a', 'gair', 'gself', 'dair', 'nair']
