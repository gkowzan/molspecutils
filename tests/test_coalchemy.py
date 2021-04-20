from pathlib import Path
from sqlalchemy.orm import Session
from sqlalchemy import create_engine
import spectroscopy.happier as h
from spectroscopy.molecule import COAlchemyMode, DiatomState

sql_path = Path(h.hitran_cache) / 'CO.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))

co_mode = COAlchemyMode(engine)
statepp = DiatomState(nu=0, j=1)
statep = DiatomState(nu=1, j=2)

params = co_mode.line_params((statepp, statep))
