# * Imports and DB
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from sqlalchemy import create_engine, select
from spectroscopy.molecule import CH3ClAlchemyMode, SymTopState
from spectroscopy.alchemy.meta import hitran_cache
plt.ion()

sql_path = Path(hitran_cache) / 'CH3Cl.sqlite3'
engine = create_engine("sqlite:///" + str(sql_path))

ch3cl_mode = CH3ClAlchemyMode(engine)

# * Plot intensities
xs = list(range(2, 16))
yr = np.array([ch3cl_mode.mu((SymTopState(0, x, 1), SymTopState(1, x+1, 1))) for x in xs])
yq = np.array([ch3cl_mode.mu((SymTopState(0, x, 1), SymTopState(1, x, 1))) for x in xs])
yp = np.array([ch3cl_mode.mu((SymTopState(0, x, 1), SymTopState(1, x-1, 1))) for x in xs])

fig, ax = plt.subplots()
ax.plot(xs, yr, label='R')
ax.plot(xs, yq, label='Q')
ax.plot(xs, yp, label='P')
ax.legend(loc='best')
