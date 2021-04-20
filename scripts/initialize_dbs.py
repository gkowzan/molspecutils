"""Initialize sqlalchemy databases for CO and CH3Cl."""
from pathlib import Path
from spectroscopy.alchemy.convert import convert
from spectroscopy.alchemy.meta import hitran_cache

Path(hitran_cache).mkdir(parents=True, exist_ok=True)

print("Fetching and converting CO data...")
convert('CO', hitran_cache, True, True)

print("Fetching and converting CH3Cl data (this might take a while)...")
convert('CH3Cl', hitran_cache, True, True)
