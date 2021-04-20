Modules for (cavity-enhanced) linear spectroscopy.

# Installation
Install the package by downloading or cloning the repo and calling the following inside main repository directory (containing `setup.py`):

``` sh
python -m pip install -e .
```

or by installing directly from the repo with pip

``` sh
python -m pip git+ssh://git@gitlab.com:allisonlab/mdcs/spectroscopy.git@master
```

To use SQLAlchemy classes which are required for `rotsim2d` molecule, one needs
to execute the `spectroscopy/scripts/initialize_dbs.py` script to download
spectroscopic data and fill the local database.

# TODO
- add more constraints to the database
