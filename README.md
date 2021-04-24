Modules for (cavity-enhanced) high-resolution molecular spectroscopy.

# Installation
Install the package by downloading or cloning the repo and calling the following inside main repository directory (containing `setup.py`):

``` sh
python -m pip install -e .
```

or by installing directly from the repo with pip

``` sh
python -m pip install git+ssh://git@gitlab.com/allisonlab/mdcs/molspecutils.git@master
```

To use SQLAlchemy classes which are required for `rotsim2d` module, one needs to
execute the `molspecutils_init` script to download spectroscopic data and fill
the local database. After installing the package with pip, the script should be
present in `$PATH`, so executing the command from shell should be sufficient.

# TODO
- add more constraints to the database
