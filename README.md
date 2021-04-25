Modules for (cavity-enhanced) high-resolution molecular spectroscopy.

# Installation
Install the package from our private GitLab repository by executing:

``` sh
pip install molspecutils --extra-index-url https://<token_name>:<token>@gitlab.com/api/v4/projects/26140156/packages/pypi
```

where `<token_name>` and `<token>` are obtained by creating a personal token
with `read_api` scope. See [Creating personal access
tokens](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token)
for details.

To use SQLAlchemy classes which are required for `rotsim2d` module, one needs to
execute the `molspecutils_init` script to download spectroscopic data and fill
the local database. After installing the package with pip, the script should be
present in `$PATH`, so executing the command from shell should be sufficient.

# Development
To make changes and supply patches, it is best to clone the repository:

``` sh
git clone git@gitlab.com:allisonlab/mdcs/spectroscopy.git
```

install [poetry](https://python-poetry.org/docs/#installation) and execute:

``` sh
poetry install
```

This will install all the dependences in a virtual environment and install the
package in development mode.

# TODO
- add more constraints to the database
