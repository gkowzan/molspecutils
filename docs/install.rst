Installation
============

**Dependencies**

- numpy, scipy, appdirs, SQLAlchemy
- `knickknacks <https://gitlab.com/allisonlab/mdcs/shed>`_

**Install**

.. highlight:: sh

Install the package from our private GitLab repository by executing::

  pip install rotsim2d --extra-index-url https://<token_name>:<token>@gitlab.com/api/v4/projects/26140156/packages/pypi

where `<token_name>` and `<token>` are obtained by creating a personal token
with `read_api`` scope. See `Creating personal access tokens
<https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html#creating-a-personal-access-token>`_
for details.

Molecular data
++++++++++++++

`molspecutils.molecule.COAlchemyMode` and
`molspecutils.molecule.CH3ClAlchemyMode` provide simple interfaces to data on
rovibrational transitions and energy levels of CO and CH3Cl :math:`\nu_3`
vibrational mode. The molecular data is not
bundled with the package and during the first run :mod:`molspecutils` will use
HAPI [1]_ to fetch transitions line lists from HITRAN [2]_ and cache them for
future use.

.. [1] R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, "HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data", J. Quant. Spectrosc. Radiat. Transfer *177*, 15-30 (2016).
.. [2] I.E. Gordon *et al.* "The HITRAN2016 molecular spectroscopic database", J. Quant. Spectrosc. Radiat. Transfer *203*, 3-69 (2017).
