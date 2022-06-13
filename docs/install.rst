Installation
============

**Dependencies**

- attrs
- numpy
- scipy
- appdirs
- SQLAlchemy

**Install**

.. highlight:: sh

The package is available on the Python Package Index and can be most easily installed with `pip`::

  pip install molspecutils

Molecular data
++++++++++++++

:class:`molspecutils.molecule.COAlchemyMode` and
:class:`molspecutils.molecule.CH3ClAlchemyMode` provide simple interfaces to data on
rovibrational transitions and energy levels of CO and CH3Cl :math:`\nu_3`
vibrational mode. The molecular data is not
bundled with the package and during the first run :mod:`molspecutils` will use
HAPI [1]_ to fetch transitions line lists from HITRAN [2]_ and cache them for
future use.

.. [1] R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, "HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data", J. Quant. Spectrosc. Radiat. Transfer *177*, 15-30 (2016).
.. [2] I.E. Gordon *et al.* "The HITRAN2016 molecular spectroscopic database", J. Quant. Spectrosc. Radiat. Transfer *203*, 3-69 (2017).
