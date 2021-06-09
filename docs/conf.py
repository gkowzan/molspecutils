# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'molspecutils'
copyright = '2020, Grzegorz Kowzan'
author = 'Grzegorz Kowzan'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'sphinx.ext.mathjax',
    # 'sphinx.ext.autosummary',
    # 'numpydoc'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinxdoc'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         'searchbox.html',
#     ]
# }
# html_show_sourcelink = True

# -- Intersphinx
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'matplotlib': ('https://matplotlib.org/', None),
    'xarray': ('http://xarray.pydata.org/en/stable/', None),
    'rotsim2d': ('/home/grz/python/packages/rotsim2d/docs/_build/html', None),
    'anytree': ('https://anytree.readthedocs.io/en/latest/', None)
}

# -- pygments
pygments_style = 'sphinx'

# -- numpydoc
# numpydoc_class_members_toctree = False
# numpydoc_show_class_members = False
# numpydoc_show_inherited_class_members = False

# -- autodoc
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}
autoclass_content = 'both'
autodoc_member_order = 'bysource'