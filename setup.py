from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    name='spectroscopy',
    version='1.0.0',
    description="Cavity-enhanced linear spectroscopy",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://gitlab.com/allisonlab/mdcs/spectroscopy',
    author='Grzegorz Kowzan',
    author_email='grzegorz@kowzan.eu',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Physics'],
    python_requires='>=3.5',
    install_requires=['numpy', 'scipy', 'xarray', 'matplotlib', 'pyfftw',
                      'shed @ git+ssh://git@gitlab.com/allisonlab/mdcs/shed.git@master'],
    packages=find_packages()
)
