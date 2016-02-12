__author__ = 'ManorLab'

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(name='empanada',
      version='0.0.2',
      classifiers=['License :: Free for non-commercial use'],
      description='EMPANADA: a tool for evidence-based assignment of genes to pathways in metagenomic data',
      long_description=(read('README.rst') + '\n\n' + read('HISTORY.rst') + '\n\n' + read('AUTHORS.rst') + '\n\n' + read('LICENSE') + '\n\n'),
      author='Ohad Manor',
      author_email='omanor@gmail.com',
      url='http://elbo.gs.washington.edu/software_empanada.html',
      packages=['empanada'],
      package_data={'empanada': ['data/*.tab', 'examples/*.tab']},
      install_requires=['NumPy >= 1.6.1', 'pandas >= 0.14'],
      scripts=['scripts/run_empanada.py', 'tests/test_empanada.py'],
      )



