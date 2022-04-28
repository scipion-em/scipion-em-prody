"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from prody2 import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-prody',  # Required
    version=__version__,  # Required
    description='Scipion plugin for ProDy, the Bahar lab Python Package for Protein Dynamics Analysis',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-em/scipion-em-prody',  # Optional
    author='James Krieger',  # Optional
    author_email='jmkrieger@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing scipion-3.0 prody-2.0',  # Optional
    packages=find_packages(),
    entry_points={'pyworkflow.plugin': 'prody2 = prody2'},
    package_data={  # Optional
       'prody2': ['icon.png', 'protocols.conf'],
    },
    install_requires=['threadpoolctl==3.1.0']
)
