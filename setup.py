# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

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
