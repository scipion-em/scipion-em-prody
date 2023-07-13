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

import os
import pwem
import pyworkflow.utils as pwutils
from pyworkflow import Config

from .constants import *


__version__ = "3.3.0"
_logo = "icon.png"
_references = ['Zhang2021']


class Plugin(pwem.Plugin):
    _supportedVersions = VERSIONS
    _url = "https://github.com/scipion-em/scipion-em-prody"

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch ProDy. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addProDyPackage(env, ver,
                                default=ver == PRODY_DEFAULT_VER_NUM)

    @classmethod
    def addProDyPackage(cls, env, version, default=False):
        PRODY_INSTALLED = 'prody_%s_installed' % version
        installCmd = []

        if version == DEVEL:
            # Use latest scipion branch of prody on my github
            installCmd.append('cd .. &&')
            clonePath = os.path.join(pwem.Config.EM_ROOT, "ProDy")
            if not os.path.exists(clonePath):
                installCmd.append('git clone -b scipion https://github.com/jamesmkrieger/ProDy.git ProDy &&')
            installCmd.append('cd ProDy &&')

        # Install downloaded code
        installCmd.append('pip install -U -e . && python setup.py build_ext --inplace --force &&')

        if version == DEVEL:
            installCmd.append('cd .. && cd prody-github &&')

        installCmd.append('python -c "import os; os.environ.setdefault(\'HOME\', \'{0}\')" &&'.format(Config.SCIPION_HOME + os.path.sep))

        # Flag installation finished
        installCmd.append('touch %s' % PRODY_INSTALLED)

        # Install PDBFixer and OpenMM for ClustENM
        OPEN_MM_INSTALLED = 'openmm_installed'
        installOpenMM = 'conda install -c conda-forge pdbfixer -y && touch %s' % OPEN_MM_INSTALLED

        prody_commands = [(" ".join(installCmd), PRODY_INSTALLED),
                          (installOpenMM, OPEN_MM_INSTALLED)]

        envHome = os.environ.get('HOME', "")
        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there, and home since prody needs it to configure
        installEnvVars = {'PATH': envPath, 'HOME': envHome} if envPath else {'HOME': envHome}

        if version == DEVEL:
            env.addPackage('prody', version=version,
                            tar='void.tgz',
                            commands=prody_commands,
                            neededProgs=cls.getDependencies(),
                            default=default,
                            vars=installEnvVars)
        else:
            env.addPackage('prody', version=version,
                           url='https://github.com/prody/ProDy/archive/refs/tags/v{0}.tar.gz'.format(version),
                           buildDir='ProDy-{0}'.format(version),
                           commands=prody_commands,
                           neededProgs=cls.getDependencies(),
                           default=default,
                           vars=installEnvVars)            

    @classmethod
    def getProgram(cls, program):
        """ Create ProDy command line. """
        fullProgram = '%s && prody %s' % (
            cls.getCondaActivationCmd(),
            program)

        return fullProgram
