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


__version__ = 'devel'
_logo = "icon.png"
_references = ['Zhang2021']


class Plugin(pwem.Plugin):
    _supportedVersions = VERSIONS
    _url = "https://github.com/scipion-em/scipion-em-prody"

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(PRODY_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)

    @classmethod
    def getProDyEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(PRODY_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

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
        ENV_NAME = 'scipion3'
        # try to get CONDA activation command
        installCmd = [cls.getCondaActivationCmd()]

        # Activate the new environment
        installCmd.append('conda activate %s;' % ENV_NAME)

        if version == DEVEL:
            # Use latest scipion branch of prody on my github
            installCmd.append('git clone -b scipion https://github.com/jamesmkrieger/ProDy.git; cd ProDy;')

        # Install downloaded code
        installCmd.append('pip install -U -e .; python setup.py build_ext --inplace --force; ')

        if version == DEVEL:
            installCmd.append('cd ..; ')

        # configure ProDy to automatically handle secondary structure information
        installCmd.append('python -c "import os; os.environ.setdefault(\'HOME\', \'{0}\')"'.format(Config.SCIPION_HOME + os.path.sep))
        installCmd.append('python -c "import prody; prody.confProDy(auto_secondary=True)"')

        # Flag installation finished
        installCmd.append('&& touch %s' % PRODY_INSTALLED)

        prody_commands = [(" ".join(installCmd), PRODY_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None

        if version == DEVEL:
            env.addPackage('ProDy', version=version, dep=['threadpoolctl'],
                            buildDir='prody2', tar='void.tgz',
                            commands=prody_commands,
                            neededProgs=cls.getDependencies(),
                            default=default,
                            vars=installEnvVars)
        else:
            env.addPackage('ProDy', version=version, dep=['threadpoolctl'],
                           url='https://github.com/prody/prody/archive/v%s.tar.gz' % version,
                           commands=prody_commands,
                           neededProgs=cls.getDependencies(),
                           default=default,
                           vars=installEnvVars)            

    @classmethod
    def getProgram(cls, program):
        """ Create ProDy command line. """
        fullProgram = '%s %s && prody %s' % (
            cls.getCondaActivationCmd(), cls.getProDyEnvActivation(),
            program)

        return fullProgram

    @classmethod
    def getActiveVersion(cls, *args):
        """ Return the env name that is currently active. """
        envVar = cls.getVar(PRODY_ENV_ACTIVATION)
        return envVar.split()[-1]

    @classmethod
    def IS_V210(cls):
        return cls.getActiveVersion().startswith(getProDyEnvName('2.1.0'))
