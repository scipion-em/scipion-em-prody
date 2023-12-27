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
_references = ['ProDy2']


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

        ENV_NAME = getProDyEnvName(version)

        installCmd = [
            cls.getCondaActivationCmd(),
            f'conda create -y -n {ENV_NAME} python=3.9 &&',
            f'conda activate {ENV_NAME} &&']

        # Install TEMPy for ClustENM fitting and threadpoolctl for control of thread pools for apps generally
        TEMPY_INSTALLED = 'tempy_installed'
        installTEMPy = installCmd.copy()
        installTEMPy.append('pip install biotempy==2.0.0 threadpoolctl && touch %s' % TEMPY_INSTALLED)
        installCmd.pop(1) # remove conda create to only do it the first time

        # Install PDBFixer and OpenMM for ClustENM
        OPENMM_INSTALLED = 'openmm_installed'
        installOpenMM = installCmd.copy()
        installOpenMM.append('conda install -c conda-forge openmm==7.6 pdbfixer -y && touch %s' % OPENMM_INSTALLED)

        prodyCommands = [(" ".join(installTEMPy), TEMPY_INSTALLED),
                         (" ".join(installOpenMM), OPENMM_INSTALLED)]

        PRODY_INSTALLED_OWN = 'prody_%s_installed_own_env' % version
        PRODY_INSTALLED_SCIPION = 'prody_%s_installed_scipion_env' % version
        for i, PRODY_INSTALLED in enumerate([PRODY_INSTALLED_OWN, PRODY_INSTALLED_SCIPION]):
            if i == 1:
                installCmd = []

            if version == DEVEL:
                # Use latest prody on github
                installCmd.append('cd .. &&')
                clonePath = os.path.join(pwem.Config.EM_ROOT, "ProDy")
                if not os.path.exists(clonePath):
                    installCmd.append('git clone -b scipion https://github.com/jamesmkrieger/ProDy.git ProDy &&')
                installCmd.append('cd ProDy &&')
                installCmd.append('git pull &&')

            # Install downloaded code
            installCmd.append('pip install -U -e . && python setup.py build_ext --inplace --force &&')

            if version == DEVEL:
                installCmd.append('cd .. && cd prody-github &&')

            installCmd.append('python -c "import os; os.environ.setdefault(\'HOME\', \'{0}\')" &&'.format(Config.SCIPION_HOME + os.path.sep))

            # Flag installation finished
            installCmd.append('touch %s' % PRODY_INSTALLED)

            prodyCommands.append((" ".join(installCmd.copy()), PRODY_INSTALLED))

        envHome = os.environ.get('HOME', "")
        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there, and home since prody needs it to configure
        installEnvVars = {'PATH': envPath, 'HOME': envHome} if envPath else {'HOME': envHome}

        if version == DEVEL:
            env.addPackage('prody', version=version,
                            tar='void.tgz',
                            commands=prodyCommands,
                            neededProgs=cls.getDependencies(),
                            default=default,
                            vars=installEnvVars)
        else:
            env.addPackage('prody', version=version,
                           url='https://github.com/prody/ProDy/archive/refs/tags/v{0}.tar.gz'.format(version),
                           buildDir='ProDy-{0}'.format(version),
                           commands=prodyCommands,
                           neededProgs=cls.getDependencies(),
                           default=default,
                           vars=installEnvVars)            

    @classmethod
    def getProgram(cls, program, script=False):
        """ Create ProDy command line. """
        if script:
            fullProgram = '%s %s && python %s' % (
                cls.getCondaActivationCmd(), cls.getEnvActivation(),
                PRODY_SCRIPTS+'/'+program)
        else:
            fullProgram = '%s %s && prody %s' % (
                cls.getCondaActivationCmd(), cls.getEnvActivation(),
                program)

        return fullProgram

    @classmethod
    def getEnvActivation(cls):
        return "conda activate %s" % getProDyEnvName(DEVEL)
