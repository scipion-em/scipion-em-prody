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
from collections import OrderedDict
import os
import pwem
import pyworkflow.utils as pwutils
from pyworkflow import Config

from .constants import *


__version__ = "3.4.0"
_logo = "icon.png"
_references = ['ProDy2']


class Plugin(pwem.Plugin):
    _supportedVersions = VERSIONS
    _url = "https://github.com/scipion-em/scipion-em-prody"

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(PRODY_ENV_ACT, "conda activate prody-{0}".format(PRODY_DEFAULT_VER_NUM))

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
                                default=(ver==PRODY_DEFAULT_VER_NUM))

    @classmethod
    def addProDyPackage(cls, env, version, default=False):

        ENV_NAME = getProDyEnvName(version)

        installCmd = [
            cls.getCondaActivationCmd(),
            f'conda create -y -n {ENV_NAME} python=3.9 &&',
            f'conda activate {ENV_NAME} &&']

        # Install TEMPy for ClustENM fitting, scikit-learn-extra for Kmedoids
        # and threadpoolctl for control of thread pools for apps generally
        TEMPY_INSTALLED = 'tempy_installed'
        installTEMPy = installCmd.copy()
        installTEMPy.append('pip install biotempy==2.0.0 scikit-learn-extra '
                            'threadpoolctl requests mdtraj pyparsing==3.1.1 && touch %s' % TEMPY_INSTALLED)
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
            if i == 0:
                
                installCmd.append('git clone https://github.com/jamesmkrieger/ProDy.git ProDy &&')
                installCmd.append('cd ProDy &&')
                installCmd.append('git fetch &&')

                installCmd.append('git remote add upstream https://github.com/prody/ProDy.git &&')
                installCmd.append('git fetch upstream &&')
                installCmd.append('git checkout upstream/main &&')
                
                installCmd.append('git checkout scipion &&')
                installCmd.append('git pull &&')

                installCmd.append('pip install -Ue . && python setup.py build_ext --inplace --force &&')
            else:
                installCmd = []
                installCmd.append('cd ProDy &&')
                installCmd.append('pip install -Ue . &&')
            
            installCmd.append('cd .. && touch %s' % PRODY_INSTALLED)

            prodyCommands.append((" ".join(installCmd.copy()), PRODY_INSTALLED))

        envHome = os.environ.get('HOME', "")
        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there, and home since prody needs it to configure
        installEnvVars = {'PATH': envPath, 'HOME': envHome} if envPath else {'HOME': envHome}

        env.addPackage('prody', version=version,
                        tar='void.tgz',
                        buildDir='ProDy',
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
        return cls.getVar(PRODY_ENV_ACT)

def parseMatchDict(cls):
    if cls.chainOrders.get() != "":
        cls.matchDic = eval(cls.chainOrders.get())
    else:
        cls.matchDic = OrderedDict()

    if not isinstance(cls.matchDic, OrderedDict):
        cls.matchDic = OrderedDict()

    cls.labels = list(cls.matchDic.keys())
