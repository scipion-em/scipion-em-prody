# **************************************************************************
# *
# * Authors:     James Krieger (jamesmkrieger@gmail.com)
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
import numpy
import Bio
from pwem.convert.atom_struct import cifToPdb
import pyworkflow.utils as pwutils

from .constants import *


__version__ = "3.4.0"
_logo = "icon.png"
_references = ['ProDy2']


file_path = os.path.abspath(__file__)
dir_path = os.path.split(os.path.split(file_path)[0])[0]

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
        ENV_YAML_PATH = os.path.join(dir_path, 'myenv.yaml')
        prodyCommands = []

        PRODY_V241_INSTALLED = 'prody_v2.4.1_installed'
        installCmd = [cls.getCondaActivationCmd()]
        installCmd.append('pip install -U ProDy==2.4.1 &&')
        installCmd.append('pip install numpy=={0} biopython=={1} &&'.format(numpy.__version__,
                                                                            Bio.__version__))
        installCmd.append('touch %s' % PRODY_V241_INSTALLED)
        prodyCommands.append((" ".join(installCmd.copy()), PRODY_V241_INSTALLED))

        PRODY_INSTALLED = 'prody_%s_installed' % version
        installProDyGithub = [
            cls.getCondaActivationCmd(),
            f'conda env create -f {ENV_YAML_PATH} -n {ENV_NAME} &&',
            f'conda activate {ENV_NAME} &&']
        installProDyGithub.append('git clone https://github.com/jamesmkrieger/ProDy.git ProDy &&')
        installProDyGithub.append('cd ProDy &&')
        installProDyGithub.append('git checkout scipion &&')
        installProDyGithub.append('git pull &&')
        installProDyGithub.append('pip install -Ue . && python setup.py build_ext --inplace --force &&')
        installProDyGithub.append('cd .. && touch %s' % PRODY_INSTALLED)
        prodyCommands.append((" ".join(installProDyGithub.copy()), PRODY_INSTALLED))

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

def copyConvertPDB(infilename, outfilename):
    from os import path, symlink
    import shutil

    extension = path.splitext(infilename)[1]
    if extension == ".pdb":
        shutil.copy(infilename, outfilename)
    elif extension == '.cif':
        cifToPdb(infilename, outfilename)
