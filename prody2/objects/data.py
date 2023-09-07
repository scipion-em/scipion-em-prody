import numpy as np
import os
import prody
from pwem.objects import (EMObject, EMSet, Pointer, Integer,
                          String, SetOfNormalModes)
import logging
logger = logging.getLogger(__name__)

try:
    from pwem.objects import TrajFrame, SetOfTrajFrames
except ImportError:

    class TrajFrame(EMObject):
        """Represents an trajectory frame object"""

        def __init__(self, location=None, **kwargs):
            """
            Params:
            :param location: Could be a valid location: (index, filename)
            or  filename
            """
            EMObject.__init__(self, **kwargs)
            # Frame location is composed by an index and a filename
            self._index = Integer(0)
            self._filename = String()
            if location:
                self.setLocation(location)

        def getIndex(self):
            return self._index.get()

        def setIndex(self, index):
            self._index.set(index)

        def getFileName(self):
            """ Use the _objValue attribute to store filename. """
            return self._filename.get()

        def setFileName(self, filename):
            """ Use the _objValue attribute to store filename. """
            self._filename.set(filename)

        def getLocation(self):
            """ This function return the frame index and filename.
            It will only differ from getFileName, when the frame
            is contained in a trajectory and the index makes sense.
            """
            return self.getIndex(), self.getFileName()

        def setLocation(self, *args):
            """ Set the frame location, see getLocation.
            Params:
                First argument can be:
                1. a tuple with (index, filename)
                2. a index, this implies a second argument with filename
            """
            first = args[0]
            t = type(first)
            if t == tuple:
                index, filename = first
            elif t == int:
                index, filename = first, args[1]
            else:
                raise TypeError('setLocation: unsupported type %s as input.' % t)

            self.setIndex(index)
            self.setFileName(filename)

        def getBaseName(self):
            return os.path.basename(self.getFileName())

        def copyInfo(self, other):
            """ Copy basic information """
            self.copyAttributes(other, '_samplingRate')

        def copyLocation(self, other):
            """ Copy location index and filename from other frame. """
            self.setIndex(other.getIndex())
            self.setFileName(other.getFileName())

        def getFiles(self):
            filePaths = set()
            filePaths.add(self.getFileName())
            return filePaths


    class SetOfTrajFrames(EMSet):
        """ Represents a set of TrajFrames """
        ITEM_TYPE = TrajFrame

        def __init__(self, **kwargs):
            EMSet.__init__(self, **kwargs)

            self._ref = Pointer()
            self.setRef(kwargs.get('ref', None))
            if self._ref.get() is not None and not isinstance(self._ref.get(), TrajFrame):
                self._ref = None
                logger.warning("Reference frame should be a TrajFrame. "
                            "Whatever else was provided will be ignored.")

        def hasRef(self):
            return self._ref.hasValue()

        def getRef(self):
            """ Returns the reference frame associated with
            this SetOfTrajFrames"""
            return self._ref.get()

        def setRef(self, ref):
            """ Set the reference frame associates with
            this set of trajectory frames.
            """
            if ref is None:
                self._ref = Pointer()
            else:
                if ref.isPointer():
                    self._ref.copy(ref)
                else:
                    self._ref.set(ref)

                if not self._ref.hasExtended():
                    logger.warning("FOR DEVELOPERS: Direct pointers to objects should be avoided. "
                                "They are problematic in complex streaming scenarios. "
                                "Pass a pointer to a protocol with extended "
                                "(e.g.: input param are this kind of pointers. Without get()!)")

        def getFiles(self):
            filePaths = set()
            uniqueFiles = self.aggregate(['count'], '_filename', ['_filename'])

            for row in uniqueFiles:
                filePaths.add(row['_filename'])
            return filePaths

        def appendFromFrames(self, framesSet):
            """ Iterate over the frames and append
            every frame that is enabled.
            """
            for frame in framesSet:
                if frame.isEnabled():
                    self.append(frame)


class ProDyNpzEnsemble(SetOfTrajFrames):
    """Handles ProDy ens.npz files
    as sets of ensemble member frames"""

    def loadEnsemble(self, orderBy='id', direction='ASC'):
        """Make a new ensemble with all the items and write a new ens.npz file"""
        filenames = list(self.getFiles())
        oldEnsembles = [prody.loadEnsemble(filename) for filename in filenames]
        for i, ens in enumerate(oldEnsembles):
            if not isinstance(ens, prody.PDBEnsemble):
                oldEnsembles[i] = prody.PDBEnsemble(ens)

        newEnsemble = prody.PDBEnsemble()
        newSizes = []
        for i, item in enumerate(self.iterItems(orderBy=orderBy,
                                                direction=direction)):
            fnameIndex = filenames.index(item.getFileName())
            confIndex = item.getIndex() - 1 # back to python

            ensemble = oldEnsembles[fnameIndex]
            coords = ensemble.getCoordsets(selected=False)[confIndex]
            label = ensemble.getLabels()[confIndex]

            weights = ensemble.getWeights(selected=False)
            if weights is not None:
                weights = weights[confIndex]
            else:
                weights = 1.

            if ensemble.getData('size') is not None:
                newSizes.append(ensemble.getData('size')[confIndex])
            else:
                newSizes.append(np.zeros(ensemble.numCoordsets())[confIndex])

            if i == 0:
                newEnsemble.setCoords(ensemble.getCoords(selected=False))
                newEnsemble.setAtoms(ensemble.getAtoms())

            newEnsemble.addCoordset(coords, weights, label)

        newEnsemble.setData('size', newSizes)

        return newEnsemble

class SetOfGnmModes(SetOfNormalModes):
    pass
