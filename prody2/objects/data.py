import numpy as np
import os
import prody
from prody2.constants import ENSEMBLE_WEIGHTS
from pwem.objects import (EMObject, EMSet, SetOfNormalModes, SetOfClasses3D,
                          Pointer, Integer, Float, String,
                          AtomStruct, SetOfAtomStructs)
from pwem.convert import AtomicStructHandler
from pyworkflow.utils import logger

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
        setattr(self, ENSEMBLE_WEIGHTS, Float(kwargs.get('weight')))

    def getIndex(self):
        return self._index.get()

    def setIndex(self, index):
        self._index.set(index)

    def getWeight(self):
        return getattr(self, ENSEMBLE_WEIGHTS)

    def setWeight(self, weight):
        setattr(self, ENSEMBLE_WEIGHTS, weight)

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
            2. an index, this implies a second argument with filename
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
        """ Set the reference frame associated with
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
                newSizes.append(np.ones(ensemble.numCoordsets())[confIndex])

            if i == 0:
                newEnsemble.setCoords(ensemble.getCoords(selected=False))
                newEnsemble.setAtoms(ensemble.getAtoms())

            newEnsemble.addCoordset(coords, weights, label)

        newEnsemble.setData('size', newSizes)

        return newEnsemble

def replaceCoordsets(oldNpzEns, coordsets, suffix=''):
    oldEnsemble = oldNpzEns.loadEnsemble()

    newEnsemble = prody.PDBEnsemble()
    newEnsemble.addCoordset(coordsets, 
                            weights=oldEnsemble.getWeights(), 
                            label=oldEnsemble.getLabels())
    newEnsemble.setAtoms(oldEnsemble.getAtoms())

    newFilename = oldNpzEns[1].getFileName().replace('.ens.npz', f'{suffix}.ens.npz')
    prody.saveEnsemble(newEnsemble, newFilename)

    newNpzEns = ProDyNpzEnsemble().create(os.path.split(newFilename)[0])

    frames = [frame.clone() for frame in oldNpzEns.iterItems()]
    for i, frame in enumerate(frames):
        frame.setLocation((i+1, newFilename))
        newNpzEns.append(frame)

    return newNpzEns


class SetOfGnmModes(SetOfNormalModes):
    pass

class SetOfLdaModes(SetOfNormalModes):
    pass

class SetOfLogisticModes(SetOfLdaModes):
    pass

class Atom(EMObject):
    """Atom information"""

    def __init__(self, location=None, **kwargs):
        EMObject.__init__(self, **kwargs)

        # Frame location is composed by an index and a filename
        self._index = Integer(0)
        self._filename = String()
        if location:
            # create with new data, else read existing data
            self.setLocation(location)

            handler = AtomicStructHandler()
            handler.read(self.getFileName())
            struct = handler.getStructure()
            self._bioAtom = list(struct.get_atoms())[self.getIndex()]

            self._serial = Integer(self._bioAtom.serial_number)
            self._name = String(self._bioAtom.name)

            description = self._bioAtom.get_full_id()
            self._chid = String(description[2])
            self._resnum = Integer(description[3][1])
            self._resname = self._bioAtom.get_parent().get_resname()

            coords = self._bioAtom.get_coord()
            self._x = Float(coords[0])
            self._y = Float(coords[1])
            self._z = Float(coords[2])

            self._occupancy = Float(self._bioAtom.occupancy)
            self._bfactor = Float(self._bioAtom.bfactor)
            self._element = String(self._bioAtom.element)

    def copyInfo(self, other):
        self.copy(other, copyId=False)

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
            2. an index, this implies a second argument with filename
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

    def getIndex(self):
        return self._index.get()

    def setIndex(self, index):
        self._index.set(index)

    def getFileName(self):
        return self._filename.get()

    def setFileName(self, index):
        self._filename.set(index)

    def getSerial(self):
        return self._serial.get()

    def setSerial(self, value):
        self._serial.set(value)

    def getName(self):
        return self._name.get()

    def setName(self, value):
        self._name.set(value)

    def getChid(self):
        return self._chid.get()

    def setChid(self, value):
        self._chid.set(value)

    def setResnum(self, value):
        self._resnum.set(value)

    def getResnum(self):
        return self._resnum.get()

    def setResname(self, value):
        self._resname.set(value)

    def getResname(self):
        return self._resname.get()

    def getX(self):
        return self._x.get()

    def setX(self, value):
        self._x.set(value)

    def setY(self, value):
        self._y.set(value)

    def getY(self):
        return self._y.get()

    def setZ(self, value):
        self._x.set(value)

    def getZ(self):
        return self._x.get()

    def getOccupancy(self):
        return self._occupancy.get()

    def setOccupancy(self, value):
        self._occupancy.set(value)

    def setBfactor(self, value):
        self._bfactor.set(value)

    def getBfactor(self):
        return self._bfactor.get()

    def setElement(self, value):
        self._element.set(value)

    def getElement(self):
        return self._element.get()

    def __str__(self):
        return "\n    serial=%s\n    name=%s\n    chid=%s\n     resname=%s\n    resnum=%s\n\n" % \
               (self._serial.get(),
                self._name.get(),
                self._chid.get(),
                self._resname.get(),
                self._resnum.get())

class SetOfAtoms(EMSet):
    """ Represents a set of Atoms """
    ITEM_TYPE = Atom

class ClassTraj(SetOfTrajFrames):
    """Class from clustering trajectories"""
    REP_TYPE = TrajFrame

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not
        _mapperPath or _size from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])

    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfClasses
        pass

class SetOfClassesTraj(SetOfClasses3D):
    """ SetOfClasses from trajectories"""
    ITEM_TYPE = ClassTraj
    REP_TYPE = AtomStruct
    REP_SET_TYPE = SetOfAtomStructs


def loadAndWriteEnsemble(cls):
    """Handle inputs to load ensemble into ProDy and write outputs"""

    if isinstance(cls.inputEnsemble, Pointer):
        inputEnsemble = [cls.inputEnsemble.get()]
    else:
        inputEnsemble = [ensemble.get() for ensemble in cls.inputEnsemble]

    for i, ensemble in enumerate(inputEnsemble):
        if isinstance(ensemble, SetOfAtomStructs):
            ags = prody.parsePDB([tarStructure.getFileName() for tarStructure in ensemble])
            ens = prody.buildPDBEnsemble(ags, match_func=prody.sameChainPos, seqid=0.,
                                        overlap=0., superpose=False, degeneracy=cls.degeneracy.get())
            # the ensemble gets built exactly as the input is setup and nothing gets rejected
        else:
            ens = inputEnsemble[i].loadEnsemble()

        if i == 0:
            cls.ens = ens
        else:
            cls.ens += ens

    cls.ens.select(cls.selstr.get())

    if os.path.exists(cls._getPath()):
        avgStruct = cls.ens.getAtoms()
        avgStruct.setCoords(cls.ens.getCoords())

        cls.pdbFileName = cls._getPath('atoms.pdb')
        prody.writePDB(cls.pdbFileName, avgStruct)
        cls.averageStructure = AtomStruct()
        cls.averageStructure.setFileName(cls.pdbFileName)

        cls.dcdFileName = cls._getPath('ensemble.dcd')
        prody.writeDCD(cls.dcdFileName, cls.ens)

        cls.npzFileName = cls._getPath('ensemble.ens.npz')
        prody.saveEnsemble(cls.ens, cls.npzFileName)
        cls.npz = ProDyNpzEnsemble().create(cls._getPath())
        for j in range(cls.ens.numConfs()):
            frame = TrajFrame((j+1, cls.npzFileName), objLabel=cls.ens.getLabels()[j])
            cls.npz.append(frame)