import prody
from pwem.objects import TrajFrame, SetOfTrajFrames

class ProDyNpzEnsemble(SetOfTrajFrames):
    """Class for handling ProDy ens.npz files
    as sets of ensemble members"""

    def writeNPZ(self, fnNPZ, orderBy='id', direction='ASC'):
        """Make a new ensemble with all the items and write a new ens.npz file"""
        filenames = list(self.getFiles())
        old_ensembles = [prody.loadEnsemble(filename) for filename in filenames]

        new_ensemble = prody.Ensemble()
        if self.hasRef():
            if isinstance(self._ref.get(), TrajFrame):
                ref = self._ref.get()
                ref_ens = prody.loadEnsemble(ref.getFileName())
                new_ensemble.setCoords(ref_ens.getCoordsets(ref.getIndex(), 
                                                            selected=False)[0])
                new_ensemble.setAtoms(ref_ens.getAtoms())

        for i, item in enumerate(self.iterItems(orderBy=orderBy,
                                                direction=direction)):
            fname_index = filenames.index(item.getFileName())
            conf_index = item.getIndex()

            ensemble = old_ensembles[fname_index]
            coords = ensemble.getCoordsets(conf_index, selected=False)[0]
            if not self.hasRef():
                new_ensemble.setCoords(coords)
                new_ensemble.setAtoms(ensemble.getAtoms())
                self.setRef(item)
            new_ensemble.addCoordset(coords)

        #return new_ensemble
        prody.saveEnsemble(new_ensemble, fnNPZ)
