if __name__ == '__main__':
    import argparse
    import prody
    from prody.utilities import ZERO
    from os.path import join
    import numpy as np

    NMA_METRIC_OVERLAP = 0
    NMA_METRIC_COV_OVERLAP = 1
    NMA_METRIC_RWSIP = 2

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputPdb', type=str, required=True)
    parser.add_argument('--inputModes', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)
    parser.add_argument('--nModes', type=int, required=True)
    parser.add_argument('--outputPdb', type=str, required=True)

    args = parser.parse_args()

    modesPath = args.folder
    pdb1 = args.inputPdbFns
    modesFn1 = args.inputModesFns

    modes = prody.parseScipionModes(modesFn1, pdb=pdb1)
    atoms = prody.parsePDB(pdb1)

    numModes = args.nModes
    mode = modes[:numModes]

    domains = prody.calcGNMDomains(mode)
    prody.writePDB(args.outputPdb, atoms, beta=domains)
