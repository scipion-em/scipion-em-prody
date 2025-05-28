if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join
    import numpy as np

    NMA_METRIC_OVERLAP = 0
    NMA_METRIC_COV_OVERLAP = 1
    NMA_METRIC_RWSIP = 2

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputPdbFns', type=str, required=True)
    parser.add_argument('--inputModesFns', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    parser.add_argument('--metric', type=str, required=True)

    parser.add_argument('--match', type=bool,
                        required=False, default=False)
    parser.add_argument('--diag', type=bool,
                        required=False, default=False)
    parser.add_argument('--norm', type=bool,
                        required=False, default=False)

    args = parser.parse_args()

    folder = args.folder
    pdb1, pdb2 = args.inputPdbFns.split()
    modesFn1, modesFn2 = args.inputModesFns.split()
    match = eval(args.match)

    modes1 = prody.parseScipionModes(modesFn1, pdb=pdb1)
    modes2 = prody.parseScipionModes(modesFn2, pdb=pdb2)

    nModes = np.max([modes1.numModes(), modes2.numModes()])
    nModesMin = np.min([modes1.numModes(), modes2.numModes()])

    if nModesMin != 1 and args.match:
        modeEns = prody.ModeEnsemble()
        modeEns.addModeSet(modes1)
        modeEns.addModeSet(modes2)
        modeEns.match()

        matchInds = prody.matchModes(modes1, modes2, index=True)

        prody.writeArray(join(folder, 'matchInds.txt'),
            np.array(matchInds, dtype=int)[1]+1,
            format='%3d')

        atoms = prody.parsePDB(pdb1)

        typeStr = str(type(modes2)).lower().split('.')[-1].split("'")[0]
        args.nmdFileName = join(folder, 'matched_modes.{0}.nmd'.format(typeStr))

        prody.writeNMD(args.nmdFileName, modeEns[1], atoms)
        prody.writeScipionModes(join(folder, ), modeEns[1], write_star=True)
    else:
        modeEns = [modes1, modes2]

    if args.metric == NMA_METRIC_OVERLAP:
        if args.norm:
            args.matrix = prody.calcOverlap(modeEns[0], modeEns[1],
                                            diag=args.diag)
        else:
            # Calculate direct dot product without vector normalisation found in calcOverlap
            args.matrix = modes1.getEigvecs().T @ modes2.getEigvecs()

        if args.matrix.ndim == 1:
            args.matrix.reshape(-1, 1)

    else:
        args.matrix = np.empty((nModes-6, 1))
        for i in range(6, nModes):
            if args.metric == NMA_METRIC_COV_OVERLAP:
                args.matrix[i-6, 0] = prody.calcEnsembleSpectralOverlaps(
                    modeEns[:, 6:i+1])[0, 1]
            else:
                args.matrix[i-6, 0] = prody.calcRWSIP(modeEns[0, 6:i+1],
                                                      modeEns[1, 6:i+1])

    prody.writeArray(join(folder, 'matrix.txt'), args.matrix,
                     format='%' + str(max([len(str(int(np.max(args.matrix)))),
                                           len(str(int(np.min(args.matrix))))]) + 4) + '.2f')
