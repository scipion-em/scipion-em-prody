if __name__ == '__main__':
    import argparse
    import prody
    import numpy as np

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--modesFn', type=str, required=True)
    parser.add_argument('--atomsFn', type=str, required=True)
    parser.add_argument('--coeffsFn', type=str, required=True)
    parser.add_argument('--numCoeffs', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)
    parser.add_argument('--nmdFileName', type=str, required=True)
    parser.add_argument('--npzFileName', type=str, required=True)

    args = parser.parse_args()

    modes = prody.parseScipionModes(args.modesFn)
    atoms = prody.parsePDB(args.atomsFn)
    coeffs = np.loadtxt(args.coeffsFn)

    vector = modes[0] * coeffs[0]

    numCoeffs = int(args.numCoeffs)
    for i in range(1, numCoeffs):
        vector += modes[i] * coeffs[i]

    outModes = prody.NMA()
    outModes.setEigens(vector.getArray().reshape(-1,1))
    prody.writeScipionModes(args.folder, outModes, write_star=True)

    prody.writeNMD(args.nmdFileName, outModes, atoms)
    prody.saveModel(outModes, args.npzFileName)
