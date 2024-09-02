if __name__ == '__main__':
    import argparse
    import prody
    import numpy as np

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--modesFn', type=str, required=True)
    parser.add_argument('--collecFn', type=str, required=True)
    parser.add_argument('--eigvalsFn', type=str, required=True)
    parser.add_argument('--gnmCheckFn', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)
    parser.add_argument('--collecThreshold', type=str, required=True)

    args = parser.parse_args()

    modes = prody.loadModel(args.modesFn)

    collectivity = prody.calcCollectivity(modes)
    if isinstance(collectivity, float):
        collectivityList = [collectivity]
    else:
        collectivityList = list(collectivity)

    np.savetxt(args.collecFn, collectivityList)

    eigvals = modes.getEigvals()
    np.savetxt(args.eigvalsFn, eigvals)
    np.savetxt(args.gnmCheckFn, [not modes.is3d()])

    idxSorted = [i[0] for i in sorted(enumerate(collectivityList), 
                                        key=lambda x: x[1], reverse=True)]
    numModes = modes.numModes()
    modeNum = modes.getIndices()

    score = []
    for _ in range(numModes):
        score.append(0)

    for i in range(numModes):
        score[idxSorted[i]] = idxSorted[i] + modeNum[i] + 2
    
    for i in range(numModes):
        score[i] = float(score[i]) / (2.0 * numModes)

    collectivityThreshold = float(args.collecThreshold)
    prody.writeScipionModes(args.folder, modes, scores=score, only_sqlite=True,
                            collectivityThreshold=collectivityThreshold)
