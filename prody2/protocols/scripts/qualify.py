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

    args = parser.parse_args()

    modes = prody.parseScipionModes(args.modesFn)

    collectivity = prody.calcCollectivity(modes)
    if isinstance(collectivity, float):
        collectivityList = [collectivity]
    else:
        collectivityList = list(collectivity)

    np.savetxt(args.collecFn, collectivityList)

    eigvals = modes.getEigvals()
    np.savetxt(args.eigvalsFn, eigvals)

    np.savetxt(args.gnmCheckFn, [not modes.is3d()])
