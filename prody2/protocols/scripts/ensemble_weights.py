if __name__ == '__main__':
    import argparse
    import numpy as np
    import os
    import prody

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--filename', type=str, required=True)
    
    args = parser.parse_args()
    ens = prody.loadEnsemble(os.path.join(args.path, args.filename))

    np.savetxt(os.path.join(args.path, 'weights.txt'), ens.getSizes())
    np.savetxt(os.path.join(args.path, 'labels.txt'), ens.getLabels(), fmt='%s')
