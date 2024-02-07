
if __name__ == '__main__':
    import argparse
    import prody
    import numpy as np
    import os

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputEns', type=str, required=True)
    parser.add_argument('--nClusters', type=int, required=True)
    parser.add_argument('--outputDir', type=str, required=True)

    args = parser.parse_args()

    ens = prody.loadEnsemble(args.inputEns)
    medoids, labels, counts = prody.calcKmedoidClusters(ens.getCoordsets(), 
                                                        args.nClusters)

    np.savetxt(os.path.join(args.outputDir, "cluster_labels.txt"), 
               labels, fmt="%d")      # save cluster labels for each frame
    np.savetxt(os.path.join(args.outputDir, "cluster_medoids.txt"), 
               medoids, fmt="%d")    # save the medoids as frame numbers
    np.savetxt(os.path.join(args.outputDir, "cluster_counts.txt"), 
               counts, fmt="%d")      # save the number of members in each cluster
