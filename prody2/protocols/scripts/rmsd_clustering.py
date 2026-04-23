
if __name__ == '__main__':
    import argparse
    import prody
    import numpy as np
    import os

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputEns', type=str, required=True)
    parser.add_argument('--rmsdThreshold', type=int, required=True)
    parser.add_argument('--outputDir', type=str, required=True)
    parser.add_argument('--treeMethod', type=str, required=True)

    args = parser.parse_args()

    ens = prody.loadEnsemble(args.inputEns)

    matrix = ens.getRMSDs(pairwise=True)
    labels = ens.getLabels()
    tree = prody.calcTree(labels, matrix, method=args.treeMethod)
    subgroups = prody.findSubgroups(tree, args.rmsdThreshold)

    counts = [len(sg) for sg in subgroups]
    classLabels = np.zeros(ens.numCoordsets(), dtype=int)
    repIdx = np.zeros(len(subgroups), dtype=int)
    sgIdx = []
    for i, sg in enumerate(subgroups):
        sgIdx.append([labels.index(label) for label in sg])
        submatrix = matrix[sgIdx[i], :][:, sgIdx[i]]
        repIdx[i] = sgIdx[i][np.argmin(np.mean(submatrix, axis=0))]
        classLabels[sgIdx[i]] = i

    reordMatrix, reordIndices = prody.reorderMatrix(labels, matrix, tree)

    np.savetxt(os.path.join(args.outputDir, "cluster_labels.txt"), 
               labels, fmt="%d")           # save cluster labels for each frame
    np.savetxt(os.path.join(args.outputDir, "cluster_reps.txt"), 
               repIdx, fmt="%d")           # save the representatives as frame numbers
    np.savetxt(os.path.join(args.outputDir, "cluster_counts.txt"), 
               counts, fmt="%d")           # save the number of members in each cluster

    np.savetxt(os.path.join(args.outputDir, "reordering_indices.txt"), 
               reordIndices, fmt="%d")     # save the reordering indices from reorderMatrix
    np.savetxt(os.path.join(args.outputDir, "reordered_rmsd_matrix.txt"), 
               reordMatrix, fmt="%d")      # save the reordered RMSD matrix
