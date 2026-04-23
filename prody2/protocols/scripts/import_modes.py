
if __name__ == '__main__':
    import argparse
    import os
    import prody

    NMD = 0
    MODES_NPZ = 1
    SCIPION = 2

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--importType', type=str, required=True)
    parser.add_argument('--filesPaths', type=str, required=True)
    parser.add_argument('--filesPattern2', type=str, required=True)
    parser.add_argument('--protocolPath', type=str, required=True)
    parser.add_argument('--pdbFileName', type=str, required=True)

    args = parser.parse_args()

    filesPaths = eval(args.filesPaths)
    pdbFilename = args.pdbFileName
    protocolPath = args.protocolPath

    if args.importType == SCIPION:
        filePath = filesPaths[0]
    else:
        folderPath = os.path.split(filesPaths[0])[0]
        filesPattern1 = os.path.split(filesPaths[0])[1]

    if pdbFilename != "None":
        outPdbfileName = os.path.join(protocolPath, os.path.basename(pdbFilename))
        if not pdbFilename.endswith('atoms.pdb'):
            outPdbfileName = outPdbfileName[:4] + '_atoms.pdb'
        prody.writePDB(outPdbfileName, prody.parsePDB(pdbFilename))

    if args.importType == NMD:
        if not filesPattern1.endswith('.nmd'):
            filesPattern1 += '.nmd'

        if filesPattern1.find('pca') != -1:
            prodyType = prody.PCA
        elif filesPattern1.find('lra') != -1:
            prodyType = prody.LRA
        elif filesPattern1.find('gnm') != -1:
            prodyType = prody.GNM
        else:
            prodyType = prody.NMA

        outModes, atoms = prody.parseNMD(os.path.join(folderPath, filesPattern1),
                                                    type=prodyType)

        if pdbFilename == "None":
            pdbFilename = prody.writePDB(os.path.join(protocolPath, 'extra/atoms'),
                                         atoms)

    elif args.importType == MODES_NPZ:
        if not filesPattern1.endswith('.npz'):
            filesPattern1 += '.npz'
        outModes = prody.loadModel(os.path.join(folderPath, filesPattern1))

    elif args.importType == SCIPION:
        outModes = prody.parseScipionModes(filePath, pdb=pdbFilename)

    else:
        pattern2 = args.filesPattern2.get()
        outModes = prody.parseGromacsModes(folderPath, eigval_fname=filesPattern1,
                                                eigvec_fname=pattern2, average_pdb=pdbFilename)

    prody.writeScipionModes(protocolPath, outModes, write_star=True)
    
    if args.importType != NMD:
        atoms = prody.parsePDB(pdbFilename)
        typeStr = str(type(outModes)).lower().split('.')[-1].split("'")[0]
        nmdFileName = os.path.join(protocolPath, 'modes.{0}.nmd'.format(typeStr))
        prody.writeNMD(nmdFileName, outModes, atoms)
