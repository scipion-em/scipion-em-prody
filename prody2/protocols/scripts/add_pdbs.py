if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join, splitext, basename

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFns', type=str, required=True)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    folder = args.folder
    inputFns = args.inputFns.split()
    ags = prody.parsePDB(inputFns,
                         unite_chains=eval(args.uniteChains))

    outAg = ags[0]
    for ag in ags[1:]:
        outAg += ag

    pdbFileName = join(folder, 'joined_atoms.pdb')
    prody.writePDB(pdbFileName, outAg)

    fo = open(join(folder, 'pdb_data.txt'), 'w')
    fo.write('\t'.join([pdbFileName, str(outAg.numAtoms()),
                        str(outAg.numResidues()),
                        str(outAg.numChains())]) + '\n')
    fo.close()
