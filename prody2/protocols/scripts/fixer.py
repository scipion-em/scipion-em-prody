
if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join, splitext, basename

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFn', type=str, required=True)
    parser.add_argument('--pH', type=float, required=True)
    parser.add_argument('--outputFn', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    filename = prody.addMissingAtoms(args.inputFn, pH=args.pH, outfile=args.outputFn,
                                     method='pdbfixer', model_residues=True)
    ag = prody.parsePDB(filename)

    fo = open(join(args.folder, 'pdb_data.txt'), 'w')
    fo.write('\t'.join([filename,
                        str(ag.numAtoms()),
                        str(ag.numResidues()),
                        str(ag.numChains())]) + '\n')
    fo.close()
