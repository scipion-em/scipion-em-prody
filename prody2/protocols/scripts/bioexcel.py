
if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--accession', type=str, required=True)
    parser.add_argument('--selection', type=str, required=True)
    parser.add_argument('--frames', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    selection = args.selection
    if selection == 'None':
        selection = None

    pdbFileName = prody.fetchBioexcelPDB(args.accession, folder=args.folder,
                                         selection=selection)
    ag = prody.parsePDB(pdbFileName)

    prody.fetchBioexcelTrajectory(args.accession, folder=args.folder,
                                  frames=args.frames, selection=selection)

    fo = open(join(args.folder, 'pdb_data.txt'), 'w')
    fo.write('\t'.join([pdbFileName, str(ag.numAtoms()),
                        str(ag.numResidues()),
                        str(ag.numChains())]) + '\n')
    fo.close()
