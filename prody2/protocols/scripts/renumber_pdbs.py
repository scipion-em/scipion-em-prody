if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join, splitext, basename

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFn', type=str, required=True)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    parser.add_argument('--selection', type=str, required=True)
    parser.add_argument('--offset', type=int, required=True)
    parser.add_argument('--chain', type=str, required=True)

    args = parser.parse_args()
    folder = args.folder
    inputFn = args.inputFn

    ag = prody.parsePDB(inputFn)

    sel = ag.select(args.selection)
    sel.setResnums(sel.getResnums() + args.offset)

    chain = args.chain
    if chain.strip() != '':
        sel.setChids(chain)

    pdbFileName = join(folder, 'renum_atoms.pdb')
    prody.writePDB(pdbFileName, ag)

    resnums = ag.getResnums()

    fo = open(join(folder, 'pdb_data.txt'), 'w')
    fo.write('\t'.join([pdbFileName, str(ag.numAtoms()),
                        str(ag.numResidues()), str(ag.numChains()),
                        str(resnums[0]), str(resnums[-1]),
                        str(max(resnums)), str(min(resnums))]) + '\n')
    fo.close()
