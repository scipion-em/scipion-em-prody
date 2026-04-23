if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join, splitext, basename

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFn', type=str, required=True)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    folder = args.folder
    inputFn = args.inputFn

    ags = prody.parsePDB(inputFn, alt='all', compressed=False,
                         biomol=True, extend_biomol=True,
                         unite_chains=eval(args.uniteChains))
    if isinstance(ags, prody.AtomGroup):
        ags = [ags]

    fo = open(join(folder, 'pdb_data.txt'), 'w')
    for i, ag in enumerate(ags):
        filename = join(folder, splitext(basename(inputFn))[0] \
                        + '_atoms_{0}.pdb'.format(i))
        prody.writePDB(filename, ag)
        fo.write('\t'.join([filename, str(ag.numAtoms()),
                           str(ag.numResidues()),
                           str(ag.numChains())]) + '\n')
    fo.close()
