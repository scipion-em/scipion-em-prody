if __name__ == '__main__':
    import argparse
    import prody
    import os

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    folder = args.folder
    prody.pathPDBFolder(folder)

    inputFn = prody.fetchPDB(args.pdb, compressed=False)
    if inputFn == None:
        inputFn = prody.fetchPDB(args.pdb, format="cif",
                                 compressed=False)

    prody.pathPDBFolder("") # release

    fo = open(os.path.join(folder, 'inputFn.txt'), 'w')
    fo.write(inputFn)
    fo.close()
