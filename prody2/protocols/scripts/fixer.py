
if __name__ == '__main__':
    import argparse
    import prody

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFn', type=str, required=True)
    parser.add_argument('--pH', type=float, required=True)
    parser.add_argument('--outputFn', type=str, required=True)

    args = parser.parse_args()

    prody.addMissingAtoms(args.inputFn, pH=args.pH, outfile=args.outputFn, 
                          method='pdbfixer', model_residues=True)
