
if __name__ == '__main__':
    import argparse
    import prody

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--accession', type=str, required=True)
    parser.add_argument('--selection', type=str, required=True)
    parser.add_argument('--frames', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    prody.fetchBioexcelPDB(args.accession, folder=args.folder,
                           selection=args.selection)
    prody.fetchBioexcelTrajectory(args.accession, folder=args.folder,
                                  frames=args.frames, selection=args.selection)
