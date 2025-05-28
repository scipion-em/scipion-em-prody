if __name__ == '__main__':
    import argparse
    import prody
    from os.path import join

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFns', type=str, required=True)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)
    parser.add_argument('--rmsd', type=float, required=True)

    parser.add_argument('--nSteps', type=int, required=True)
    parser.add_argument('--pos', type=bool, required=True)
    parser.add_argument('--neg', type=bool, required=True)


    args = parser.parse_args()

    folder = args.folder
    mobFn, tarFn = args.inputFns.split()


    mob = prody.parsePDB(mobFn, alt='all')
    tar = prody.parsePDB(tarFn, alt='all')

    rmsd = args.rmsd
    if rmsd == 0:
        rmsd = prody.calcRMSD(mob, tar)      

    defvec = prody.calcDeformVector(mob, tar)

    outModes = prody.NMA('defvec')
    outModes.setEigens(defvec.getArray().reshape(-1, 1))
    prody.writeScipionModes(folder, outModes, write_star=True)
    prody.writeNMD(join(folder, 'modes.nmd'), outModes, mob)

    animationsDir = join(folder, 'extra', 'animations')
    fnAnimation = join(animationsDir, "animated_mode_001")

    outAtoms = prody.traverseMode(defvec, mob, rmsd=rmsd,
                                  n_steps=args.nSteps,
                                  pos=args.pos,
                                  neg=args.neg)
    prody.writePDB(fnAnimation+".pdb", outAtoms)
