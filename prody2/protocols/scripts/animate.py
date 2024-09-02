if __name__ == '__main__':
    import argparse
    import prody
    import numpy as np
    from os.path import join

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--modesFn', type=str, required=True)
    parser.add_argument('--atomsFn', type=str, required=True)
    parser.add_argument('--animationsDir', type=str, required=True)

    parser.add_argument('--nzero', type=str, required=True)
    parser.add_argument('--rmsd', type=str, required=True)
    parser.add_argument('--nSteps', type=str, required=True)
    parser.add_argument('--pos', type=str, required=True)
    parser.add_argument('--neg', type=str, required=True)

    args = parser.parse_args()

    modes = prody.loadModel(args.modesFn)
    atoms = prody.parsePDB(args.atomsFn)

    nzero = int(args.nzero)
    rmsd = float(args.rmsd)
    nSteps = int(args.nSteps)
    pos = bool(args.pos)
    neg = bool(args.neg)

    for i, mode in enumerate(modes[nzero:]):
        modenum = i+nzero+1
        fnAnimation = join(args.animationsDir, "animated_mode_%03d"
                        % modenum)
        prody.writePDB(fnAnimation+".pdb", 
                        prody.traverseMode(mode, atoms, rmsd=rmsd, n_steps=nSteps,
                                            pos=pos, neg=neg)
                        )

        fhCmd=open(fnAnimation+".vmd",'w')
        fhCmd.write("mol new %s.pdb\n" % fnAnimation)
        fhCmd.write("animate style Rock\n")
        fhCmd.write("display projection Orthographic\n")
        fhCmd.write("mol modcolor 0 0 Index\n")

        numAtomsP = numAtomsCA = 0
        if atoms.select('name P') is not None:
            numAtomsP = atoms.select('name P').numAtoms()
        if atoms.ca is not None:
            numAtomsCA = atoms.ca.numAtoms()

        numAtomsRep = numAtomsCA + numAtomsP
        if numAtomsRep == atoms.numAtoms():
            fhCmd.write("mol modstyle 0 0 Beads 2.000000 8.000000\n")
            # fhCmd.write("mol modstyle 0 0 Beads 1.800000 6.000000 "
            #         "2.600000 0\n")
        else:
            fhCmd.write("mol modstyle 0 0 NewRibbons 1.800000 6.000000 "
                    "2.600000 0\n")
        fhCmd.write("animate speed 0.5\n")
        fhCmd.write("animate forward\n")
        fhCmd.close()
