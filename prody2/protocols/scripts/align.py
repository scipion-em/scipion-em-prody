
# chain matching methods
BEST_MATCH = 0
SAME_CHID = 1
SAME_POS = 2
CUSTOM = 3

# residue mapping methods
NOTHING = 0 # stop trivial mapping if trivial mapping fails
PWALIGN = 1 # biopython pwalign local pairwise sequence alignment after trivial mapping
CEALIGN = 2 # combinatorial extension (CE) as in PyMOL
DEFAULT = 3 # try pwalign then CE

NOT_DUMMY_SELSTR = "not dummy"

if __name__ == '__main__':
    import argparse
    import numpy as np
    from os.path import join
    import prody

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--mobFn', type=str, required=True)
    parser.add_argument('--tarFn', type=str, required=False)
    parser.add_argument('--transformationFn', type=str, required=False)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    parser.add_argument('--matchFunc', type=str, required=False)
    parser.add_argument('--chmapFn', type=str, required=False)
    parser.add_argument('--mapping', type=str, required=False)
    parser.add_argument('--seqid', type=str, required=False)
    parser.add_argument('--overlap', type=str, required=False)
    parser.add_argument('--rmsdReject', type=str, required=False)
    parser.add_argument('--keepMismatching', type=str, required=False)

    args = parser.parse_args()

    folder = args.folder
    has_trans = args.transformationFn is not None

    mob = prody.parsePDB(args.mobFn, alt='all',
                         unite_chains=bool(args.uniteChains))

    if args.tarFn is not None:
        tar = prody.parsePDB(args.tarFn, alt='all',
                             unite_chains=bool(args.uniteChains))

        matchFuncId = int(args.matchFunc)
        if matchFuncId == BEST_MATCH:
            matchFunc = prody.bestMatch
        elif matchFuncId == SAME_CHID:
            matchFunc = prody.sameChid
        elif matchFuncId == SAME_POS:
            matchFunc = prody.sameChainPos
        else:
            fi = open(args.chmapFn, 'r')
            chmap = fi.readlines()
            fi.close()
            matchFunc = lambda chain1, chain2: prody.userDefined(chain1, chain2, chmap)

        mappingId = args.mapping
        if mappingId == DEFAULT:
            mapping = 'auto'
        elif mappingId == PWALIGN:
            mapping = 'pwalign'
        elif mappingId == CEALIGN:
            mapping = 'ce'
        else:
            mapping = False

        seqid = float(args.seqid)
        overlap = float(args.overlap)
        rmsdReject = float(args.rmsdReject)

        mobAmapList = prody.alignChains(mob.protein, tar.protein,
                                        seqid=seqid,
                                        overlap=overlap,
                                        match_func=matchFunc,
                                        mapping=mapping,
                                        rmsd_reject=rmsdReject)
        if len(mobAmapList):
            mobAmap = mobAmapList[0]
            mobSel = mobAmap.select(NOT_DUMMY_SELSTR).copy()
            mobSel.setTitle(mob.getTitle())

            tarAmapList = prody.alignChains(tar.protein, mobSel,
                                            seqid=seqid,
                                            overlap=overlap,
                                            match_func=matchFunc,
                                            mapping=mapping,
                                            rmsd_reject=rmsdReject)
            if len(tarAmapList):
                tarAmap = tarAmapList[0]
                tarSel = tarAmap.select(NOT_DUMMY_SELSTR).copy()
                tarSel.setTitle(tar.getTitle())

                if mobSel.numAtoms != tarSel.numAtoms():
                    mobAmapList = prody.alignChains(mobSel, tarSel,
                                                    seqid=seqid,
                                                    overlap=overlap,
                                                    match_func=matchFunc,
                                                    mapping=mapping,
                                                    rmsd_reject=rmsdReject)
                    if len(mobAmapList):
                        mobAmap = mobAmapList[0]
                        mobSel = mobAmap.select(NOT_DUMMY_SELSTR).copy()
                        mobSel.setTitle(mob.getTitle())

                    tarAmapList = prody.alignChains(tarSel, mobSel,
                                                    seqid=seqid,
                                                    overlap=overlap,
                                                    match_func=matchFunc,
                                                    mapping=mapping,
                                                    rmsd_reject=rmsdReject)
                    if len(tarAmapList):
                        tarAmap = tarAmapList[0]
                        tarSel = tarAmap.select(NOT_DUMMY_SELSTR).copy()
                        tarSel.setTitle(tar.getTitle())

                if has_trans is False:
                    T = prody.calcTransformation(mobSel, tarSel)
                else:
                    transformationMatrix = np.loadtxt(args.transformationFn)
                    T = prody.Transformation(transformationMatrix)

                keepMismatching = bool(args.keepMismatching)
                if keepMismatching:
                    alg = prody.applyTransformation(T, mob)
                    tarSel = tar
                else:
                    alg = prody.applyTransformation(T, mobSel)

                prody.writePDB(join(folder, 'mobile.pdb'), alg)
                prody.writePDB(join(folder, 'target.pdb'), tarSel)

                matrixFileName = join(folder, 'transformation.txt')
                prody.writeArray(matrixFileName, T.getMatrix())
    else:
        if has_trans:
            transformationMatrix = np.loadtxt(args.transformationFn)
            T = prody.Transformation(transformationMatrix)
            alg = prody.applyTransformation(T, mob)

            prody.writePDB(join(folder, 'mobile.pdb'), alg)
