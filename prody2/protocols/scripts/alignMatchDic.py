if __name__ == '__main__':
    import argparse
    from collections import OrderedDict
    from os.path import join
    import prody

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--mobFn', type=str, required=True)
    parser.add_argument('--tarFn', type=str, required=True)
    parser.add_argument('--uniteChains', type=str, required=True)
    parser.add_argument('--chainOrders', type=str, required=True)
    parser.add_argument('--customOrder', type=str, required=True)
    parser.add_argument('--index', type=str, required=True)
    parser.add_argument('--folder', type=str, required=True)

    args = parser.parse_args()

    def createMatchDic(mobFn, tarFn, uniteChains, chainOrders, customOrder, index):
        index = int(index)
        mob = prody.parsePDB(mobFn, alt='all',
                                  unite_chains=uniteChains)
        tar = prody.parsePDB(tarFn, alt='all',
                                  unite_chains=uniteChains)
        
        try:
            matchDic = eval(chainOrders)
            _ = matchDic.keys()
        except (AttributeError, TypeError):
            matchDic = OrderedDict()
            matchDic[mob.getTitle()] = getInitialChainOrder(mob)
            matchDic[tar.getTitle()] = getInitialChainOrder(tar)
            
        if index == 0:
            label = mob.getTitle()
            if customOrder == '':
                matchDic[label] = getInitialChainOrder(mob)
            else:
                matchDic[label] = customOrder
        else:
            label = tar.getTitle()
            if customOrder == '':
                matchDic[label] = getInitialChainOrder(tar)
            else:
                matchDic[label] = customOrder

        return matchDic
    
    def getInitialChainOrder(struct):
        return ''.join([ch.getChid() for ch in struct.iterChains()])

    matchDic = createMatchDic(args.mobFn, args.tarFn, args.uniteChains, args.chainOrders, 
                              args.customOrder, args.index)
    fo = open(join(args.folder, 'matchDic.txt'))
    fo.write(matchDic)
    fo.close()
