import numpy as np

class Block:
    '''
    input:
    f = (1,) array or (N,) array or list
    
    use:
    aver, std = smda.Block().block(f, 5)
    '''
    def __init__(self):
        pass

    def block(self, f, nblocks):
        f = np.array(f)
        f = f[(len(f)%nblocks):]
        sf = np.array(np.split(f, nblocks))
        blocks = np.average(sf, axis=1)
        blockaverage = np.average(blocks, axis=0)
        blockstd = np.std(blocks, axis=0)
        return blockaverage, blockstd
        
    
    
