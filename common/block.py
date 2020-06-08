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
        sf = np.array_split(f, nblocks)
        blocks = [np.average(block, axis=0) for block in sf]
        
        blockavg = np.average(blocks, axis=0)
        blockstd = np.std(blocks, axis=0)
        
        return blockavg, blockstd
        
    
    
