import numpy as np

class Block:
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
        
    
    
