import numpy as np
cimport numpy as np
import cython

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function

def density_help(np.ndarray[np.int64_t,   ndim=1] indices,
                 np.ndarray[np.float64_t, ndim=1] masses,
                 np.ndarray[np.float64_t, ndim=1] density):
    
    cdef int i = 0
    cdef int s = masses.shape[0]
    
    for i in range(0, s):
        density[indices[i]] += masses[i]

    return density


