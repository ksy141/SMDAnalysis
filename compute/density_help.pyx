import numpy as np
cimport numpy as cnp
import cython

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function

def density_help(cnp.ndarray[long, ndim=1] indices,
                 cnp.ndarray[double, ndim=1] masses,
                 cnp.ndarray[double, ndim=1] density):
    
    cdef int i = 0
    cdef long index
    cdef int s = masses.shape[0]
    
    for i in range(0, s):
        density[indices[i]] += masses[i]

    return density


