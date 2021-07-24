from   MDAnalysis.analysis import align
import numpy as np

class Fit:
    '''
    input: numpy position array that has a shape of (-1, 3)
    mob_sel: the selected position will be compared against ref_sel
    mob: the wanted position

    use:
    pos, rmsd = smda.Fit(mob_sel, ref_sel, mob)
    numpy array with a shape of (-1, 3)
    '''

    def __init__(self, mob_sel, ref_sel, mob):
        # Remove the center of geometry
        mob_sel0 = mob_sel - np.average(mob_sel, axis=0)
        ref_sel0 = ref_sel - np.average(ref_sel, axis=0)
        mob0     = mob     - np.average(mob_sel, axis=0)

        R, rmsd  = align.rotation_matrix(mob_sel0, ref_sel0)

        mob1     = np.matmul(R, mob0.T).T
        mob2     = mob1 + np.average(ref_sel, axis=0)

        return mob2, rmsd

