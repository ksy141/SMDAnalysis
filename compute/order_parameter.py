from __future__ import absolute_import, division, print_function
from ..common.block import Block
from ..common.frame import Frame
from MDAnalysis import Universe
import numpy as np

class OrderParameters:
    def __init__(self):
        pass
    
    def compute_OP(self, u, atomlists = [], 
                   resname = False, add_sel = False, sel_update = False, 
                   nblocks = 5, b=0, e=1000000):

        assert isinstance(u, Universe), 'provide a proper universe'
        assert resname, 'provide a resname'
        
        ### Get the closest begin and end frames corresponding to b, e
        bframe, eframe = Frame().frame(u, b, e)
        C_numbers = []
        Cs = []
        Hs = []
        repeat = []
        for atoms in atomlists:
            C_number = atoms[0][2:]
            C_numbers.append(int(C_number))
            
            Cs.append(atoms[0])
            Hs.append(atoms[1:])
            repeat.append(len(atoms)-1) ## How many hydrogen atoms per center atom


        Hs_f = [item for sublist in Hs for item in sublist]
        assert int(np.sum(repeat)) == len(Hs_f), "wrong in repeats"
        ## Cs =   [C31, C32, C33, ..., C318]
        ## Hs =   [[O31, O32], [H2X, H2Y], [H3X, H3Y], ... [H18X, H18Y, H18Z]]
        ## Hs_f = [O31, O32, H2X, H2Y, H3X, ... , H18Z]
        ## repeat = [2, 2, 2, ..., 3]

        g1 = "resname %s" %resname + " and name " + " ".join(Cs)
        g2 = "resname %s" %resname + " and name " + " ".join(Hs_f)
        if add_sel:
            g1 += " and %s" %add_sel
            g2 += " and %s" %add_sel

        if sel_update:
            group1 = u.select_atoms(g1, updating=True)
            group2 = u.select_atoms(g2, updating=True)
 
        else:
            group1 = u.select_atoms(g1)
            group2 = u.select_atoms(g2)
        
        natoms       = len(Cs) ## How many center atoms
        nmols        = int(len(group1.positions)/natoms) ## How many molecules
        repeats      = repeat * nmols ## [2, 2, 2, ..., 3, | 2, 2, 2, ..., 3, | ...]
        splits       = np.cumsum(repeats)

        print('# of mols: %d' %nmols)
        print('# of Carbons per molecule: %d' %natoms)

        output = []
        for ts in u.trajectory[bframe:eframe]:
            if int(ts.time) % 100000 == 0:
                print("analyzing %d ns" %(int(ts.time/100000) * 100))
            p1 = np.repeat(group1.positions, repeats, axis=0)
            p2 = group2.positions
            dp = p2 - p1
            norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
            cos_theta = dp[...,2]/norm
            S = -0.5 * (3 * np.square(cos_theta) - 1)
            # out = np.split(S, splits)

            new_S = self._repeat(S, repeats)
            new_S.shape = (nmols, natoms)
            results = np.average(new_S, axis=0)
            output.append(results)
     
            #results = []
            #for mol in range(nmols):
            #    each_mol = []
            #    for i in range(natoms):
            #        index = mol*natoms + i
            #        each_mol.append(np.average(out[index]))
            #    results.append(each_mol)
            #output.append(np.average(results, axis=0))
        
        avg, std = Block().block(output, nblocks)

        return np.transpose([C_numbers, avg, std])



    def compute_OP_frame(self, pos1, pos2, n):
        p1 = np.repeat(pos1, n-1, axis=0) #carbon
        p2 = pos2 #hydrogens
        dp = p2 - p1
        norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
        cos_theta = dp[...,2]/norm
        S = -0.5 * (3 * np.square(cos_theta) - 1)
        order_param = np.average(S)
        return order_param
 
    
    def _repeat(self, x, reps):
        ''' 
        x = [1, 3, 5, 7, 9]
        reps = [2, 2, 1]
        out = [2, 6, 9]
        '''
        assert len(x) == int(np.sum(reps)), 'repeats wrong'

        i = 0
        out = []
        for rep in reps:
            tmp = []
            for r in range(rep):
                tmp.append(x[i])
                i += 1
            out.append(np.average(tmp))
    
        return np.array(out)


