from __future__ import absolute_import, division, print_function
from ..common.traj_to_numpy import TrajectoryToNumpy
from MDAnalysis import Universe
import numpy as np

class OrderParameters:
    def __init__(self):
        pass
    
    def compute_OP(self, u, atoms = [], resname = False, add_sel = False, sel_update = False, nblocks = 5):
        assert isinstance(u, Universe)
        numframes = len(u.trajectory)
        size = numframes/nblocks

        if len(atoms) < 2:
            print('please provide more than one atoms')
        
        if not resname:
            print('please provide resname')
         
        natoms = len(atoms)
        center = atoms[0]

        g1 = "resname %s and name %s" %(resname, center)
        g2 = "resname %s and (name " %resname + " ".join(atoms[1:]) + ")"

        if add_sel:
            g1 += " and %s" %add_sel
            g2 += " and %s" %add_sel
        
        if sel_update:
            group1 = u.select_atoms(g1, updating=True)
            group2 = u.select_atoms(g2, updating=True)
        else:
            group1 = u.select_atoms(g1)
            group2 = u.select_atoms(g2)
        
        
        blocks = []
        for block in range(nblocks):
            Scd = []
            for ts in u.trajectory[int(block*size):int((block+1)*size)]:
                p1 = np.repeat(group1.positions, natoms-1, axis=0)
                p2 = group2.positions
                dp = p2 - p1
                norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
                cos_theta = dp[...,2]/norm
                S = -0.5 * (3 * np.square(cos_theta) -1 )
                order_param = np.average(S)
                Scd.append(order_param)
            blocks.append(np.average(Scd))
        
        blockaverage = np.average(blocks) 
        blockstd = np.std(blocks)

        return blockaverage, blockstd



POPC1 = [['C31', 'O31', 'O32'],
         ['C32', 'H2X', 'H2Y'], 
         ['C33', 'H3X', 'H3Y'], 
         ['C34', 'H4X', 'H4Y'], 
         ['C35', 'H5X', 'H5Y'], 
         ['C36', 'H6X', 'H6Y'], 
         ['C37', 'H7X', 'H7Y'], 
         ['C38', 'H8X', 'H8Y'], 
         ['C39', 'H9X', 'H9Y'],
         ['C310', 'H10X', 'H10Y'], 
         ['C311', 'H11X', 'H11Y'],
         ['C312', 'H12X', 'H12Y'],
         ['C313', 'H13X', 'H13Y'],
         ['C314', 'H14X', 'H14Y'],
         ['C315', 'H15X', 'H15Y'],
         ['C316', 'H16X', 'H16Y', 'H16Z']]
 
POPC2 = [['C21', 'O21', 'O22'],
         ['C22', 'H2R', 'H2S'], 
         ['C23', 'H3R', 'H3S'], 
         ['C24', 'H4R', 'H4S'], 
         ['C25', 'H5R', 'H5S'], 
         ['C26', 'H6R', 'H6S'], 
         ['C27', 'H7R', 'H7S'], 
         ['C28', 'H8R', 'H8S'], 
         ['C29', 'H91'],
         ['C210', 'H101'], 
         ['C211', 'H11R', 'H11S'],
         ['C212', 'H12R', 'H12S'],
         ['C213', 'H13R', 'H13S'],
         ['C214', 'H14R', 'H14S'],
         ['C215', 'H15R', 'H15S'],
         ['C216', 'H16R', 'H16S'],
         ['C217', 'H17R', 'H17S'],
         ['C218', 'H18R', 'H18S', 'H18T']]


DOPE1 = [['C31', 'O31', 'O32'],
         ['C32', 'H2X', 'H2Y'], 
         ['C33', 'H3X', 'H3Y'], 
         ['C34', 'H4X', 'H4Y'], 
         ['C35', 'H5X', 'H5Y'], 
         ['C36', 'H6X', 'H6Y'], 
         ['C37', 'H7X', 'H7Y'], 
         ['C38', 'H8X', 'H8Y'], 
         ['C39', 'H9X'],
         ['C310', 'H10X'], 
         ['C311', 'H11X', 'H11Y'],
         ['C312', 'H12X', 'H12Y'],
         ['C313', 'H13X', 'H13Y'],
         ['C314', 'H14X', 'H14Y'],
         ['C315', 'H15X', 'H15Y'],
         ['C316', 'H16X', 'H16Y'], 
         ['C317', 'H17X', 'H17Y'],
         ['C318', 'H18X', 'H18Y', 'H18Z']]
 
DOPE2 = [['C21', 'O21', 'O22'],
         ['C22', 'H2R', 'H2S'], 
         ['C23', 'H3R', 'H3S'], 
         ['C24', 'H4R', 'H4S'], 
         ['C25', 'H5R', 'H5S'], 
         ['C26', 'H6R', 'H6S'], 
         ['C27', 'H7R', 'H7S'], 
         ['C28', 'H8R', 'H8S'], 
         ['C29', 'H9R'],
         ['C210', 'H10R'], 
         ['C211', 'H11R', 'H11S'],
         ['C212', 'H12R', 'H12S'],
         ['C213', 'H13R', 'H13S'],
         ['C214', 'H14R', 'H14S'],
         ['C215', 'H15R', 'H15S'],
         ['C216', 'H16R', 'H16S'],
         ['C217', 'H17R', 'H17S'],
         ['C218', 'H18R', 'H18S', 'H18T']]


SAPI1 = [['C31', 'O31', 'O32'],
         ['C32', 'H2X', 'H2Y'], 
         ['C33', 'H3X', 'H3Y'], 
         ['C34', 'H4X', 'H4Y'], 
         ['C35', 'H5X', 'H5Y'], 
         ['C36', 'H6X', 'H6Y'], 
         ['C37', 'H7X', 'H7Y'], 
         ['C38', 'H8X', 'H8Y'], 
         ['C39', 'H9X', 'H9Y'],
         ['C310', 'H10X', 'H10Y'], 
         ['C311', 'H11X', 'H11Y'],
         ['C312', 'H12X', 'H12Y'],
         ['C313', 'H13X', 'H13Y'],
         ['C314', 'H14X', 'H14Y'],
         ['C315', 'H15X', 'H15Y'],
         ['C316', 'H16X', 'H16Y'],
         ['C317', 'H17X', 'H17Y'],
         ['C318', 'H18X', 'H18Y', 'H18Z']]

SAPI2 = [['C21', 'O21', 'O22'],
         ['C22', 'H2R', 'H2S'], 
         ['C23', 'H3R', 'H3S'], 
         ['C24', 'H4R', 'H4S'], 
         ['C25', 'H5R'], 
         ['C26', 'H6R'], 
         ['C27', 'H7R', 'H7S'], 
         ['C28', 'H8R'], 
         ['C29', 'H9R'],
         ['C210', 'H10R', 'H10S'], 
         ['C211', 'H11R'],
         ['C212', 'H12R'],
         ['C213', 'H13R', 'H13S'],
         ['C214', 'H14R'],
         ['C215', 'H15R'],
         ['C216', 'H16R', 'H16S'],
         ['C217', 'H17R', 'H17S'],
         ['C218', 'H18R', 'H18S'],
         ['C219', 'H19R', 'H19S'],
         ['C220', 'H20R', 'H20S', 'H20T']]

DPPC1 = [['C31', 'O31', 'O32'],
         ['C32', 'H2X', 'H2Y'], 
         ['C33', 'H3X', 'H3Y'], 
         ['C34', 'H4X', 'H4Y'], 
         ['C35', 'H5X', 'H5Y'], 
         ['C36', 'H6X', 'H6Y'], 
         ['C37', 'H7X', 'H7Y'], 
         ['C38', 'H8X', 'H8Y'], 
         ['C39', 'H9X', 'H9Y'],
         ['C310', 'H10X', 'H10Y'], 
         ['C311', 'H11X', 'H11Y'],
         ['C312', 'H12X', 'H12Y'],
         ['C313', 'H13X', 'H13Y'],
         ['C314', 'H14X', 'H14Y'],
         ['C315', 'H15X', 'H15Y'],
         ['C316', 'H16X', 'H16Y', 'H16Z']]

DPPC2 = [['C21', 'O21', 'O22'],
         ['C22', 'H2R', 'H2S'], 
         ['C23', 'H3R', 'H3S'], 
         ['C24', 'H4R', 'H4S'], 
         ['C25', 'H5R', 'H5S'], 
         ['C26', 'H6R', 'H6S'], 
         ['C27', 'H7R', 'H7S'], 
         ['C28', 'H8R', 'H8S'], 
         ['C29', 'H9R', 'H9S'],
         ['C210', 'H10R', 'H10S'], 
         ['C211', 'H11R', 'H11S'],
         ['C212', 'H12R', 'H12S'],
         ['C213', 'H13R', 'H13S'],
         ['C214', 'H14R', 'H14S'],
         ['C215', 'H15R', 'H15S'],
         ['C216', 'H16R', 'H16S', 'H16T']]

