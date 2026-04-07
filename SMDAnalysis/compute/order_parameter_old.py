from __future__ import absolute_import, division, print_function
from ..common.traj_to_numpy import TrajectoryToNumpy
from ..common.block import Block
from ..common.frame import Frame
from MDAnalysis import Universe
import numpy as np

class OrderParameters:
    def __init__(self):
        pass
    
    def compute_OP(self, u, atoms = [], 
                   resname = False, add_sel = False, sel_update = False, 
                   nblocks = 5, b=0, e=1000000):

        assert isinstance(u, Universe)
        
        ### Get the closest begin and end frame corresponding to b, e
        bframe, eframe = Frame().frame(u, b, e)
        #print("frame starts at: %d" %bframe)
        #print("frame ends   at: %d" %eframe)
 
        natoms = len(atoms)
        center = atoms[0]
 
        if len(atoms) < 2:
            print('please provide more than one atoms')
        
        if not resname:
            print('please provide resname')
        
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

        Scd = []
        for ts in u.trajectory[bframe:eframe]:
            p1 = np.repeat(group1.positions, natoms-1, axis=0)
            p2 = group2.positions
            dp = p2 - p1
            norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
            cos_theta = dp[...,2]/norm
            S = -0.5 * (3 * np.square(cos_theta) -1 )
            order_param = np.average(S)
            Scd.append(order_param)
            
        average, std = Block().block(Scd, nblocks)

        return average, std




