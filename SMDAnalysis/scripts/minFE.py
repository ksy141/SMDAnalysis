import MDAnalysis as mda
import SMDAnalysis as smda
from MDAnalysis.analysis import align

u = mda.Universe('../protein.whole.gro', '../protein.whole.xtc')
trp = u.select_atoms('bynum 281:304')
ala = u.select_atoms('bynum 924:933')
protein = u.select_atoms('protein')

ref = mda.Universe('../hairpin.ld.pdb')

odist = open('distance.dat', 'w')
ormsd  = open('rmsd.dat', 'w')

with mda.Writer('minFE.xtc', protein.n_atoms) as W:
    for ts in u.trajectory:
        if ts.frame % 10 == 0:
            print('frame: ', ts.frame)
        ### distance between W166-A203
        trpc = trp.center_of_geometry()
        alac = ala.center_of_geometry()
        pbc  = u.dimensions[0:3]
        d = smda.Distance(trpc, alac, pbc).distance(pbc=True)
        odist.write('{:5d} {:6.3f}\n'.format(ts.frame, d))
        
        if 6.65 < d < 7.1:
            ### align and then save
            selection = "protein and name CA and resid 15:57"
            old, new = align.alignto(u, ref, selection, weights="mass")
            ormsd.write('{:5d} {:6.3f}\n'.format(ts.frame, new))
            W.write(protein)
 
odist.close()
ormsd.close()


