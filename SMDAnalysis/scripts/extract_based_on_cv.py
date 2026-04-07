import MDAnalysis as mda
import SMDAnalysis as smda

protein_index = "bynum 7 9 11 15 16 17"
u = mda.Universe("run.gro", "run.xtc")

protein = u.select_atoms("protein")
prot = u.select_atoms(protein_index)
memb = u.select_atoms("name P")

file = open('cv', 'w')
with mda.Writer('Desktop/interest.xtc', u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        smda.Wholemolecules(protein)
        dz = smda.COMDistance(memb, prot).scale_distance(pbc=True)[2]
        file.write("{: .3f} {: .3f}\n".format(u.trajectory.time, dz))

        if -0.1 < dz < 0.1:
            W.write(u.atoms)

file.close()
