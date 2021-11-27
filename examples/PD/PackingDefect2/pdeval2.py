import MDAnalysis as mda
import SMDAnalysis as smda
import os

u = mda.Universe('../md.gro', '../md.xtc')

pd = smda.PackingDefect2()
lipid = '~/SMDAnalysis/FF/CHARMM/toppar/top_all36_lipid.rtf'
TRIO  = '~/SMDAnalysis/FF/CHARMM/toppar/trio.str'
SAPI  = '~/SMDAnalysis/FF/CHARMM/toppar/toppar_all36_lipid_inositol.str'

radii = {'POPC': pd.read_top('POPC', os.path.expanduser(lipid)),
         'DOPE': pd.read_top('DOPE', os.path.expanduser(lipid)),
         'SAPI': pd.read_top('SAPI', os.path.expanduser(SAPI)),
         'TRIO': pd.read_top('TRIO', os.path.expanduser(TRIO))}

if __name__ == '__main__':
    MEMB = u.select_atoms('resname POPC DOPE SAPI TRIO')
    pdPMDA = smda.PackingDefect2PMDA([MEMB], radii)
    pdPMDA.N = 5000 # Change to a larger number for a big system
    #pdPMDA.prob = False
    pdPMDA.run(n_jobs=-1)

