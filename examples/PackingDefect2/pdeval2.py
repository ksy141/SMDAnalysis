import MDAnalysis as mda
import SMDAnalysis as smda
import os

u = mda.Universe('md.gro')

pd = smda.PackingDefect2()
lipid = '../../FF/CHARMM/toppar/top_all36_lipid.rtf'
TRIO  = '../../FF/CHARMM/toppar/TRIO.str'
SAPI  = '../../FF/CHARMM/toppar/toppar_all36_lipid_inositol.str'

radii = {'POPC': pd.read_top('POPC', lipid),
         'DOPE': pd.read_top('DOPE', lipid),
         'SAPI': pd.read_top('SAPI', SAPI),
         'TRIO': pd.read_top('TRIO', TRIO)}

if __name__ == '__main__':
    MEMB = u.select_atoms('resname POPC DOPE SAPI TRIO')
    pdPMDA = smda.PackingDefect2PMDA([MEMB], radii)
    #pdPMDA.prefix = 'md.'
    pdPMDA.N = 6000 # Change to a larger number for a big system
    pdPMDA.run()

