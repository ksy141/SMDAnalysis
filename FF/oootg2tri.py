import MDAnalysis as mda

u = mda.Universe('CHARMMGUI_TOP/step6.0_minimization.gro')
resname = 'TRIV'


### CHANGE THE RESIDUE NAME
ag = u.select_atoms('resname OOOTG')
ag.residues.resnames = resname


### CHANGE THE ORDER OF O12 and C11
O12  = u.select_atoms('resname ' + resname + ' and name O12')
O12p = O12.positions

C11  = u.select_atoms('resname ' + resname + ' and name C11')
C11p = C11.positions

O12.names = 'C11'
O12.positions = C11p

C11.names = 'O12'
C11.positions = O12p


### CHANGE PQO to ABC (chain-1)
for i in range(2, 19):
    ag = u.select_atoms('resname ' + resname + ' and name H' + str(i) + 'P')
    ag.names = 'H' + str(i) + 'A'

    ag = u.select_atoms('resname ' + resname + ' and name H' + str(i) + 'Q')
    ag.names = 'H' + str(i) + 'B'

    ag = u.select_atoms('resname ' + resname + ' and name H' + str(i) + 'O')
    ag.names = 'H' + str(i) + 'C'


### SAVE THE STRUCTURE
u.atoms.write('step6.0_minimization.gro')


