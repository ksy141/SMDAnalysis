import MDAnalysis as mda
import numpy as np

class CGMapping:
    """
    Map all-atom trajectories to CG trajectories
    """

    def __init__(self, mappings=None):
        """
        Set up a CGMapping
        Parameters
        ----------
        mappings : {'resname': {'CGsite': ['AA atom names']}}
        
        >>> mappings['POPC']['PO4'] = ['P','O11','O12','O13','O14']
        >>> mappings['278_ALA']['CA'] = ['CA']
        >>> mappings['278_ALA']['SC'] = ['CB']
        """

        self.mappings = mappings


    def run(self, u):
        self.analyze_universe(u)
        uCG = self.create_universe(u)
        uCG.trajectory[0].dt = u.trajectory[0].dt
        self.analyze_universe(uCG)

        for i, ts in enumerate(u.trajectory):
            if i%100 == 0:
                print("processing %d/%d frames" %(i, len(u.trajectory)))
            uCG.trajectory[i].dimensions = u.trajectory[i].dimensions

            XCG = []
            FCG = []

            for residue in u.atoms.residues:
                MEMBresname = residue.resname
                PROTresname = str(residue.resid) + '_' + residue.resname
    
                if MEMBresname in self.mappings.keys():
                    res = MEMBresname
                elif PROTresname in self.mappings.keys():
                    res = PROTresname
                else:
                    continue
                    
                AAatoms = residue.atoms
                for CGname, AAnames in self.mappings[res].items():
                    bA = np.isin(AAatoms.names, AAnames)
                    selAAatoms = AAatoms[bA]
                    
                    pos = selAAatoms.positions
                    w   = selAAatoms.masses
                    X = np.average(pos, weights=w, axis=0)
                    XCG.append(X)

                    force = selAAatoms.forces
                    F = np.sum(force, axis=0)
                    FCG.append(F)

            uCG.trajectory[i].positions = XCG
            uCG.trajectory[i].forces    = FCG
        
        return uCG
        

    def create_universe(self, u):
        '''
        >>> n_residues = 1000
        >>> n_atoms    = 3 * 1000
        >>> resindices = [0,0,0,1,1,1,...,999,999,999]
        >>> resindices = np.repeat(np.arange(n_residues), 3)
        >>> segindices = [0] * n_residues
        
        len(resindices) = n_atoms
        len(segindices) = n_residues
        
        Note resindices MUST start from 0 and
        segindices is optional.

        >>> sol = mda.Universe.empty(n_atoms, 
        ...                          n_residues = n_residues,
        ...                          atom_resindex = resindices,
        ...                          residue_segindex = segindices,
        ...                          trajectory = True)

        >>> sol.add_TopologyAttr('resname', ['SOL'] * n_residues)
        >>> sol.add_TopologyAttr('resid',   np.arange(n_residues) + 1)
        >>> sol.add_TopologyAttr('name',    ['O','H1','H2'] * n_residues) 
        >>> sol.add_TopologyAttr('type',    ['O','H', 'H']  * n_residues)
        >>> sol.add_TopologyAttr('segid',   ['SOL'])

        len(resname) = n_residues
        len(resid)   = n_residues
        len(name)    = n_atoms
        len(type)    = n_atoms
        len(segid)   = n_segids
        '''

        try:
            u.atoms.forces
        except:
            for i, ts in enumerate(u.trajectory):
                u.trajectory[i].forces = np.zeros((u.atoms.n_atoms, 3))
        
        nframes = len(u.trajectory)
        

        attr = {'resname': [],
                'resid': [],
                'resindex': [],
                'name': [],
                'mass': []}
        
        n_residues = 0
        n_atoms = 0
        
        printed = []
        for residue in u.atoms.residues:
            # POPC for MEMBresname
            # 275_TRP for PROTresname
            MEMBresname = residue.resname
            PROTresname = str(residue.resid) + '_' + residue.resname 

            if MEMBresname in self.mappings.keys():
                res = MEMBresname
            elif PROTresname in self.mappings.keys():
                res = PROTresname
            else:
                continue
                
            AAatoms = residue.atoms
            n_residues += 1
            for CGname, AAnames in self.mappings[res].items():
                n_atoms += 1
                bA = np.isin(AAatoms.names, AAnames)

                comp = '%s+%s' %(res, CGname)
                if not comp in printed:
                    printed.append(comp)
                    print(comp + ': ' + ' '.join(AAatoms[bA].names))

                attr['resindex'].append(n_residues - 1)
                attr['name'].append(CGname)
                attr['mass'].append(np.sum(AAatoms[bA].masses))

            attr['resid'].append(residue.resid)
            attr['resname'].append(residue.resname)

        
        uCG = mda.Universe.empty(n_atoms = n_atoms,
                                 n_residues = n_residues,
                                 atom_resindex = attr['resindex'],
                                 trajectory = True,
                                 forces = True)
        
        uCG.add_TopologyAttr('resid', attr['resid'])
        uCG.add_TopologyAttr('resname', attr['resname'])
        uCG.add_TopologyAttr('name', attr['name'])
        uCG.add_TopologyAttr('mass', attr['mass'])

        fac = np.zeros((nframes, n_atoms, 3))
        uCG.load_new(fac, forces=fac, order='fac')
    
        return uCG


    def analyze_universe(self, universe):
        natoms  = universe.atoms.n_atoms
        nresids = len(universe.residues)
        nframes = len(universe.trajectory)
        print("\n---------------------------")
        print("UNIVERSE INFO")
        print("N_atoms:  {:d}".format(natoms))
        print("N_resids: {:d}".format(nresids))
        print("N_frames: {:d}".format(nframes))
        print("---------------------------\n")
    
    

