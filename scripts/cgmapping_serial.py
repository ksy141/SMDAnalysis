from MDAnalysis.coordinates.TRR import TRRWriter
import MDAnalysis as mda
import numpy as np

class CGMappingSerial():
    """
    Map all-atom trajectories to CG trajectories
    mappings['POPC']['PO4'] = ['P', 'O11', 'O12', 'O13', 'O14']
    mappings['1_TRP']['CA'] = ['CA']
    mappings['1_TRP']['SC'] = ['CB', ...]

    >>> cg  = smda.CGMappingSerial(mappings)
    >>> uCG = cg.run(u)
    >>> cg.write('CGTraj', uCG)
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
        uCG = self.create_universe(u)
        AAma, CGma, N_res = self.ma(u, uCG)
        
        def calculate(AAag, N):
            p = np.split(AAag.positions, N)
            f = np.split(AAag.forces, N)
            w = np.split(AAag.masses, N)[0]
            
            XCG = np.average(p, weights=w, axis=1)
            FCG = np.sum(f, axis=1)
            return XCG, FCG

        for i, ts in enumerate(u.trajectory):
            if i%100 == 0:
                print("processing %d/%d frames" %(i, len(u.trajectory)))
            uCG.trajectory[i].dimensions = u.trajectory[i].dimensions

            for resname in self.mappings.keys():
                for atn in self.mappings[resname].keys():
                    XCG, FCG = calculate(AAma[resname][atn], N_res[resname])
                    CGma[resname][atn].positions = XCG
                    CGma[resname][atn].forces    = FCG
        return uCG


    def ma(self, u, uCG):
        """
        Parameters
        ----------
        u:   AA universe
        uCG: CG universe
        
        Returns
        -------
        AAma, CGma, N_res

        >>> mappings = {'MeO': {'CH3': ['C', 'HA', 'HB', 'HC'],
        ...                     'OH':  ['O', 'H']}}

        >>> AAma = {'MeO': {'CH3': u.select_atoms('resname MeO and name C HA HB HC'),
        ...                 'OH':  u.select_atoms('resname MeO and name O H')}}

        >>> CGma = {'MeO': {'CH3': uCG.select_atoms('resname MeO and name CH3'),
        ...                 'OH':  uCG.select_atoms('resname MeO and name OH')}}
 
        >>> N_res = {'MeO':   u.select_atoms('resname MeO').n_residues =
        ...                 uCG.select_atoms('resname MeO').n_residues}
        """

        ### mappings to mapping AtomGroups  
        AAma = {}
        for resname in self.mappings.keys():
            AAma[resname] = {}
            sresname = resname.split('_')
            
            if len(sresname) == 1: # mappings['POPC']['NC3']
                for atn in self.mappings[resname].keys():
                    txt = 'resname ' + resname
                    txt += ' and name '
                    txt += ' '.join(self.mappings[resname][atn])
                    AAma[resname][atn] = u.select_atoms(txt)

            else: # mappings['278_TRP']['CA']
                for atn in self.mappings[resname].keys():
                    txt =  'resname ' + sresname[1]
                    txt += ' and resid ' + sresname[0]
                    txt += ' and name '
                    txt += ' '.join(self.mappings[resname][atn])
                    AAma[resname][atn] = u.select_atoms(txt)

        
        ### cg mappings to cg ma
        CGma = {}
        for resname in self.mappings.keys():
            CGma[resname] = {}

            for atn in self.mappings[resname].keys():
                sresname = resname.split('_')
                if len(sresname) == 1:
                    txt = 'resname %s and name %s' %(resname, atn)
                    CGma[resname][atn] = uCG.select_atoms(txt)
                else:
                    txt = 'resname %s and resid %s and name %s' %(sresname[1], sresname[0], atn) 
                    CGma[resname][atn] = uCG.select_atoms(txt)


        ### N
        N_res = {}
        for resname in self.mappings.keys():
            sresname = resname.split('_')
            if len(sresname) == 1:
                txt = 'resname %s' %resname
            else:
                txt = 'resname %s and resid %s' %(sresname[1], sresname[0])
            N_res[resname] = u.select_atoms(txt).n_residues

        return AAma, CGma, N_res


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

        self.analyze_universe(u)

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
        uCG.add_TopologyAttr('type', attr['name'])
        uCG.add_TopologyAttr('mass', attr['mass'])

        fac = np.zeros((nframes, n_atoms, 3))
        uCG.load_new(fac, forces=fac, order='fac')
        self.analyze_universe(uCG)
        
        uCG.trajectory[0].dt = u.trajectory[0].dt
        for i, ts in enumerate(u.trajectory):
            uCG.trajectory[i].dimensions = u.trajectory[i].dimensions
     
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
    
    
    def name2type(self, universe, types={}):
        """ 
        When types is given: Use for - if loop.
        >>> types = {'SOL': 1, 'SOD': 2, 'CLA': 3}

        When types is not given:
        It will assign the number based on atom name.
        Use for - if - else loop.

        When types is given for some atoms:
        Problematic!
        """

        t = 1
        for atom in universe.atoms:
            if atom.name in types.keys():
                atom.type = types[atom.name]
            else:
                atom.type = t
                types[atom.name] = t
                t += 1
        print(types)


    def write(self, fname, universe, units=None):
        universe.trajectory[-1]
        universe.atoms.write(fname + '.gro')
        
        trr = TRRWriter(fname + '.trr', universe.atoms.n_atoms)
        if units:
            trr.units = units

        for ts in universe.trajectory:
            trr.write(ts)
        trr.close()


