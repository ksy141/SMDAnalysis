from pmda.parallel import ParallelAnalysisBase
from MDAnalysis.coordinates.TRR import TRRWriter
import MDAnalysis as mda
import numpy as np

class CGMapping():
    """
    >>> cging = smda.CGMapping()
    >>> ags, CGnames = cging.create_ags(u, mappings)
    >>> m = smda.CGMappingPMDA(ags)
    >>> m.run(n_jobs = 4)
    >>> uCG = cging.assign(m._results, CGnames, u, mappings)
    >>> cging.write('CGTraj', uCG)
    """
 
    def __init__(self):
        pass

    def create_universe(self, u, mappings):
        '''
        Parameters
        ----------
        mappings : {'resname': {'CGsite': ['AA atom names']}}
        
        >>> mappings['POPC']['PO4'] = ['P','O11','O12','O13','O14']
        >>> mappings['278_ALA']['CA'] = ['CA']
        >>> mappings['278_ALA']['SC'] = ['CB']
        
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

            if MEMBresname in mappings.keys():
                res = MEMBresname
            elif PROTresname in mappings.keys():
                res = PROTresname
            else:
                continue
                
            AAatoms = residue.atoms
            n_residues += 1
            for CGname, AAnames in mappings[res].items():
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

    
    def create_ags(self, u, mappings):
        AAma    = []
        CGnames = []
        for resname in mappings.keys():
            sresname = resname.split('_')
            
            if len(sresname) == 1: # mappings['POPC']['NC3']
                for atn in mappings[resname].keys():
                    txt = 'resname ' + resname
                    txt += ' and name '
                    txt += ' '.join(mappings[resname][atn])
                    AAma.append(u.select_atoms(txt))
                    CGnames.append(resname + '_' + atn)

            elif len(sresname) == 2: #mappings['278_TRP']['CA']
                for atn in mappings[resname].keys():
                    txt =  'resname ' + sresname[1]
                    txt += ' and resid ' + sresname[0]
                    txt += ' and name '
                    txt += ' '.join(mappings[resname][atn])
                    AAma.append(u.select_atoms(txt))
                    CGnames.append(resname + '_' + atn)

        return AAma, CGnames


    def assign(self, XF, CGnames, u, mappings):
        uCG = self.create_universe(u, mappings)
        XF = np.vstack(XF)
        assert len(XF) == uCG.trajectory.n_frames, 'frame number?'
        
        CGma = []
        for CGname in CGnames:
            s = CGname.split('_')
            
            if len(s) == 2: #'POPC_NC3'
                txt = 'resname ' + s[0]
                txt += ' and name ' + s[1]
                CGma.append(uCG.select_atoms(txt))

            elif len(s) == 3: # '278_TRP_CA'
                txt = 'resid ' + s[0]
                txt += ' and resname ' + s[1]
                txt += ' and name ' + s[2]
                CGma.append(uCG.select_atoms(txt))
            
        for frame, ts in enumerate(uCG.trajectory):
            for i, ag in enumerate(CGma):
                ag.positions = XF[frame][0][i]
                ag.forces    = XF[frame][1][i]
        
        return uCG

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

    
    def write(self, fname, universe, units=None):
        universe.trajectory[-1]
        universe.atoms.write(fname + '.gro')
        
        trr = TRRWriter(fname + '.trr', universe.atoms.n_atoms)
        if units:
            trr.units = units

        for ts in universe.trajectory:
            trr.write(ts)
        trr.close()


class CGMappingPMDA(ParallelAnalysisBase):
    """
    Assist CGMapping
    """

    def __init__(self, atomgroups):
        u = atomgroups[0].universe
        self.N = [ag.n_residues for ag in atomgroups]
        super(CGMappingPMDA, self).__init__(u, atomgroups)

    def _prepare(self):
        pass

    def _single_frame(self, ts, atomgroups):
        if (ts.time / 1000) % 10 == 0:
            print("processing %d ns" %(int(ts.time/1000)))

        XCG = []; FCG = []
        for ag, N in zip(atomgroups, self.N):
            p = np.split(ag.positions, N)
            f = np.split(ag.forces, N)
            w = np.split(ag.masses, N)[0]

            X = np.average(p, weights=w, axis=1)
            F = np.sum(f, axis=1)

            XCG.append(X)
            FCG.append(F)
        return [XCG, FCG]

    def _conclude(self):
        pass


