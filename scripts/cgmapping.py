import MDAnalysis as mda
import numpy as np

class CGMapping:
    """Map all-atom trajectories to CG trajectories
    Input mapping array and type (int number)

    Example
    -------
    >>> mappings = {'MeO': {'CH3': ['C', 'HA', 'HB', 'HC'],
    ...                     'OH':  ['O', 'H']}}

    >>> names2types = {'CH3': 'CG1', 'OH': 'CG2'}

    >>> ma = {'MeO': {'CH3': u.select_atoms('resname MeO and name C HA HB HC'),
    ...               'OH':  u.select_atoms('resname MeO and name O H')}}

    >>> u = mda.Universe('methanol_1728.pdb', 'traj.whole.trr')
    >>> uCG = cging.run(u)
    
    Make bonds
    CH3 = [0, 2, 4, ...] (CG atom indices of CH3)
    OH  = [1, 3, 5, ...] (CG atom indices of OH)
    np.concatenate([CH3, OH]) = [0, 2, 4, ..., 1, 3, 5...]
    np.concatenate([CH3, OH]).reshape(2, -1) = [[0, 2, 4, ...], 
                                                [1, 3, 5, ...]]
    np.concatenate([CH3, OH]).reshape(2, -1).T = [[0, 1], 
                                                  [2, 3],
                                                  [4, 5], ...]

    >>> CH3   = uCG.select_atoms('resname MeO and name CH3').indices 
    >>> OH    = uCG.select_atoms('resname MeO and name OH').indices
    >>> bonds = np.concatenate([CH3, OH]).reshape(2, -1).T 
    >>> bonds = tuple(map(tuple, bonds)) 
    >>> uCG.add_TopologyAttr('bonds', bonds)

    >>> cgatoms2nums = CGMapping().cgatoms2nums
    >>> cgma         = CGMapping().cgma

    cgatoms2nums = {'CH3': '1', 'OH': '2'}
    cgma         = {'MeO': {'CH3': uCG.select_atoms('resname MeO and name CH3'),
                            'OH':  uCG.select_atoms('resname MeO and name OH')}}
    """

    def __init__(self, mappings=None, names2types=None):
        """Set up a CGMapping
        Parameters
        ----------
        mappings : {'resname': {'CGsite': ['AA atom names']}}
        names2types : {'First CG site': 'CG1', 'Second CG site': 'CG2', ...}
        """

        self.mappings = mappings
        self.names2types = names2types
        pass

    def run(self, u):
        self.analyze_universe(u)
        
        ### mappings to mapping AtomGroups  
        ma = {}
        for resname in self.mappings.keys():
            ma[resname] = {}
            
            for atn in self.mappings[resname].keys():
                txt = 'resname ' + resname
                txt += ' and name '
                txt += ' '.join(self.mappings[resname][atn])
                print(resname + '+' + atn + ': ' + txt)
                ma[resname][atn] = u.select_atoms(txt)

       
        ### Create blank Universe
        uCG = self.create_universe(ma)
        
        ### Put data into blank universe
        for i, ts in enumerate(u.trajectory):
            if i%100 == 0:
                print("processing %d/%d frames" %(i,len(u.trajectory)))
            uCG.trajectory[i].dimensions = u.dimensions

            XCG = np.empty((0, 3))
            FCG = np.empty((0, 3))

            for resname in ma.keys():
                X = []
                F = []

                for ag in ma[resname].values():
                    a, b = np.unique(ag.names, return_counts = True)
                    natoms = len(a)

                    w = ag.masses[:natoms]
                    X1 = self._block_avg(ag.positions, w, natoms)
                    F1 = self._block_sum(ag.forces, natoms)
                
                    X.append(X1)
                    F.append(F1)
        
                XCG = np.append(XCG, self._mix(X), axis=0)
                FCG = np.append(FCG, self._mix(F), axis=0)
        
            uCG.trajectory[i].positions = XCG
            uCG.trajectory[i].forces = FCG

        self.analyze_universe(uCG)
 

        ### cgatoms to cg atom numbers
        cgatoms2nums = {}
        for resname in self.mappings.keys():
            #cgatoms2nums[resname] = {}
            n = 1

            for atn in self.mappings[resname].keys():
                #cgatoms2nums[resname][atn] = str(n)
                cgatoms2nums[atn] = str(n)
                n += 1

        ### cg mappings to cg ma
        cgma = {}
        for resname in self.mappings.keys():
            cgma[resname] = {}

            for atn in self.mappings[resname].keys():
                cgma[resname][atn] = uCG.select_atoms('resname %s and name %s' %(resname, atn))
        
        self.cgatoms2nums = cgatoms2nums
        self.cgma         = cgma
        self.uCG          = uCG
        
        return uCG
    

    def create_universe(self, ma):
        nres_all   = 0
        natoms_all = 0
        last_resid = 0
        
        attrs = {'resnames': [],
                 'resids':   [],
                 'names':    [],
                 'masses':   [],
                 'types':    []}
    
        for resname in ma.keys():
            ws = []
            names = []
            resnames = []
            resids = []
            types  = []
    
            for cgname, ag in ma[resname].items():
                a, b = np.unique(ag.names, return_counts = True)

                natoms = len(a)            
                nres = b[0]
    
                natoms_all += nres
    
                ws.append(self._block_1d_sum(ag.masses, natoms))
                names.append([cgname] * nres)
                types.append([self.names2types[cgname][2:]] * nres)
                resnames = [resname] * nres
                resids.append(np.arange(last_resid, last_resid + nres))
            
            nframes = len(ag.universe.trajectory)
            last_resid += nres    
            nres_all += nres
    
            attrs['resnames'].append(resnames)
            attrs['resids'].append(  self._mix_1d(resids))
            attrs['names'].append(   self._mix_1d(names))
            attrs['masses'].append(  self._mix_1d(ws))
            attrs['types'].append(   self._mix_1d(types))
    
        
        for key, value in attrs.items():
            attrs[key] = [item for sublist in value for item in sublist]
    
        uCG = mda.Universe.empty(n_atoms = natoms_all,
                                 n_residues = nres_all,
                                 atom_resindex = attrs['resids'],
                                 trajectory = True,
                                 forces = True)
        
        uCG.transfer_to_memory()
        uCG.add_TopologyAttr('masses')
        uCG.add_TopologyAttr('resnames')
        uCG.add_TopologyAttr('names')
        uCG.add_TopologyAttr('resids')
        uCG.add_TopologyAttr('types')
    
    
        ag = uCG.select_atoms('all')
        ag.masses = attrs['masses']
        ag.residues.resnames = attrs['resnames']
        ag.names  = attrs['names']
        ag.types  = attrs['types']
        
        fac = np.zeros((nframes, natoms_all, 3))
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
    
    
    def _mix_1d(self, arrays):
        """Parameters
        ----------
        arrays : an array of 1D arrays; lengths of 1D arrays are the same. 

        Examples
        --------
        >>> arrays = [[A, A, A], [B, B, B], [C, C, C]]
        >>> _mix_1d(arrays) = [A, B, C, A, B, C, A, B, C]
        """

        narray = len(arrays)
        nlen   = len(arrays[0])
        out = []
    
        for i in range(nlen):
            for j in range(narray):
                out.append(arrays[j][i])
    
        return np.array(out)

    
    def _mix(self, arrays):
        """Mix CG positions/forces
        Examples CG sites of CH3 and OH = CG atoms of CH3 and OH
        -------------------------------
        arrays = [np.array(
                  [[px(CH3_1), py(CH3_1), pz(CH3_1)],
                   [px(CH3_2), py(CH3_2), pz(CH3_2)],
                   ...]),

                  np.array(
                  [[px(OH1), py(OH1), pz(OH1)],
                   [px(OH2), py(OH2), pz(OH2)],
                   ...])]
        
        np.concatenate(arrays, axis=1) =
        [[px(CH3_1), py(CH3_1), pz(CH3_1), px(OH1), py(OH1), pz(OH1)],
         [px(CH3_2), py(CH3_2), pz(CH3_2), px(OH2), py(OH2), pz(OH2)],
         ...]

        _mix(arrays) = 
        [[px(CH3_1), py(CH3_1), pz(CH3_1)],
         [px(OH1),   py(OH1),   pz(OH1)], 
         [px(CH3_2), py(CH3_2), pz(CH3_2)],
         [px(OH2),   py(OH2),   pz(OH2)],
         ...]
        """

        dim = arrays[0].shape[1] #3
        out = np.concatenate(arrays, axis=1).reshape(-1, dim)
        return out
    
    
   
    def _block_avg(self, f, w, natoms):
        """Map AA positions to CG positions for each CG type for each molecule type
        Examples CG site OH = AA atoms of O and H
        -------------------
        natoms for OH = 2
        w for OH = [16, 1]
        
        arrays for OH = 
        [[px(O1), py(O1), pz(O1)], 
         [px(H1), py(H1), pz(H1)],
         [px(O2), py(O2), pz(O2)],
         [px(H2), py(H2), pz(H2)], ...]

        x  = [px(O1), px(H1), px(O2), px(H2), ...]
        
        x.reshape(-1, natoms) = 
        [[px(O1), px(H1)],
         [px(O2), px(H2)], ...]
        
        np.average(x.reshape(-1, natoms), weights=w, axis=1) = 
        [16/17 * fx(O1) + 1/17 * fx(H1), 16/17 * fx(O2) + 1/17 * fx(H2), ...]
        """

        x = f[:,0]
        y = f[:,1]
        z = f[:,2]
    
        xt = np.average(x.reshape(-1, natoms), weights=w, axis=1)
        yt = np.average(y.reshape(-1, natoms), weights=w, axis=1)
        zt = np.average(z.reshape(-1, natoms), weights=w, axis=1)
    
        return np.transpose([xt, yt, zt])
    
    def _block_sum(self, f, natoms):
        """Map AA forces to CG forces for each CG type for each molecule type
        Examples CG site OH = AA atoms of O and H
        -------------------
        natoms for OH = 2
        
        arrays for OH = 
        [[fx(O1), fy(O1), fz(O1)], 
         [fx(H1), fy(H1), fz(H1)],
         [fx(O2), fy(O2), fz(O2)],
         [fx(H2), fy(H2), fz(H2)], ...]

        x  = [fx(O1), fx(H1), fx(O2), fx(H2), ...]
        
        x.reshape(-1, natoms) = 
        [[fx(O1), fx(H1)],
         [fx(O2), fx(H2)], ...]
        
        np.sum(x.reshape(-1, natoms), axis=1) = 
        [fx(O1) + fx(H1), fx(O2) + fx(H2), ...]
        """

        x = f[:,0]
        y = f[:,1]
        z = f[:,2]
    
        xt = np.sum(x.reshape(-1, natoms), axis=1)
        yt = np.sum(y.reshape(-1, natoms), axis=1)
        zt = np.sum(z.reshape(-1, natoms), axis=1)
    
        return np.transpose([xt, yt, zt])
    
    
    def _block_1d_sum(self, f, natoms):
        f = np.array(f)
        return np.sum(f.reshape(-1, natoms), axis=1)


    def bond_sort(self, blist):
        """ sort bonds so that the bonds' order is the same with data file.
        First sort rows and then columns
        
        blists = {'DOPC': [['HG', 'MG'], ['MG', 'MT2'], ['MG', 'MT3'], ['MT2', 'ET2'], ['MT3', 'ET3']]}
        blist  = blists['DOPC'] = [['HG', 'MG'], ['MG', 'MT2'], ['MG', 'MT3'], ['MT2', 'ET2'], ['MT3', 'ET3']]
        
        cgatoms2nums = {'HG': '1', 'MG': '2', 'MT2': '3', 'ET2': '4', 'MT3': '5', 'ET3': '6'}
        vfunc(blist) =            [['1',  '2'],  ['2',  '3'],   ['2',  '4'],   ['3',    '4'],  ['5', '6']]
        """
        
        func  = lambda t: self.cgatoms2nums[t]
        vfunc = np.vectorize(func)

        t = []
        lastaxis = vfunc(blist).astype(int)
        for (i, j), (at1, at2) in zip(lastaxis, blist):
            if i > j:
                t.append([at2, at1])
            else:
                t.append([at1, at2])

        t = np.array(t)
        s = vfunc(t).astype(int)
        return t[np.lexsort(s[:, ::-1].T)]


    def angle_sort(self, alist):
        """ sort angles so that the angles' order is the same with data file.
        First sort rows (compare 0 and 2 elements in each row) 
        and then columns (based on the 0th index)

        alists = {'DOPC': [['HG',  'MG', 'MT2'], 
                           ['MG', 'MT2', 'ET2'], 
                           ['HG',  'MG', 'MT3'], 
                           ['MG', 'MT3', 'ET3'], 
                           ['MT2', 'MG', 'MT3']]}
        
        alist = alists['DOPC']
        vfunc(alist) = [['1', '2', '3'], 
                        ['2', '3', '4'],
                        ['1', '2', '5'],
                        ['2', '5', '6'],
                        ['3', '2', '4']]
        """
        
        func  = lambda t: self.cgatoms2nums[t]
        vfunc = np.vectorize(func)

        t = []
        lastaxis = vfunc(alist).astype(int)
        for (i, j, k), (at1, at2, at3) in zip(lastaxis, alist):
            if i > k:
                t.append([at3, at2, at1])
            else:
                t.append([at1, at2, at3])
    
        t = np.array(t)
        s = vfunc(t).astype(int)
        return t[np.lexsort(s[:, ::-1].T)]
        

    def write_sys(self, fname, blists, alists):
        """names2types = {'HG': '1',  'MG': '2', 
                          'MT2': '3', 'ET2': '4', 
                          'MT3': '3', 'ET3': '4'}
        list(set(names2types.values())) = ['1', '2', '3', '4']
        """

        def type_order(str1, str2):
            ret = None
            if int(str1[2:]) > int(str2[2:]):
                ret = (str2, str1)
            else:
                ret = (str1, str2)
            return ret

        ofile = open(fname, 'w')
        v = sorted(list(set(self.names2types.values())))
        
        ### Pair coeff
        for t0 in v:
            i = int(t0[2:])
            for t1 in v:
                j = int(t1[2:])
                if not i > j:
                    ofile.write('pair_coeff {i:d} {j:d} {it:s}_{jt:s}.table {it:s}_{jt:s}\n'.format(
                        i=i, j=j, it=t0, jt=t1))
            ofile.write('\n')
        ofile.write('\n')
        

        ### Bond coeff
        saved_blists = []
        i = 0
        
        for blist in blists.values():
            for b in blist:
                t0, t1 = type_order(self.names2types[b[0]], self.names2types[b[1]])
                if not [t0, t1] in saved_blists:
                    i += 1
                    saved_blists.append([t0, t1])
                    ofile.write('bond_coeff {i:d} {t0:s}_{t1:s}_bon.table {t0:s}_{t1:s}_bon\n'.format(
                        i=i, t0=t0, t1=t1))
        ofile.write('\n')


        ### Angle coeff
        saved_alists = []
        i = 0

        for alist in alists.values():
            for a in alist:
                t1 = self.names2types[a[1]]
                t0, t2 = type_order(self.names2types[a[0]], self.names2types[a[2]])
                if not [t0, t1, t2] in saved_alists:
                    i += 1
                    saved_alists.append([t0, t1, t2])
                    ofile.write('angle_coeff {i:d} {t1:s}_{t0:s}_{t2:s}_ang.table {t1:s}_{t0:s}_{t2:s}_ang\n'.format(
                        i=i, t0=t0, t1=t1, t2=t2))
            
        ofile.write('\n')
        ofile.close()


    def write_top(self, fname, blists, alists):
        ofile = open(fname, 'w')
        ofile.write('cgsites %d\n' %(self.uCG.atoms.n_atoms))
        
        v = list(self.names2types.values())
        cgtypes = 'cgtypes %d\n' %(len(set(v)))
        for typename in sorted(set(v), key=lambda x: v.index(x)):
            cgtypes += typename + '\n'
        
        ofile.write(cgtypes)
        ofile.write('moltypes %d\n' %(len(self.mappings.keys())))
        
        for resname in self.mappings.keys():
            mol = 'mol %d -1\nsitetypes\n' %(len(self.mappings[resname].keys()))
            for atomname in self.mappings[resname].keys():
                mol += self.names2types[atomname][2:] + '\n'
            ofile.write(mol)
            
            bonds_txt = 'bonds %d\n' %(len(blists[resname]))
            for b0, b1 in blists[resname]:
                bonds_txt += '{} {}\n'.format(self.cgatoms2nums[b0],
                                              self.cgatoms2nums[b1])
            ofile.write(bonds_txt)
            
            angs_txt = 'angles %d 0\n' %(len(alists[resname]))
            for a0, a1, a2 in alists[resname]:
                angs_txt += '{} {} {}\n'.format(self.cgatoms2nums[a1],
                                                self.cgatoms2nums[a0],
                                                self.cgatoms2nums[a2])
            ofile.write(angs_txt)
            ofile.write('dihedrals 0 0\n')
        
        # SYSTEM
        systems = 'system %d\n' %(len(self.mappings.keys()))
        for Mt, resname in enumerate(self.mappings.keys(), 1):
            Mn = len(self.cgma[resname][list(self.mappings[resname].keys())[0]])
            systems += '%d %d\n' %(Mt, Mn)
        ofile.write(systems)
        
        ofile.close()

