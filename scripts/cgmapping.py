import MDAnalysis as mda
import numpy as np

class CGMapping:
    def __init__(self, mappings=None, names2types=None):
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
                types.append([self.names2types[cgname]] * nres)
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
            attrs[key] = np.array(value).flatten()
    
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
        narray = len(arrays)
        nlen   = len(arrays[0])
        out = []
    
        for i in range(nlen):
            for j in range(narray):
                out.append(arrays[j][i])
    
        return np.array(out)
    
    
    def _mix(self, arrays):
        shape = arrays[0].shape[1]
        out = np.concatenate(arrays, axis=1).reshape(-1, shape)
        return out
    
    
   
    def _block_avg(self, f, w, natoms):
        x = f[:,0]
        y = f[:,1]
        z = f[:,2]
    
        xt = np.average(x.reshape(-1, natoms), weights=w, axis=1)
        yt = np.average(y.reshape(-1, natoms), weights=w, axis=1)
        zt = np.average(z.reshape(-1, natoms), weights=w, axis=1)
    
        return np.transpose([xt, yt, zt])
    
    def _block_sum(self, f, natoms):
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
    
    
    
