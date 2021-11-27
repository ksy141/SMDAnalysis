from pmda.parallel import ParallelAnalysisBase
import numpy as np
from .radii import types_radii
from MDAnalysis import Universe
import MDAnalysis as mda

class PackingDefectOld:
    def __init__(self):
        pass

    def read_top(self, resname, topology_file):
        """
        Examples
        --------
        pd = smda.PackingDefectOld()
        ff = os.getenv('HOME') + '/Dropbox/ff/charmm36.ff/'
        lipid = ff + 'top_all36_lipid.rtf'
        TRIO  = ff + 'TRIO.rtf'
        SAPI  = ff + 'stream/lipid/toppar_all36_lipid_inositol.str'

        radii = {'POPC': pd.read_top('POPC', lipid),
                 'DOPE': pd.read_top('DOPE', lipid),
                 'SAPI': pd.read_top('SAPI', SAPI),
                 'TRIO': pd.read_top('TRIO', TRIO)}

        Return
        ------
        radii['POPC']['P'] = [P radius, 'n']
        """

        tails = []
        tails += ['C2%d' %i for i in range(2, 23)]
        tails += ['C3%d' %i for i in range(2, 23)]
        tails += ['H%dR' %i for i in range(2, 23)]
        tails += ['H%dS' %i for i in range(2, 23)]
        tails += ['H%dX' %i for i in range(2, 23)]
        tails += ['H%dY' %i for i in range(2, 23)]
        tails += ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']
        
        top = open(topology_file)
        startread = False
        output = {}
        for line in top:
            if line.startswith('!'):
                continue
            if line.startswith('RESI {}'.format(resname)):
                startread = True
            if startread and line.startswith('BOND'):
                startread = False
            if startread and line.startswith('ATOM'):
                sl = line.split()
                if sl[1] in tails:
                    acyl = 'a'
                else:
                    acyl = 'n'
                output[sl[1]] = [types_radii[sl[2]], acyl]
        return output


    def defect_size(self, matrices, nbins=600, bin_max=150, 
            prob=True, fname='defect_histogram.dat'):

        rdf_settings = {'bins': nbins, 'range': (0, bin_max)}
        _, edges = np.histogram([-1], **rdf_settings)
        bins = 0.5 * (edges[1:] + edges[:-1])
    
        hist = np.zeros(nbins)
        for matrix in matrices:
            defects = []
            graph = self._make_graph(matrix)
            visited = set([])
            for n in graph:
                if n not in visited:
                    defect_loc = self._dfs(graph, n)
                    visited = visited.union(defect_loc)
                    #defects.append(len(defect_loc) * 0.01) #A2 to nm2
                    defects.append(len(defect_loc))
            
            hist += np.histogram(defects, **rdf_settings)[0]
        
        if np.sum(hist) == 0:
            return

        if prob:
            hist /= np.sum(hist)
        np.savetxt(fname, np.transpose([bins, hist]), fmt="%8.5f")


    def _dfs(self, graph, start):
        visited, stack = set(), [start]
        while stack:
            vertex = stack.pop()
            if vertex not in visited:
                visited.add(vertex)
                stack.extend(graph[vertex] - visited)
        return visited
    
    def _make_graph(self, matrix):
        graph = {}
        xis, yis = matrix.shape
    
        for (xi, yi), value in np.ndenumerate(matrix):
            if value == 0:
                continue
    
            n = xi * yis + yi
            nlist = []
    
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    # 8-neighbors
                    x = divmod(xi + dx, xis)[1]
                    y = divmod(yi + dy, yis)[1]
                    if matrix[x, y] == 1:
                        ndn = x * yis + y
                        nlist.append(ndn)
    
                   # 4 neighbors
                   # if dx * dy == 0:
                   #     x = divmod(xi + dx, xis)[1]
                   #     y = divmod(yi + dy, yis)[1]
                   #     if matrix[x, y] == 1:
                   #         ndn = x * yis + y
                   #         nlist.append(ndn)
    
            graph[n] = set(nlist) - set([n])
        return graph


C2 = ' '.join(['C2%d' %i for i in range(2, 23)])
C3 = ' '.join(['C3%d' %i for i in range(2, 23)])
class PackingDefectPMDAOld(ParallelAnalysisBase):
    def __init__(self, atomgroups, radii, nbins=600, bin_max=150):
        u = atomgroups[0].universe
        self.N  = 3000 #The maximum number of defects
        self.dt = u.trajectory[0].dt
        self.dx = 1
        self.dy = 1
        self.dz = 1
        self.zcut = 1
        self.radii = radii
        self.nbins = nbins
        self.bin_max = bin_max
        super(PackingDefectPMDAOld, self).__init__(u, atomgroups)

    def _prepare(self):
        pass

    def _single_frame(self, ts, atomgroups):
        if (ts.time / 1000) % 10 == 0:
            print("processing %d ns" %(int(ts.time/1000)))
        
        dim = ts.dimensions.copy()
        pbc = dim[0:3]
        print(pbc)
        pbc_xy0 = np.array([pbc[0],pbc[1],0])
        pbc_xyz = np.array([pbc[0],pbc[1],pbc[2]])
        #hz  = pbc[2]/2 #half z
        xarray = np.arange(0, pbc[0], self.dx)
        yarray = np.arange(0, pbc[1], self.dy)
        xx, yy = np.meshgrid(xarray, yarray)

        Mup = np.zeros_like(xx)
        Mdw = np.zeros_like(xx)

        PL = atomgroups[0]
        hz = np.average(PL.select_atoms('name P').positions[:,2])
        aa = PL.universe.atoms
        aa.positions -= pbc_xy0 * np.floor(aa.positions / pbc_xyz)

        zlimup = np.max(PL.positions[:,2])
        zlimdw = np.min(PL.positions[:,2])

        PL_up = PL.select_atoms('name P and prop z > %f' %hz).residues.atoms
        PL_dw = PL.select_atoms('name P and prop z < %f' %hz).residues.atoms
        umemb = PL_up.select_atoms('name ' + C2 + C3)
        lmemb = PL_dw.select_atoms('name ' + C2 + C3)
        utz   = np.average(umemb.positions[:,2])
        ltz   = np.average(lmemb.positions[:,2])
        
        for atom in PL:
            xatom, yatom, zatom = atom.position
            catom = atom.residue.atoms.select_atoms('name C2')[0].position[2]
            patom = atom.residue.atoms.select_atoms('name P')[0].position[2]
            
            if patom > hz: #up
                if zatom < catom - 6*self.zcut:
                    continue
                l = 'up'
                zarray = np.arange(catom-1, zlimup+1, self.dz)

            elif patom < hz: #dw
                if zatom > catom + 6*self.zcut:
                    continue
                l = 'dw'
                zarray = np.arange(catom+1, zlimdw-1, -self.dz)

            radius, acyl = self.radii[atom.resname][atom.name]
            if acyl == 'a':
                dM = 1e3
            elif acyl == 'n':
                dM = 1e6
 
            dxx = xx - xatom
            dxx -= pbc[0] * np.around(dxx/pbc[0])
                
            dyy = yy - yatom
            dyy -= pbc[1] * np.around(dyy/pbc[1])
            
            dzz = zarray - zatom
                
            # in order to improve efficiency
            # look around only the neighboring grids 
            # search distance limit = (1+r)**2
            dist_lim = (1 + radius)**2
            bx = dxx**2 < dist_lim 
            by = dyy**2 < dist_lim
            bz = dzz**2 < dist_lim

            new_xarray = xarray[bx[0]]
            new_yarray = yarray[by[:,0]]
            new_zarray = zarray[bz]
            new_marray = np.meshgrid(new_xarray, new_yarray, new_zarray)
            new_grid = np.vstack(new_marray).reshape(3,-1).T #grid.positions
            
            dist_meet = (np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)/2 + radius)**2
            dr2 = np.sum((new_grid - atom.position)**2, axis=1)
            br  = dr2 < dist_meet

            if np.sum(br) == 0:
                continue
            
            # new_grid[br] -> grid < dist_meet
            # new_grid[br][:,:-1] -> exclude the last column (z)
            # np.unique -> same x, y but different z
            ind = (np.unique(new_grid[br][:,:-1], axis=0)).astype(int)
            #ind = ind.astype(int)

            if l == 'up':
                Mup[ind[:,0], ind[:,1]] += dM
            elif l == 'dw':
                Mdw[ind[:,0], ind[:,1]] += dM
        
        if len(atomgroups) == 1:
            return Mup, Mdw, zlimup, zlimdw, dim
        

        ### TG
        TG = atomgroups[1]
        glycerols = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32',
                     'C1', 'C2', 'C3', 'C11', 'C21', 'C31',
                     'HA', 'HB', 'HS', 'HX', 'HY']
        
        # loop over TG atoms efficiently.. not all TG atoms.
        newTG = TG.select_atoms('prop z>%f or prop z<%f' %(utz-3, ltz+3))
        for atom in newTG:
            xatom, yatom, zatom = atom.position

            if zatom > utz-3:
                l = 'up'
                zarray = np.arange(utz, zlimup+1, self.dz)
            elif zatom < ltz+3:
                l = 'dw'
                zarray = np.arange(ltz, zlimdw-1, -self.dz)
            else:
                assert 1==0, 'TG atom?'
            
            radius, acyl = self.radii[atom.resname][atom.name]
            if atom.name in glycerols:
                dM = 1
            else:
                dM = 1e-3
 
            dxx = xx - xatom
            dxx -= pbc[0] * np.around(dxx/pbc[0])
                
            dyy = yy - yatom
            dyy -= pbc[1] * np.around(dyy/pbc[1])
            
            dzz = zarray - zatom
                
            # in order to improve efficiency
            # look around only the neighboring grids 
            # search distance limit = (1+r)**2
            dist_lim = (1 + radius)**2
            bx = dxx**2 < dist_lim 
            by = dyy**2 < dist_lim
            bz = dzz**2 < dist_lim

            new_xarray = xarray[bx[0]]
            new_yarray = yarray[by[:,0]]
            new_zarray = zarray[bz]
            new_marray = np.meshgrid(new_xarray, new_yarray, new_zarray)
            new_grid = np.vstack(new_marray).reshape(3,-1).T #grid.positions
            
            dist_meet = (np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)/2 + radius)**2
            dr2 = np.sum((new_grid - atom.position)**2, axis=1)
            br  = dr2 < dist_meet

            if np.sum(br) == 0:
                continue
            
            # new_grid[br] -> grid < dist_meet
            # new_grid[br][:,:-1] -> exclude the last column (z)
            # np.unique -> same x, y but different z
            ind = (np.unique(new_grid[br][:,:-1], axis=0)).astype(int)
            #ind = ind.astype(int)

            if l == 'up':
                Mup[ind[:,0], ind[:,1]] += dM
            elif l == 'dw':
                Mdw[ind[:,0], ind[:,1]] += dM
 
        return Mup, Mdw, zlimup, zlimdw, dim

    
    def _conclude(self):
        print("Concluding...")
        results = np.vstack(self._results)
        Mup    = results[:,0]
        Mdw    = results[:,1]
        zlimup = results[:,2].flatten()
        zlimdw = results[:,3].flatten()
        dim    = np.vstack(results[:,4])

        N  = self.N
        df = Universe.empty(n_atoms = N, 
                            n_residues = N, 
                            atom_resindex = np.arange(N), 
                            trajectory=True)

        df.add_TopologyAttr('resname', ['O'] * N)
        df.add_TopologyAttr('name',    ['O'] * N)
        df.add_TopologyAttr('resid', np.arange(N) + 1)

        nframes = len(dim)
        fac = np.zeros((nframes, N, 3))
        df.load_new(fac, order='fac')
        df.trajectory[0].dt = self.dt

        for i, ts in enumerate(df.trajectory):
            df.trajectory[i].dimensions = dim[i]

        PLacyl = df.copy()
        Deep   = df.copy()
        TGacyl = df.copy()
        TGglyc = df.copy()

        xxyy = []
        for i in range(len(df.trajectory)):
            xarray = np.arange(0, dim[i][0], self.dx)
            yarray = np.arange(0, dim[i][1], self.dy)
            xx, yy = np.meshgrid(xarray, yarray)
            xxyy.append([xx, yy])

        ### PL acyl
        PLacyls = []
        for i, ts in enumerate(PLacyl.trajectory):
            num = 0
            bA = (1e3 <= Mup[i]) & (Mup[i] < 1e6)
            PLacyls.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                PLacyl.atoms[num].position = np.array([y1, x1, zlimup[i]])
                num += 1
                
            bA = (1e3 <= Mdw[i]) & (Mdw[i] < 1e6)
            PLacyls.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                PLacyl.atoms[num].position = np.array([y1, x1, zlimdw[i]])
                num += 1

        ### Deep
        Deeps = []
        for i, ts in enumerate(Deep.trajectory):
            num = 0
            bA = (Mup[i] == 0)
            Deeps.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                Deep.atoms[num].position = np.array([y1, x1, zlimup[i]])
                num += 1
                
            bA = (Mdw[i] == 0)
            Deeps.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                Deep.atoms[num].position = np.array([y1, x1, zlimdw[i]])
                num += 1

        ### TG glycerol
        TGglycs = []
        for i, ts in enumerate(TGglyc.trajectory):
            num = 0
            bA = (1 <= Mup[i]) & (Mup[i] < 1e3)
            TGglycs.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                TGglyc.atoms[num].position = np.array([y1, x1, zlimup[i]])
                num += 1
                
            bA = (1 <= Mdw[i]) & (Mdw[i] < 1e3)
            TGglycs.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                TGglyc.atoms[num].position = np.array([y1, x1, zlimdw[i]])
                num += 1

        ### TG acyl
        TGacyls = []
        for i, ts in enumerate(TGacyl.trajectory):
            num = 0
            bA = (0 < Mup[i]) & (Mup[i] < 1)
            TGacyls.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                TGacyl.atoms[num].position = np.array([y1, x1, zlimup[i]])
                num += 1
                
            bA = (0 < Mdw[i]) & (Mdw[i] < 1)
            TGacyls.append(bA.astype(int))
            for x1, y1 in zip(xxyy[i][0][bA], xxyy[i][1][bA]):
                TGacyl.atoms[num].position = np.array([y1, x1, zlimdw[i]])
                num += 1


        def write(fname, u):
            u.trajectory[-1]
            u.atoms.write(fname + '.gro')
            with mda.Writer(fname + '.xtc', u.atoms.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(u.atoms)
        
        write('PLacyl', PLacyl)
        write('Deep',   Deep)
        write('TGglyc', TGglyc)
        write('TGacyl', TGacyl)

        PD = PackingDefectOld()
        PD.defect_size(PLacyls, fname='PLacyl.dat', nbins=self.nbins, bin_max=self.bin_max)
        PD.defect_size(Deeps,   fname='Deep.dat',   nbins=self.nbins, bin_max=self.bin_max)
        PD.defect_size(TGglycs, fname='TGglyc.dat', nbins=self.nbins, bin_max=self.bin_max)
        PD.defect_size(TGacyls, fname='TGacyl.dat', nbins=self.nbins, bin_max=self.bin_max)


