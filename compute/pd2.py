from pmda.parallel import ParallelAnalysisBase
import numpy as np
from .radii import types_radii
from MDAnalysis import Universe
import MDAnalysis as mda
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

class PackingDefect2:
    def __init__(self):
        pass

    def read_top(self, resname, topology_file):
        """
        Examples
        --------
        pd = smda.PackingDefect2()
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
        TGglyc = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32',
                  'C1', 'C2', 'C3', 'C11', 'C21', 'C31',
                  'HA', 'HB', 'HS', 'HX', 'HY']

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

                ### PL headgroup: -1
                ### PL acyl chains: 1
                ### TG glycerol: 2
                ### TG tail cahins: 3
                if resname != 'TRIO':
                    if sl[1] in tails:
                        acyl = 1
                    else:
                        acyl = -1
                else:
                    if sl[1] in TGglyc:
                        acyl = 2
                    else:
                        acyl = 3

                output[sl[1]] = [types_radii[sl[2]], acyl]
        return output


    def apply_top(self, ag, output):
        u = ag.universe
        assert hasattr(ag, 'radii'), 'radii not defined!'
        assert hasattr(ag, 'bfactors'), 'bfactors not defined!'
        first_residue = ag.residues[0]
        n_residues = ag.residues.n_residues

        radii = []; types = []
        for atom in first_residue.atoms:
            r, t = output[atom.name]
            radii.append(r)
            if t == 'a':
                types.append(0)
            else:
                types.append(1)

        ag.radii = radii * n_residues
        ag.bfactors = types * n_residues


    def defect_size(self, matrices, nbins=400, bin_max=150,
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


class PackingDefect2PMDA(ParallelAnalysisBase):
    def __init__(self, atomgroups, radii, nbins=400, bin_max=150):
        u = atomgroups[0].universe
        self.N  = 3000 #The maximum number of defects
        self.dt = u.trajectory[0].dt
        self.dx = 1
        self.dy = 1
        self.radii = radii
        self.nbins = nbins
        self.bin_max = bin_max
        super(PackingDefect2PMDA, self).__init__(u, atomgroups)

    def _prepare(self):
        pass

    def _single_frame(self, ts, atomgroups):
        ag = atomgroups[0]
        dim = ts.dimensions.copy()
        pbc = dim[0:3]
        print('time: {:.3f}    pbc: {:.3f} {:.3f} {:.3f}'.format(ts.time/1000, pbc[0], pbc[1], pbc[2]))
        pbc_xy0 = np.array([pbc[0],pbc[1],0])
        pbc_xyz = np.array([pbc[0],pbc[1],pbc[2]])
        aa = ag.universe.atoms
        aa.positions -= pbc_xy0 * np.floor(aa.positions / pbc_xyz)
        hz = np.average(ag.select_atoms('name P').positions[:,2])

        xarray = np.arange(0, pbc[0], self.dx)
        yarray = np.arange(0, pbc[1], self.dy)
        xx, yy = np.meshgrid(xarray, yarray)

        M = {}; M['up'] = np.zeros_like(xx); M['dw'] = np.zeros_like(xx)
        Z = {}; Z['up'] = np.zeros_like(xx); Z['dw'] = np.zeros_like(xx)
        Z['up'] += hz; Z['dw'] += hz

        zlim = {}
        zlim['up'] = np.max(ag.positions[:,2])
        zlim['dw'] = np.min(ag.positions[:,2])

        PL = {}
        PL['up'] = ag.select_atoms('name P and prop z > %f' %hz).center_of_mass()[2]
        PL['dw'] = ag.select_atoms('name P and prop z < %f' %hz).center_of_mass()[2]

        atoms = {}
        atoms['up'] = ag.select_atoms('prop z > %f' %(PL['up'] - 20))
        atoms['dw'] = ag.select_atoms('prop z < %f' %(PL['dw'] + 20))

        for l in ['up', 'dw']:
            for atom in atoms[l]:
                xatom, yatom, zatom = atom.position

                if l == 'up': assert zatom > PL[l] - 20, 'check Z pos'
                if l == 'dw': assert zatom < PL[l] + 20, 'check Z pos'

                radius, acyl = self.radii[atom.resname][atom.name]

                dxx =  xx - xatom
                dxx -= pbc[0] * np.around(dxx/pbc[0])

                dyy =  yy - yatom
                dyy -= pbc[1] * np.around(dyy/pbc[1])

                dist_meet = (np.sqrt(self.dx**2 + self.dy**2)/2 + radius)**2
                bAr = dxx ** 2 + dyy ** 2 < dist_meet

                if acyl == -1:
                    M[l][bAr] = acyl
                    continue

                else:
                    bAnP = M[l] >= 0

                    if l == 'up':
                        baZ = zatom > Z[l]
                    else:
                        baZ = zatom < Z[l]

                    bA = bAr & bAnP & baZ
                    M[l][bA] = acyl
                    Z[l][bA] = zatom

        return M['up'], M['dw'], PL['up']+5, PL['dw']-5, dim
        #return M['up'], M['dw'], zlim['up'], zlim['dw'], dim


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
                            residue_segindex = [0] * N,
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

        defects = ['Deep', 'PLacyl', 'TGglyc', 'TGacyl']

        defect_uni = {}; defect_clu = {}
        for d in defects: defect_uni[d] = df.copy()
        for d in defects: defect_clu[d] = []
        defect_thr = {'Deep': 0,  'PLacyl': 1,  'TGglyc': 2,  'TGacyl': 3}

        for d in defects:
            for i, ts in enumerate(defect_uni[d].trajectory):
                num = 0

                bA  = (Mup[i] == defect_thr[d])
                defect_clu[d].append(bA.astype(int))
                ind = np.where(bA)
                xs  = ind[1]; ys = ind[0]

                for x1, y1 in zip(xs, ys):
                    pos = np.array([x1, y1, zlimup[i]])
                    defect_uni[d].atoms[num].position = pos
                    num += 1

                bA = (Mdw[i] == defect_thr[d])
                defect_clu[d].append(bA.astype(int))
                ind = np.where(bA)
                xs = ind[1]; ys = ind[0]

                for x1, y1 in zip(xs, ys):
                    pos = np.array([x1, y1, zlimdw[i]])
                    defect_uni[d].atoms[num].position = pos
                    num += 1

        ### DEFECT LOCALIZATION
        for d in defects:
            u = defect_uni[d]
            u.trajectory[-1]
            u.atoms.write(d + '2.gro')
            with mda.Writer(d + '2.xtc', u.atoms.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(u.atoms)

        ### DEFECT CLUSTER
        PD = PackingDefect2()
        for d in defects:
            PD.defect_size(defect_clu[d], fname=d + '2.dat', nbins=self.nbins, bin_max=self.bin_max)

