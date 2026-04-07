import numpy as np
from MDAnalysis import Universe
import MDAnalysis as mda
from .radii import types_radii


class PackingDefectBase:
    """Shared utilities for packing defect analysis.

    Subclasses supply the atom-type encoding scheme via read_top(), and the
    frame-level detection logic via _single_frame().  Everything else — DFS
    clustering, size-distribution histogramming, and defect-universe I/O — is
    shared here.
    """

    # ------------------------------------------------------------------ #
    # Topology parsing                                                     #
    # ------------------------------------------------------------------ #

    def _parse_topology(self, resname, topology_file,
                        pl_acyl, pl_head, tg_glyc, tg_acyl):
        """Parse a CHARMM RTF/STR file and return {atom_name: [radius, label]}.

        Parameters
        ----------
        resname : str
        topology_file : str
        pl_acyl, pl_head, tg_glyc, tg_acyl : float | int
            Encoding values used by the calling subclass.
        """
        tails = []
        tails += ['C2%d' % i for i in range(2, 23)]
        tails += ['C3%d' % i for i in range(2, 23)]
        tails += ['H%dR' % i for i in range(2, 23)]
        tails += ['H%dS' % i for i in range(2, 23)]
        tails += ['H%dX' % i for i in range(2, 23)]
        tails += ['H%dY' % i for i in range(2, 23)]
        tails += ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']
        TGglyc = ['O11', 'O21', 'O31', 'O12', 'O22', 'O32',
                  'C1', 'C2', 'C3', 'C11', 'C21', 'C31',
                  'HA', 'HB', 'HS', 'HX', 'HY']

        output = {}
        startread = False
        with open(topology_file) as top:
            for line in top:
                if line.startswith('!'):
                    continue
                if line.startswith('RESI {}'.format(resname)):
                    startread = True
                if startread and line.startswith('BOND'):
                    startread = False
                if startread and line.startswith('ATOM'):
                    sl = line.split()
                    if resname in  ['TRIO', 'TRIV', 'OOOTG']:
                        label = tg_glyc if sl[1] in TGglyc else tg_acyl
                    else:
                        label = pl_acyl if sl[1] in tails else pl_head
                    output[sl[1]] = [types_radii[sl[2]], label]
        return output

    # ------------------------------------------------------------------ #
    # DFS clustering (with periodic boundary conditions)                  #
    # ------------------------------------------------------------------ #

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
                    x = divmod(xi + dx, xis)[1]
                    y = divmod(yi + dy, yis)[1]
                    if matrix[x, y] == 1:
                        nlist.append(x * yis + y)
            graph[n] = set(nlist) - {n}
        return graph

    def defect_size(self, matrices, nbins, bin_max, fname, prob=True):
        bins = np.linspace(0, bin_max, nbins)
        defects = []
        for matrix in matrices:
            graph = self._make_graph(matrix)
            visited = set()
            for n in graph:
                if n not in visited:
                    defect_loc = self._dfs(graph, n)
                    visited |= defect_loc
                    defects.append(len(defect_loc))

        hist, edges = np.histogram(defects, bins)
        hist = hist.astype(np.float64)
        if np.sum(hist) == 0:
            return
        if prob:
            hist /= np.sum(hist)
        binp = 0.5 * (edges[1:] + edges[:-1])
        np.savetxt(fname, np.transpose([binp, hist]), fmt="%8.5f")

    # ------------------------------------------------------------------ #
    # Shared _conclude helpers                                             #
    # ------------------------------------------------------------------ #

    def _build_defect_universe(self, N, nframes, dt, dims):
        """Return an empty MDAnalysis Universe pre-loaded with *nframes* frames."""
        df = Universe.empty(
            n_atoms=N,
            n_residues=N,
            atom_resindex=np.arange(N),
            residue_segindex=[0] * N,
            trajectory=True,
        )
        df.add_TopologyAttr('resname', ['O'] * N)
        df.add_TopologyAttr('name',    ['O'] * N)
        df.add_TopologyAttr('resid', np.arange(N) + 1)

        fac = np.zeros((nframes, N, 3))
        df.load_new(fac, order='fac')
        df.trajectory[0].dt = dt
        for i, ts in enumerate(df.trajectory):
            df.trajectory[i].dimensions = dims[i]
        return df

    def _write_defect_trajectories(self, defect_uni, defects, prefix, suffix=''):
        """Write .gro and .xtc files for each defect type."""
        for d in defects:
            u = defect_uni[d]
            u.trajectory[-1]
            u.atoms.write('{}{}{}.gro'.format(prefix, d, suffix))
            with mda.Writer('{}{}{}.xtc'.format(prefix, d, suffix),
                            u.atoms.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(u.atoms)
