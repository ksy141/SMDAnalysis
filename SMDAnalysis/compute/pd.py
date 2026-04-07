from pmda.parallel import ParallelAnalysisBase
import numpy as np
import warnings
warnings.filterwarnings("ignore")

from ._base import PackingDefectBase

_MODES = ('accumulative', 'topmost')


class PackingDefect(PackingDefectBase):
    """Packing defect topology helper.

    Parameters
    ----------
    mode : {'accumulative', 'topmost'}
        ``'accumulative'`` — float-weighted encoding; measures packing density.
        ``'topmost'``      — integer categorical encoding; identifies surface lipid type.

    Examples
    --------
    import os
    pd = PackingDefect(mode='accumulative')
    ff = os.environ['CHARMM_TOPPAR']
    radii = {
        'POPC': pd.read_top('POPC', ff + '/top_all36_lipid.rtf'),
        'TRIO': pd.read_top('TRIO', ff + '/trio.str'),
    }
    ana = PackingDefectPMDA(atomgroups, radii, mode='accumulative')
    """

    def __init__(self, mode='accumulative'):
        if mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}, got {mode!r}")
        self.mode = mode

    def read_top(self, resname, topology_file):
        """Return {atom_name: [radius, label]} for *resname*.

        Encoding depends on mode:

        accumulative — PL headgroup=1e6, PL acyl=1e3, TG glycerol=1, TG acyl=1e-3
        topmost      — PL headgroup=-1,  PL acyl=1,   TG glycerol=2, TG acyl=3
        """
        if self.mode == 'accumulative':
            return self._parse_topology(
                resname, topology_file,
                pl_acyl=1e3, pl_head=1e6,
                tg_glyc=1,   tg_acyl=1e-3,
            )
        else:
            return self._parse_topology(
                resname, topology_file,
                pl_acyl=1,  pl_head=-1,
                tg_glyc=2,  tg_acyl=3,
            )


class PackingDefectPMDA(PackingDefectBase, ParallelAnalysisBase):
    """Parallel packing defect analysis.

    Parameters
    ----------
    atomgroups : list[AtomGroup]
    radii : dict
        Built by ``PackingDefect(mode=mode).read_top()``.
    mode : {'accumulative', 'topmost'}
        Must match the mode used to build *radii*.
    nbins : int
        Number of histogram bins for cluster-size distribution.
    bin_max : float
        Upper edge of histogram (Å²).
    prefix : str
        Output file prefix / directory.
    prob : bool
        Normalise cluster-size histogram to probability.
    """

    def __init__(self, atomgroups, radii, mode='accumulative',
                 nbins=151, bin_max=150, prefix='./', prob=True):
        if mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}, got {mode!r}")
        u = atomgroups[0].universe
        self.N       = 0  # set dynamically after grid size is known
        self.dt      = u.trajectory[0].dt
        self.dx      = 1
        self.dy      = 1
        self.dz      = 1  # used only in accumulative mode
        self.mode    = mode
        self.radii   = radii
        self.nbins   = nbins
        self.bin_max = bin_max
        self.prefix  = prefix
        self.prob    = prob
        # Scan trajectory for max box dimensions so all frames return same-shape
        # M arrays (required by PMDA) while still using per-frame pbc for wrapping.
        print("Scanning trajectory for maximum box dimensions...")
        max_x, max_y = 0.0, 0.0
        for ts in u.trajectory:
            max_x = max(max_x, ts.dimensions[0])
            max_y = max(max_y, ts.dimensions[1])
        u.trajectory[0]  # rewind
        self.nx = int(np.ceil(max_x / self.dx))
        self.ny = int(np.ceil(max_y / self.dy))
        self.N  = self.nx * self.ny * 2  # upper + lower leaflet
        super(PackingDefectPMDA, self).__init__(u, atomgroups)

    def run(self, start=None, stop=None, step=None,
            n_jobs=1, n_blocks=None):
        """Override PMDA run to avoid np.asarray on unequal-sized blocks.

        numpy ≥ 1.24 raises ValueError when blocks contain different numbers
        of frames.  We store ``_results`` as a plain list instead.
        """
        import time as _time
        from pmda.parallel import make_balanced_slices
        from dask import delayed
        from pmda.util import timeit

        if n_blocks is None:
            n_blocks = n_jobs

        start, stop, step = self._trajectory.check_slice_indices(
            start, stop, step)
        n_frames = len(range(start, stop, step))
        self.n_frames = n_frames

        slices = make_balanced_slices(n_frames, n_blocks,
                                      start=start, stop=stop, step=step)

        with timeit() as total:
            with timeit() as prepare:
                self._prepare()
            time_prepare = prepare.elapsed

            blocks = []
            _blocks = []
            with self.readonly_attributes():
                for bslice in slices:
                    task = delayed(
                        self._dask_helper, pure=False)(
                            bslice, self._indices,
                            self._top, self._traj)
                    blocks.append(task)
                    _blocks.append(range(bslice.start,
                                         bslice.stop, bslice.step))
                blocks = delayed(blocks)

                scheduler_kwargs = {'scheduler': 'multiprocessing',
                                    'num_workers': n_jobs}
                wait_start = _time.time()
                res = blocks.compute(**scheduler_kwargs)

            if len(res) == 0:
                res = [([], [], [], 0, wait_start, 0, 0)]

            with timeit() as conclude:
                # Store as list of arrays — avoids np.asarray on
                # blocks with different frame counts.
                self._results = [el[0] for el in res]
                self._blocks = _blocks
                self._conclude()

        from pmda.parallel import Timing
        self.timing = Timing(
            np.hstack([el[1] for el in res]),
            np.hstack([el[2] for el in res]), total.elapsed,
            np.array([el[3] for el in res]), time_prepare,
            conclude.elapsed,
            np.array([el[4] - wait_start for el in res]),
            np.array([el[5] for el in res]),
            np.array([el[6] for el in res]))
        return self

    def _prepare(self):
        pass

    # ------------------------------------------------------------------ #
    # Frame-level logic                                                    #
    # ------------------------------------------------------------------ #

    def _pack_result(self, Mup, Mdw, zlim_up, zlim_dw, dim):
        """Pack frame results into a single fixed-shape 2D array for PMDA."""
        result = np.zeros((2 * self.ny + 1, self.nx))
        result[0:self.ny, :]          = Mup
        result[self.ny:2*self.ny, :]  = Mdw
        result[2*self.ny, 0]          = zlim_up
        result[2*self.ny, 1]          = zlim_dw
        result[2*self.ny, 2:8]        = dim[0:6]
        return result

    def _unpack_result(self, result):
        """Unpack a single frame result from _pack_result."""
        Mup     = result[0:self.ny, :]
        Mdw     = result[self.ny:2*self.ny, :]
        zlim_up = result[2*self.ny, 0]
        zlim_dw = result[2*self.ny, 1]
        dim     = result[2*self.ny, 2:8]
        return Mup, Mdw, zlim_up, zlim_dw, dim

    def _single_frame(self, ts, atomgroups):
        if self.mode == 'accumulative':
            return self._frame_accumulative(ts, atomgroups)
        else:
            return self._frame_topmost(ts, atomgroups)

    def _frame_accumulative(self, ts, atomgroups):
        """3-D accumulative density scan (original pd.py)."""
        ag = atomgroups[0]
        dim = ts.dimensions.copy()
        pbc = dim[0:3]
        print('time: {:.3f}    pbc: {:.3f} {:.3f} {:.3f}'.format(
            ts.time / 1000, pbc[0], pbc[1], pbc[2]))

        pbc_xy0 = np.array([pbc[0], pbc[1], 0])
        pbc_xyz = np.array([pbc[0], pbc[1], pbc[2]])
        ag.universe.atoms.positions -= pbc_xy0 * np.floor(
            ag.universe.atoms.positions / pbc_xyz)

        hz = np.average(ag.select_atoms('name P').positions[:, 2])

        xarray = np.arange(self.nx) * self.dx
        yarray = np.arange(self.ny) * self.dy
        xx, yy = np.meshgrid(xarray, yarray)

        M    = {'up': np.zeros((self.ny, self.nx)), 'dw': np.zeros((self.ny, self.nx))}
        zlim = {'up': np.max(ag.positions[:, 2]),
                'dw': np.min(ag.positions[:, 2])}
        PL   = {
            'up': ag.select_atoms('name P and prop z > %f' % hz).center_of_mass()[2],
            'dw': ag.select_atoms('name P and prop z < %f' % hz).center_of_mass()[2],
        }

        C2 = ' '.join(['C2%d' % i for i in range(2, 23)])
        C3 = ' '.join(['C3%d' % i for i in range(2, 23)])
        memb = {
            'up': ag.select_atoms('(byres (name P and prop z > %f)) and name ' % hz + C2 + ' ' + C3),
            'dw': ag.select_atoms('(byres (name P and prop z < %f)) and name ' % hz + C2 + ' ' + C3),
        }
        utz = np.average(memb['up'].positions[:, 2])
        ltz = np.average(memb['dw'].positions[:, 2])

        atoms = {
            'up': ag.select_atoms('prop z > %f' % (utz - 3)),
            'dw': ag.select_atoms('prop z < %f' % (ltz + 3)),
        }

        # Cache C2 z-positions per residue to avoid repeated select_atoms calls
        c2_z = {res.resid: res.atoms.select_atoms('name C2')[0].position[2]
                for res in ag.residues
                if len(res.atoms.select_atoms('name C2')) > 0}

        dist_meet_3d = (np.sqrt(self.dx**2 + self.dy**2 + self.dz**2) / 2)

        for l in ['up', 'dw']:
            for atom in atoms[l]:
                xatom, yatom, zatom = atom.position
                radius, acyl = self.radii[atom.resname][atom.name]

                dxx = xx - xatom
                dxx -= pbc[0] * np.around(dxx / pbc[0])
                dyy = yy - yatom
                dyy -= pbc[1] * np.around(dyy / pbc[1])

                if acyl < 1e3 and l == 'up':    # TG
                    zarray = np.arange(utz, zlim['up'] + 1, self.dz)
                elif acyl < 1e3 and l == 'dw':  # TG
                    zarray = np.arange(ltz, zlim['dw'] - 1, -self.dz)
                elif acyl >= 1e3 and l == 'up': # PL
                    catom = c2_z[atom.resid]
                    zarray = np.arange(catom - 1, zlim['up'] + 1, self.dz)
                elif acyl >= 1e3 and l == 'dw': # PL
                    catom = c2_z[atom.resid]
                    zarray = np.arange(catom + 1, zlim['dw'] - 1, -self.dz)

                dzz = zarray - zatom
                dist_meet = (dist_meet_3d + radius)**2
                dist_lim  = (1 + radius)**2

                xind, yind = np.where((dxx**2 < dist_lim) & (dyy**2 < dist_lim))
                zind = np.where(dzz**2 < dist_lim)[0]

                if len(xind) == 0 or len(zind) == 0:
                    continue

                # Vectorized: broadcast (n_xy,) against (n_z,) → (n_xy, n_z)
                # Sum hits per (x,y) — matches original which added acyl per z-slice
                dx_vals = dxx[xind, yind]
                dy_vals = dyy[xind, yind]
                dist2 = (dx_vals[:, None]**2
                         + dy_vals[:, None]**2
                         + dzz[zind][None, :]**2)
                n_hits = (dist2 <= dist_meet).sum(axis=1)
                mask = n_hits > 0
                M[l][xind[mask], yind[mask]] += n_hits[mask] * acyl

        return self._pack_result(M['up'], M['dw'], PL['up'] + 5, PL['dw'] - 5, dim)

    def _frame_topmost(self, ts, atomgroups):
        """2-D topmost-atom scan (original pd2.py)."""
        ag = atomgroups[0]
        dim = ts.dimensions.copy()
        pbc = dim[0:3]
        print('time: {:.3f}    pbc: {:.3f} {:.3f} {:.3f}'.format(
            ts.time / 1000, pbc[0], pbc[1], pbc[2]))

        pbc_xy0 = np.array([pbc[0], pbc[1], 0])
        pbc_xyz = np.array([pbc[0], pbc[1], pbc[2]])
        ag.universe.atoms.positions -= pbc_xy0 * np.floor(
            ag.universe.atoms.positions / pbc_xyz)

        hz = np.average(ag.select_atoms('name P').positions[:, 2])

        xarray = np.arange(self.nx) * self.dx
        yarray = np.arange(self.ny) * self.dy
        xx, yy = np.meshgrid(xarray, yarray)

        M = {'up': np.zeros((self.ny, self.nx)), 'dw': np.zeros((self.ny, self.nx))}
        Z = {'up': np.zeros((self.ny, self.nx)) + hz, 'dw': np.zeros((self.ny, self.nx)) + hz}
        PL = {
            'up': ag.select_atoms('name P and prop z > %f' % hz).center_of_mass()[2],
            'dw': ag.select_atoms('name P and prop z < %f' % hz).center_of_mass()[2],
        }

        atoms = {
            'up': ag.select_atoms('prop z > %f' % (PL['up'] - 20)),
            'dw': ag.select_atoms('prop z < %f' % (PL['dw'] + 20)),
        }

        for l in ['up', 'dw']:
            for atom in atoms[l]:
                xatom, yatom, zatom = atom.position
                radius, acyl = self.radii[atom.resname][atom.name]

                dxx = xx - xatom
                dxx -= pbc[0] * np.around(dxx / pbc[0])
                dyy = yy - yatom
                dyy -= pbc[1] * np.around(dyy / pbc[1])

                dist_meet = (np.sqrt(self.dx**2 + self.dy**2) / 2 + radius)**2
                bAr = dxx**2 + dyy**2 < dist_meet

                if acyl == -1:  # headgroup: always overwrites
                    M[l][bAr] = acyl
                    continue

                bAnP = M[l] >= 0
                baZ  = zatom > Z[l] if l == 'up' else zatom < Z[l]
                bA   = bAr & bAnP & baZ
                M[l][bA] = acyl
                Z[l][bA] = zatom

        return self._pack_result(M['up'], M['dw'], PL['up'] + 5, PL['dw'] - 5, dim)

    # ------------------------------------------------------------------ #
    # Conclusion                                                           #
    # ------------------------------------------------------------------ #

    def _conclude(self):
        print("Concluding...")
        Mup = []; Mdw = []; zlimup = []; zlimdw = []; dim = []
        for r in self._results:
            for packed in r:
                mu, md, zu, zd, d = self._unpack_result(packed)
                Mup.append(mu);    Mdw.append(md)
                zlimup.append(zu); zlimdw.append(zd)
                dim.append(d)

        df = self._build_defect_universe(self.N, len(dim), self.dt, dim)

        defects = ['Deep', 'PLacyl', 'TGglyc', 'TGacyl']
        if self.mode == 'accumulative':
            defect_thr = {'Deep': [0, 1e-5], 'TGacyl': [1e-5, 1],
                          'TGglyc': [1, 1e3], 'PLacyl': [1e3, 1e6]}
        else:
            defect_thr = {'Deep': 0, 'PLacyl': 1, 'TGglyc': 2, 'TGacyl': 3}

        defect_uni = {d: df.copy() for d in defects}
        defect_clu = {d: []        for d in defects}

        for d in defects:
            thr = defect_thr[d]
            for i, _ in enumerate(defect_uni[d].trajectory):
                num = 0
                # Crop to actual per-frame box size before classification.
                # Phantom cells beyond the real box edge receive duplicated
                # PBC-wrapped contributions and must be excluded.
                actual_nx = int(np.round(dim[i][0] / self.dx))
                actual_ny = int(np.round(dim[i][1] / self.dy))
                for M_full, zlim in [(Mup[i], zlimup[i]), (Mdw[i], zlimdw[i])]:
                    M = M_full[:actual_ny, :actual_nx]
                    if self.mode == 'accumulative':
                        bA = (thr[0] <= M) & (M < thr[1])
                    else:
                        bA = (M == thr)
                    defect_clu[d].append(bA.astype(int))
                    ys, xs = np.where(bA)
                    for x1, y1 in zip(xs, ys):
                        defect_uni[d].atoms[num].position = [x1, y1, zlim]
                        num += 1

        self._write_defect_trajectories(defect_uni, defects, self.prefix)

        for d in defects:
            self.defect_size(defect_clu[d],
                             fname=self.prefix + d + '.dat',
                             nbins=self.nbins,
                             bin_max=self.bin_max,
                             prob=self.prob)
