from __future__ import absolute_import, division, print_function
import time
import subprocess
import numpy as np
from ..common.distance import Distance
from MDAnalysis import Universe
from MDAnalysis.lib.nsgrid import FastNS
import multiprocessing
from multiprocessing import Pool


class PackingDefects:
    def __init__(self, u, grid_size = 2, debug = False, np = multiprocessing.cpu_count()):
        assert isinstance(u, Universe)
        self.u = u
        self.debug = debug
        self.np = np
        print("Number of processes: ", np)
        self.grid_size = grid_size
        self.radii = {'C': 1.7, 'H': 1.2, 'O': 1.52, 'N': 1.55, 'P': 1.8}
        self.heads, self.glycerols, self.tails = self._selection()
        self.trios = list(set(u.select_atoms("resname TRIO").names))
        self.membranes = u.select_atoms("resname POPC DOPE SAPI TRIO")
        self.c22s = u.select_atoms("name C22 and resname POPC DOPE SAPI")
        pass


    def test(self, num):
        print(num * num)


    def compute_packing_defects(self, b=0, e=10000, geo='top_g.xyz', che='top_c.xyz', all='top_a.xyz'):
        arr = np.arange(10)
        pool = Pool(processes=4)
        print(arr, pool)
        pool.map(unwrap_self, zip([self] * len(arr), arr))


        ### b and e in ns ###
        start = time.time()

        gfile = open(geo, 'w')
        cfile = open(che, 'w')
        afile = open(all, 'w')
        gfile.close()
        cfile.close()
        afile.close()

        geo_matrices = []
        che_matrices = []
        all_matrices = []

        dx = dy = self.dz = 1
        self.cr = 1.4142 / 2 * dx
        edge = 0


        for ts in self.u.trajectory:

            t = self.u.trajectory.time / 1000 # t in ns
            if t > e:
                break
            elif b <= t <= e:
                pass
            elif t < b:
                continue

            print("Starts: %d frame" %ts.frame)

            gfile = open(geo, 'a')
            cfile = open(che, 'a')
            afile = open(all, 'a')
            gfile.write('NUM\n')
            cfile.write('NUM\n')
            afile.write('NUM\n')
            gfile.write('frame: %d       time: %.3f ns\n' %(ts.frame, t))
            cfile.write('frame: %d       time: %.3f ns\n' % (ts.frame, t))
            afile.write('frame: %d       time: %.3f ns\n' % (ts.frame, t))

            self.pbc = self.u.dimensions[0:3]

            self.ns = FastNS(self.grid_size, self.membranes.positions, self.u.dimensions)

            self.zt = np.max(self.membranes.positions[:, 2])
            self.zb = np.min(self.membranes.positions[:, 2])

            xis = int((self.pbc[0] - 2 * edge) / dx)
            yis = int((self.pbc[1] - 2 * edge) / dy)
            geo_matrix = np.zeros((xis, yis))
            che_matrix = np.zeros((xis, yis))
            all_matrix = np.zeros((xis, yis))

            grid = [[x, y] for x in range(xis) for y in range(yis)]
            print(grid[0])


            with Pool(self.np) as p:
                result = p.map(self._defect_check, grid)


            c = g = a = 0
            for coord, check in result:
                x = coord[0]
                y = coord[1]
                z = coord[2]
                if check == 1:
                    che_matrix[x][y] = 1
                    all_matrix[x][y] = 1
                    c += 1
                    a += 1
                    cfile.write('H %.3f %.3f %.3f\n' % (x, y, z))
                    afile.write('H %.3f %.3f %.3f\n' % (x, y, z))

                elif check == 2:
                    geo_matrix[x][y] = 1
                    all_matrix[x][y] = 1
                    g += 1
                    a += 1
                    gfile.write('H %.3f %.3f %.3f\n' % (x, y, z))
                    afile.write('H %.3f %.3f %.3f\n' % (x, y, z))

            all_matrices.append(all_matrix)
            geo_matrices.append(geo_matrix)
            che_matrices.append(che_matrix)

            gfile.close()
            cfile.close()
            afile.close()


            subprocess.call(['sed', '-i.bak', 's/NUM/{:d}/'.format(a), all])
            subprocess.call(['sed', '-i.bak', 's/NUM/{:d}/'.format(c), che])
            subprocess.call(['sed', '-i.bak', 's/NUM/{:d}/'.format(g), geo])

            subprocess.call('rm -rf *.bak', shell=True)

        end = time.time()
        timelength = (end - start) / 60
        print("time spent: %.2f min" % timelength)



    def _defect_check(self, position):
        assert len(position) == 2, "Please provide [x, y] to _defect_check"
        glyatom_z = 0
        x = position[0]
        y = position[1]
        z = self.zt
        check = 0
        # check = 0: no defect
        # check = 1: chemical defect
        # check = 2: geometrical defect

        while (z - glyatom_z) > -1:
            if self.debug:
                print("Z is %.3f" % z)
            r = np.array([x, y, z])
            r_cell = np.zeros((1, 3))
            r_cell[0] = r - np.array([0, 0, 0])

            indices = self.ns.search(r_cell).get_indices()
            indices = list(set([item for sublist in indices for item in sublist]))

            dist2 = np.array([Distance(r, self.membranes[i].position, self.pbc).distance2(pbc=True) for i in indices])

            if len(dist2) == 0:
                z -= self.dz
                continue

            dist_index = np.argmin(dist2)
            memb_index = indices[dist_index]
            min_atname = self.membranes[memb_index].name
            min_resname = self.membranes[memb_index].resname
            min_dist = np.sqrt(dist2[dist_index])

            dist2 = np.array(Distance(r, self.c22s.positions, self.pbc).distance2(pbc=True))
            c22s_cindex = np.argmin(dist2)
            d_gly = np.sqrt(dist2[c22s_cindex])
            glyatom_z = self.c22s[c22s_cindex].position[2]

            if min_dist < (self.cr + self.radii[min_atname[0]]):
                if self.debug:
                    print("overlap happen ")
                if min_atname in self.heads and min_resname != 'TRIO':
                    if self.debug:
                        print("with HEAD")
                    pass
                elif min_atname in self.glycerols and min_resname != 'TRIO':
                    if self.debug:
                        print("with GLYCEROL")
                    pass
                elif min_atname in self.tails and min_resname != 'TRIO':
                    if self.debug:
                        print("with TAIL")
                    check = 1
                elif min_atname in self.trios and min_resname == 'TRIO':
                    if self.debug:
                        print("wtih TRIO")
                    check = 1
                    pass
                break

            z -= self.dz

            if (z - glyatom_z) <= -1:
                if self.debug:
                    print("geometrical defect")
                check = 2

        return [[x, y, z], check]



    def _dfs(self, graph, start):
        visited, stack = set(), [start]
        while stack:
            vertex = stack.pop()
            if vertex not in visited:
                visited.add(vertex)
                stack.extend(graph[vertex] - visited)
        return visited


    def _make_graph(self, frame_matrix):
        graph = {}
        xis, yis = frame_matrix.shape

        for xi in range(xis):
            for yi in range(yis):
                if frame_matrix[xi][yi] == 0:
                    continue

                n = xi * yis + yi
                nlist = []

                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if dx * dy == 0:
                            x = divmod(xi + dx, xis)[1]
                            y = divmod(yi + dy, yis)[1]
                            if frame_matrix[x, y] == 1:
                                ndn = x * yis + y
                                nlist.append(ndn)

                graph[n] = set(nlist) - set([n])
        return graph



    def defect_size(self, matrices, nblocks = 5, nbins=100, bin_max=100, density=False, out='defect_histogram.dat'):
        bins = np.linspace(0, 150, nbins)
        num_matrix = len(matrices)
        hist = np.zeros((nblocks, nbins-1))

        for i in range(nblocks):
            start = int(num_matrix/nblocks) * i
            end   = int(num_matrix/nblocks) * (i+1)

            defects = []
            for matrix in matrices[start:end]:
                graph = self._make_graph(matrix)

                visited = set([])
                for n in graph:
                    if n not in visited:
                        defect_loc = self._dfs(graph, n)
                        visited = visited.union(defect_loc)
                        defects.append(len(defect_loc))

            hist[i], bin_edges = np.histogram(defects, bins=bins, density=density)

        average = np.average(hist, axis=0)
        std = np.std(hist, axis=0)

        file = open(out, 'w')
        for i in range(nbins-1):
            file.write('{: .3f} {: .3f} {: .3f}\n'.format(bins[i], average[i], std[i]))
        file.close()



    def _selection(self):
        head = "name P N HP32 HP42 HP52 O14 "
        for i in range(1, 7):
            head += "C1%d H1%d H%d " % (i, i, i)
            if i > 1:
                head += "O%d HO%d P%d " % (i, i, i)
            if i < 6:
                head += "H1%dA H1%dB H1%dC " % (i, i, i)

        for i in range(1, 5):
            head += "HN%d H%d OC%d HO%d " % (i, i, i, i)
            if i > 1:
                head += "OP3%d OP4%d OP5%d " % (i, i, i)

        for i in range(1, 4):
            head += "O1%d " % i

        glycerol = "name C1 C2 C3 HA HB HC HR HS HT HX HY HZ O21 O22 O31 O32 C21 C31 "

        tail = "name H91 H101 "
        for i in range(2, 23):
            tail += "C2%d C3%d H%dR H%dS H%dT H%dX H%dY H%dZ " % (i, i, i, i, i, i, i, i)

        heads = head.split()[1:]
        glycerols = glycerol.split()[1:]
        tails = tail.split()[1:]

        return heads, glycerols, tails
