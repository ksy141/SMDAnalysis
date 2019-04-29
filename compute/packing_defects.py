from __future__ import absolute_import, division, print_function
import time
import subprocess
import numpy as np
from ..common.distance import Distance
from MDAnalysis import Universe
from MDAnalysis.lib.nsgrid import FastNS


class PackingDefects:
    def __init__(self, u, grid_size = 2, debug = False):
        assert isinstance(u, Universe)
        self.u = u
        self.debug = debug
        self.grid_size = grid_size
        self.radii = {'C': 1.7, 'H': 1.2, 'O': 1.52, 'N': 1.55, 'P': 1.8, 'S': 1.8}
        self.heads, self.glycerols, self.tails = self._selection()
        self.trios = list(set(u.select_atoms("resname TRIO").names))
        self.membranes = u.select_atoms("resname POPC DOPE SAPI TRIO")
        self.c22s = u.select_atoms("name C22 and resname POPC DOPE SAPI")
        self.proteins = u.select_atoms("protein")
        self.dx = self.dy = 1
        self.dz = 0.1
        self.cr = 1.4142 / 2 * self.dx
        self.edge = 0
        self.pbc = None
        self.cut = -1

        self.glyatom_z = None
        self.start = time.time()

    def save_to_xyz(self, matrices, file):
        ofile = open(file, 'w')
        iframe = 0
        for matrix in matrices:
            num = len(matrix[matrix!=0])
            ofile.write('%d\n' %num)
            ofile.write('frame: %d\n' %iframe)
            iframe += 1

            for (x, y), z in np.ndenumerate(matrix):
                if z != 0:
                    ofile.write('H %.3f %.3f %.3f\n' % (x, y, z))
        ofile.close()


    def make_grid(self):
        ttime = self.u.trajectory.time/1000
        if ttime % 10 == 0:
            end = time.time()
            dt = (end - self.start)/60
            print("simulation time: %.1f ns  |  time spent: %d min" %(ttime, dt))
        self.pbc = self.u.dimensions[0:3]
        self.ns = FastNS(self.grid_size, self.membranes.positions, self.u.dimensions)
        try:
            self.nsp= FastNS(self.grid_size, self.proteins.positions, self.u.dimensions)
        except:
            pass
        self.zt = np.max(self.membranes.positions[:, 2])
        self.zb = np.min(self.membranes.positions[:, 2])

        self.xis = int((self.pbc[0] - 2 * self.edge) / self.dx)
        self.yis = int((self.pbc[1] - 2 * self.edge) / self.dy)

        grid = [[x, y] for x in range(self.xis) for y in range(self.yis)]
        return grid


    def store_defect(self, result):
        xis = self.xis
        yis = self.yis
        che_matrix = np.zeros((xis, yis))
        geo_matrix = np.zeros((xis, yis))

        c = g = 0
        for coord, check in result:
            x = coord[0]
            y = coord[1]
            z = coord[2]
            if check == 1:
                che_matrix[x][y] = z
                c += 1

            if check == 2:
                geo_matrix[x][y] = z
                che_matrix[x][y] = z
                g += 1
                c += 1

        return che_matrix, geo_matrix


    def defect_check(self, position):
        # closest C22 atom
        assert len(position) == 2, "Please provide [x, y] to defect_check"
        glyatom_z = 0
        x = position[0]
        y = position[1]
        z = self.zt
        check = 0
        # check = 0: no defect
        # check = 1: chemical defect
        # check = 2: geometrical defect

        while (z - glyatom_z) > self.cut:
            # if self.debug:
            #     print("Z is %.3f" % z)
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
                # if self.debug:
                #     print("overlap happen ")
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

            if (z - glyatom_z) <= self.cut:
                if self.debug:
                    print("geometrical defect")
                check = 2

        return [[x, y, z], check]

    

    def defect_check_protein(self, position):
        # closest C22 atom
        assert len(position) == 2, "Please provide [x, y] to defect_check"
        glyatom_z = 0
        x = position[0]
        y = position[1]
        z = self.zt
        check = 0
        # check = 0: no defect
        # check = 1: chemical defect
        # check = 2: geometrical defect

        while (z - glyatom_z) > self.cut:
            # if self.debug:
            #     print("Z is %.3f" % z)
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
            

            ### ADD HERE ###
            protein_check = False
            indicesp = self.nsp.search(r_cell).get_indices()
            indicesp = list(set([item for sublist in indicesp for item in sublist]))
            dist2p = np.array([Distance(r, self.proteins[i].position, self.pbc).distance2(pbc=True) for i in indicesp])
            try:
                dist_indexp = np.argmin(dist2p)
                prot_index = indicesp[dist_indexp]
                min_atnamep = self.proteins[prot_index].name
                min_distp = np.sqrt(dist2p[dist_indexp])
    
                if min_distp < min_dist:
                    min_dist = min_distp
                    min_atname = min_atnamep
                    protein_check = True
            except:
                pass
            ### ADD ABOVE ###

            dist2 = np.array(Distance(r, self.c22s.positions, self.pbc).distance2(pbc=True))
            c22s_cindex = np.argmin(dist2)
            d_gly = np.sqrt(dist2[c22s_cindex])
            glyatom_z = self.c22s[c22s_cindex].position[2]

            if min_dist < (self.cr + self.radii[min_atname[0]]):
                # if self.debug:
                #     print("overlap happen ")
                
                if protein_check:
                    if self.debug:
                        print("with protein ")
                    pass

                elif min_atname in self.heads and min_resname != 'TRIO':
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

            if (z - glyatom_z) <= self.cut:
                if self.debug:
                    print("geometrical defect")
                check = 2

        return [[x, y, z], check]




    def defect_check2(self, position):
        # Average glycerol z
        assert len(position) == 2, "Please provide [x, y] to defect_check"
        x = position[0]
        y = position[1]
        z = self.zt
        check = 0
        # check = 0: no defect
        # check = 1: chemical defect
        # check = 2: geometrical defect

        while (z - self.glyatom_z) > -1:
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

            if min_dist < (self.cr + self.radii[min_atname[0]]):
                if min_atname in self.tails and min_resname != 'TRIO':
                    check = 1
                elif min_atname in self.trios and min_resname == 'TRIO':
                    check = 2
                break

            z -= self.dz

            if (z - self.glyatom_z) <= -1:
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
                    if dx * dy == 0:
                        x = divmod(xi + dx, xis)[1]
                        y = divmod(yi + dy, yis)[1]
                        if matrix[x, y] == 1:
                            ndn = x * yis + y
                            nlist.append(ndn)

            graph[n] = set(nlist) - set([n])
        return graph



    def defect_size(self, matrices, nblocks = 5, nbins=100, bin_max=2.5, prob=False, file='defect_histogram.dat'):
        nbins += 1
        bins = np.linspace(0, bin_max, nbins)
        dbin = bins[1] - bins[0]
        num_matrix = len(matrices)
        hist = np.zeros((nblocks, nbins-1))

        for matrix in matrices:
            matrix[matrix != 0] = 1

#        for i in range(nblocks):
#            start = int(num_matrix/nblocks) * i
#            end   = int(num_matrix/nblocks) * (i+1)
#            period = end - start
#
#            defects = []
#            for matrix in matrices[start:end]:
#                graph = self._make_graph(matrix)
#
#                visited = set([])
#                for n in graph:
#                    if n not in visited:
#                        defect_loc = self._dfs(graph, n)
#                        visited = visited.union(defect_loc)
#                        defects.append(len(defect_loc))
#
#            hist[i], bin_edges = np.histogram(defects, bins=bins, density=density)
#            hist[i] *= dbin
#       
#        hist /= period
#        average = np.average(hist, axis=0)
#        std = np.std(hist, axis=0)
#
        
        hist = []
        for matrix in matrices:
            defects = []
            graph = self._make_graph(matrix)
            visited = set([])
            for n in graph:
                if n not in visited:
                    defect_loc = self._dfs(graph, n)
                    visited = visited.union(defect_loc)
                    defects.append(len(defect_loc) * 0.01) #A2 to nm2
            
            tmphist, bin_edges = np.histogram(defects, bins=bins, density=prob)
            if prob:
                tmphist *= dbin
            hist.append(tmphist)
        
        hist    = np.array(hist)
        average = np.average(hist, axis=0)
        std     = np.std(hist, axis=0)

        ofile = open(file, 'w')
        for i in range(nbins-1):
            ofile.write('{: .3f} {: .8f} {: .8f}\n'.format(bins[i], average[i], std[i]))
        ofile.close()



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
