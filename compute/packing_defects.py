from __future__ import absolute_import, division, print_function
import time
import subprocess
import numpy as np
from ..common.distance import Distance
from MDAnalysis import Universe
from MDAnalysis.lib.nsgrid import FastNS


class PackingDefects:
    def __init__(self, u, debug = False, grid_size = 3):
        assert isinstance(u, Universe)
        self.u = u
        self.debug = debug
        self.grid_size = grid_size
        pass


    def compute_packing_defects(self, b=0, e=10000):
        ### b and e in ns ###
        start = time.time()

        file = open('top_mdanalysis.xyz', 'w')
        file.close()
        self.matrix = []

        radii = {'C': 1.7, 'H': 1.2, 'O': 1.52, 'N': 1.55, 'P': 1.8}
        heads, glycerols, tails = self._selection()
        trios = list(set(self.u.select_atoms("resname TRIO").names))

        membranes = self.u.select_atoms("resname POPC DOPE SAPI TRIO")
        c22s = self.u.select_atoms("name C22 and resname POPC DOPE SAPI")


        for ts in self.u.trajectory:

            t = self.u.trajectory.time / 1000 # t in ns
            if t > e:
                break
            elif b <= t <= e:
                pass
            elif t < b:
                continue

            print("Starts: %d frame" %ts.frame)
            file = open('top_mdanalysis.xyz', 'a')
            file.write('NUM\n')
            file.write('frame: %d       time: %.3f ns\n' %(ts.frame, t))
            pbc = self.u.dimensions[0:3]

            ns  = FastNS(self.grid_size, membranes.positions, self.u.dimensions)

            zt = np.max(membranes.positions[:,2])
            zb = np.min(membranes.positions[:,2])

            dx = dy = dz = 1
            cr = 1.4142 / 2 * dx
            edge = 0
            xis = int((pbc[0] - 2 * edge) / dx)
            yis = int((pbc[1] - 2 * edge) / dy)
            matrix = np.zeros((xis, yis))

            num = 0

            for xi in range(xis):
                for yi in range(yis):
                    x = edge + dx * xi
                    y = edge + dy * yi
                    z = zt


                    if self.debug:
                        if abs(x - 37.000) > 0.1: continue
                        if abs(y - 77.000) > 0.1: continue
                        print("NOX X, Y IS %.3f %.3f" % (x, y))


                    glyatom_z = 0
                    while (z - glyatom_z) > -1:
                        if self.debug:
                            print("Z is %.3f" % z)
                        r = np.array([x, y, z])
                        r_cell = np.zeros((1,3))
                        r_cell[0] = r - np.array([0, 0, 0])

                        indices = ns.search(r_cell).get_indices()
                        indices = list(set([item for sublist in indices for item in sublist]))

                        dist2 = np.array([Distance(r, membranes[i].position, pbc).distance2(pbc=True) for i in indices])

                        if len(dist2) == 0:
                            z -= dz
                            continue

                        dist_index = np.argmin(dist2)
                        memb_index = indices[dist_index]
                        min_atname = membranes[memb_index].name
                        min_resname = membranes[memb_index].resname
                        min_dist = np.sqrt(dist2[dist_index])

                        if self.debug:
                            new_dist2 = np.array(Distance(r, membranes.positions, pbc).distance2(pbc=True))
                            new_minindex = np.argmin(new_dist2)
                            new_atname = membranes[new_minindex].name
                            new_dist = np.sqrt(new_dist2[new_minindex])
                            print("coarsed: nearast atom: %6s  %6s  %6.3f " %(min_atname, memb_index, min_dist))
                            print("atomsel: nearest atom: %6s  %6s  %6.3f " %(new_atname, new_minindex, new_dist))


                        dist2 = np.array(Distance(r, c22s.positions, pbc).distance2(pbc=True))
                        c22s_cindex = np.argmin(dist2)
                        d_gly = np.sqrt(dist2[c22s_cindex])
                        glyatom_z = c22s[c22s_cindex].position[2]


                        if self.debug:
                            new_dist2 = np.array(Distance(r, c22s.positions, pbc).distance2(pbc=True))
                            new_minindex = np.argmin(new_dist2)
                            new_dist = np.sqrt(new_dist2[new_minindex])
                            print("coarsed: nearast glycerol atom: %5d  %6.3f " % (c22s_cindex, d_gly))
                            print("atomsel: nearest glycerol atom: %5d  %6.3f " % (new_minindex, new_dist))



                        if min_dist < (cr + radii[min_atname[0]]):
                            if self.debug:
                                print("overlap happen ")
                            if min_atname in heads and min_resname != 'TRIO':
                                if self.debug:
                                    print("with HEAD")
                                pass
                            elif min_atname in glycerols and min_resname != 'TRIO':
                                if self.debug:
                                    print("with GLYCEROL")
                                pass
                            elif min_atname in tails and min_resname != 'TRIO':
                                if self.debug:
                                    print("with TAIL")
                                num += 1
                                file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                                matrix[xi][yi] = 1

                            elif min_atname in trios and min_resname == 'TRIO':
                                if self.debug:
                                    print("wtih TRIO")
                                num += 1
                                file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                                matrix[xi][yi] = 1

                                pass

                            break


                        z -= dz

                        if (z - glyatom_z) <= -1:
                            if self.debug:
                                print("geometrical defect")
                                matrix[xi][yi] = 1
                            file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                            num += 1


            self.matrix.append(matrix)
            file.close()
            subprocess.call(['sed', '-i.bak', 's/NUM/{:d}/'.format(num), 'top_mdanalysis.xyz'])
            subprocess.call(['rm', '-rf', 'top_mdanalysis.xyz.bak'])

        end = time.time()
        timelength = (end - start) / 60
        print("time spent: %.2f min" % timelength)



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



    def defect_size(self, nblocks = 10, nbins=20):
        l = [k for i in self.matrix for j in i for k in j]
        amin = min(l)
        amax = max(l)
        del l
        print(amin, amax)
        bins = np.linspace(amin, amax, nbins)
        num_matrix = len(self.matrix)
        hist = np.zeros((nblocks, nbins-1))

        for i in range(nblocks):
            start = int(num_matrix/nblocks) * i
            end   = int(num_matrix/nblocks) * (i+1)

            defects = []
            for matrix in self.matrix[start:end]:
                graph = self._make_graph(matrix)

                visited = set([])
                for n in graph:
                    if n not in visited:
                        defect_loc = self._dfs(graph, n)
                        visited = visited.union(defect_loc)
                        defects.append(len(defect_loc))

            hist[i], bin_edges = np.histogram(defects, bins=bins, density=True)

        average = np.average(hist, axis=0)
        std = np.std(hist, axis=0)

        file = open('defect_histogram.dat', 'w')
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

