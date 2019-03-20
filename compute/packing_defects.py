from __future__ import absolute_import, division, print_function
import time
import subprocess
import numpy as np
from ..common.distance import Distance
from MDAnalysis import Universe
from MDAnalysis import AtomGroup
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

        radii = {'C': 1.7, 'H': 1.2, 'O': 1.52, 'N': 1.55, 'P': 1.8}
        heads, glycerols, tails = self._selection()
        trios = list(set(self.u.select_atoms("resname TRIO").names))
        membranes = self.u.select_atoms("resname POPC DOPE SAPI")

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
            ns = FastNS(self.grid_size, membranes.positions, self.u.dimensions)

            zt = np.max(membranes.positions[:,2])
            zb = np.min(membranes.positions[:,2])

            dx = dy = dz = 1
            cr = 1.4142 / 2 * dx
            #cr = 1.7320 / 2 * dx
            edge = 0
            #cell_size = 5
            xis = int((pbc[0] - 2 * edge) / dx)
            yis = int((pbc[1] - 2 * edge) / dy)

            num = 0

            for xi in range(xis):
                for yi in range(yis):
                    x = edge + dx * xi
                    y = edge + dy * yi
                    z = zt
                    #if self.debug:
                    #    if abs(y-76.000) > 0.1: continue
                    #    if abs(x-40.000) > 0.1: continue
                    #print("NOX X, Y IS %.3f %.3f" % (x, y))
                    #xx = np.arange(x - cell_size, x + cell_size, dx)
                    #yy = np.arange(y - cell_size, y + cell_size, dy)

                    glyatom_z = 0
                    while (z - glyatom_z) > -0.5:
                        if self.debug:
                            print("z =", z)
                        # while z > pbc[2]/2:
                        r = np.array([x, y, z])
                        #print("Z is %.3f" % z)
                        # zz = np.zeros(len(xx)) + (z - 1)
                        #cells = np.array(np.meshgrid(xx, yy, zz)).T.reshape(-1, 3)
                        cells = np.array([r])

                        ### TIME CONSUMING PART ###
                        indices = ns.search(cells).get_indices()
                        #indices.extend(ns.self_search(cells).get_indices())
                        indices = list(set([item for sublist in indices for item in sublist]))

                        d = 10000
                        min_atname = None
                        min_index = None
                        min_dist = 10000
                        d_gly = 10000
                        min_glyindex = None
                        found_gly = False


                       # dist2 = Distance(r, membranes.positions, pbc).distance2(pbc = True)
                       # min_index = np.argmin(dist2)
                       # min_atom  = membranes[min_index]
                       # min_dist  = np.sqrt(np.dist2[min_index])


                        for ai in indices:
                            atom = membranes[ai]
                            atname = atom.name
                            d = Distance(r, atom.position, pbc).distance(pbc=True)
                            if d < min_dist:
                                min_atname = atname
                                min_dist = d
                                min_index = ai

                            if atname == 'C22' and d < d_gly:
                                min_glyindex = ai
                                d_gly = d
                                found_gly = True

                        if self.debug:
                            print("d_gly = ", d_gly)
                            try:
                                print("gly coordinate:", membranes[min_glyindex].position)
                            except: pass

                        # print(min_atname, min_dist)
                        if not found_gly:
                            sel = self.u.select_atoms("name C22")
                            dist2 = Distance(r, sel.positions, pbc).distance2(pbc = True)
                            iindex = np.argmin(dist2)
                            glyatom_z =  sel[iindex].position[2]
                            d_gly = np.sqrt(dist2[iindex])
                            if self.debug:
                                print("gly z (1) =", glyatom_z)
                        else:
                            glyatom = membranes[min_glyindex]
                            glyatom_z = glyatom.position[2]
                            if self.debug:
                                print("gly z (2) =", glyatom_z)


                        if min_atname == None:
                            z -= dz
                            continue

                        if min_dist < (cr + radii[min_atname[0]]):
                            if self.debug:
                                print("overlap happen ")
                            if min_atname in heads:
                                if self.debug:
                                    print("with HEAD")
                                pass
                            elif min_atname in glycerols:
                                if self.debug:
                                    print("with GLYCEROL")
                                pass
                            elif min_atname in tails:
                                if self.debug:
                                    print("with TAIL")
                                num += 1
                                file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                            elif min_atname in trios:
                                if self.debug:
                                    print("wtih TRIO")
                                num += 1
                                file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                                pass
                            #else:
                            #    print("SOMETHING WEIRD!")
                            break

                        #if min_glyindex == None:
                        #    z -= dz
                        #    continue

                        z -= dz

                        if (z - glyatom_z) <= -0.5:
                            # if z < pbc[2]/2 + 1:
                            if self.debug:
                                print("geometrical defect")
                            file.write('H %.3f %.3f %.3f\n' % (x, y, z))
                            num += 1

            file.close()
            subprocess.call(['sed', '-i.bak', 's/NUM/{:d}/'.format(num), 'top_mdanalysis.xyz'])
            subprocess.call(['rm', '-rf', 'top_mdanalysis.xyz.bak'])

        end = time.time()
        timelength = (end - start) / 60
        print("time spent: %.2f min" % timelength)

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

