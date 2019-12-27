import numpy as np
from .radii import types_radii
from MDAnalysis import Universe

tails = []
tails += ['C2%d' %i for i in range(2, 23)]
tails += ['C3%d' %i for i in range(2, 23)]
tails += ['H%dR' %i for i in range(2, 23)]
tails += ['H%dS' %i for i in range(2, 23)]
tails += ['H%dX' %i for i in range(2, 23)]
tails += ['H%dY' %i for i in range(2, 23)]
tails += ['H16Z', 'H18T', 'H91', 'H101', 'H18Z', 'H20T']
#print(tails)

class PackingDefects:
    def __init__(self):
        pass

    def read_top(self, resname, topology_file):
        ''' {atomname1: [radius, acyl or not], ...}'''
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
    

    def array_to_gro(self, fname, atomname, x, y, z):
        ofile = open(fname, 'w')
        ofile.write('Title\n')
        ofile.write('%d\n' %len(x))
        i = 0
        for x1, y1 in zip(x, y):
            ofile.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(i + 1, 'TIP3', atomname, i + 1, x[i]/10, y[i]/10, z/10))
            i += 1
        ofile.write('   1.00000   1.00000   1.00000')
        ofile.close()

        a = Universe(fname)
        g = a.select_atoms('all')
        g.write(fname)

    def IsCircleIntersectsSquare(self, radius, dxx, dyy, dx, dy):
        half_dx = dx/2
        half_dy = dy/2
        SD = np.zeros(np.shape(dxx))
        
        t = dxx + half_dx
        bA = t < 0
        SD += bA.astype(int) * t**2
        
        bA2 = np.invert(bA)
        t = dxx - half_dx
        bA3 = t > 0
        bA4 = (bA2 & bA3)
        SD += bA4.astype(int) * t**2

        t = dyy + half_dy
        bA = t < 0
        SD += bA.astype(int) * t**2
        
        bA2 = np.invert(bA)
        t = dyy - half_dy
        bA3 = t > 0
        bA4 = (bA2 & bA3)
        SD += bA4.astype(int) * t**2

        return SD <= radius * radius
        

    def defect_size(self, matrices, nbins=100, bin_max=150, prob=True, fname='defect_histogram.dat'):
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
    
    
