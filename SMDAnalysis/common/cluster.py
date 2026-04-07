import numpy as np
# from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

class Cluster:
    def __init__(self):
        pass

    def dd_to_iid(self, dd, N):
        ### dd
        # dd[i, j] = distance between index i and j
        # dd.shape = N * N

        ### N
        # N can be retrieved by the following method
        #N = (1 + np.sqrt(1 + 4 * len(dd)**2))/2
        #N = int(N)

        ### iid
        # [[index1, index2, distance], ...]
        # sorted by distance
        # iid.shape = N(N-1)/2 x 3
        
        iid = []
        for i in range(N):
            for j in range(i, N):
                if i == j:
                    continue
                iid.append([i, j, dd[i, j]])
    
        iid = sorted(iid, key=lambda x: x[2])
        return iid


    def iid_to_linkage(self, iid, N):
        indices = np.arange(N, dtype=np.int64)
        ### indices
        # Each element describes at which cluster it belongs to
        # Updated every loop
        # initially, [0, 1, 2, ..., N-1]
        # finally, [2N-2, 2N-2, 2N-2, ... ]
    
        Z = []
        ### Z 
        # [[index1, index2, distance, # of points], ...]
        # Z.shape = (N-1) x 4
        # All indices idx >= len(X) refer to the cluster formed in Z[idx - len(X)].
        # Z = linkage(X, 'simple')
    
        i = 0
        j = 0
        while len(Z) < N-1:
            new_ind = N + j
            ind1, ind2, d = iid[i]
            if indices[ind1] == indices[ind2]:
                i += 1
                continue
            tmp = sorted([indices[ind1], indices[ind2]]) + [d]
    
            indices[indices == indices[ind1]] = new_ind
            indices[indices == indices[ind2]] = new_ind
            Z.append(tmp + [np.sum(indices == new_ind)])
            i += 1
            j += 1
        Z = np.array(Z)
        return Z
    
