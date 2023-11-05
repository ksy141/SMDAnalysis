import argparse
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from   scipy.cluster.hierarchy import single, fcluster, dendrogram, ward
from   MDAnalysis.analysis.distances import self_distance_array
from   collections import Counter
import warnings
warnings.filterwarnings("ignore")

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # top
    parser.add_argument("--top", required=True, type=str, nargs='*')

    # trj
    parser.add_argument("--trj", type=str, nargs='*', default=[])
    
    # b
    parser.add_argument("--b", type=int, default=None, help="u.trajectory[b:e:s]")

    # e
    parser.add_argument("--e", type=int, default=None, help="u.trajectory[b:e:s]")

    # s
    parser.add_argument("--s", type=int, default=None, help="u.trajectory[b:e:s]")

    # selection
    parser.add_argument("--sel", type=str, required=True)

    # number of clusters
    parser.add_argument("--ncluster", type=int, default=1, help="number of clusters")

    # cutoff
    parser.add_argument("--cutoff", type=float, default=20.0, help="cutoff distance (A)")

    # csv
    parser.add_argument("--csv", type=str, default=None)

    # png
    parser.add_argument("--png", type=str, default=None)

    # percentage
    parser.add_argument("--percentage",    action="store_true",  dest="percentage", help="use percentage")
    parser.add_argument("--no-percentage", action="store_false", dest="percentage", help="use number")
    parser.set_defaults(percentage=True)
    
    # inertia tensor
    parser.add_argument("--inertia",    action="store_true",  dest="inertia", help="calculate inertia tensor")
    parser.add_argument("--no-inertia", action="store_false", dest="inertia", help="do not calculate inertia tensor")
    parser.set_defaults(inertia=True)
    return parser.parse_args()


def main():
    args = parse_args()
    u    = mda.Universe(*args.top, *args.trj)
    print("TOP: ", *args.top)
    print("TRJ: ", *args.trj)
    ag   = u.select_atoms(args.sel)
    Ntot = ag.residues.n_residues
    
    frames = [] #frame number
    data   = {} #nuc
    iner   = {} #tensor
    for i in range(args.ncluster): data[i] = []
    for i in range(args.ncluster): iner[i] = []

    for i, ts in enumerate(u.trajectory[args.b : args.e : args.s]):
        frames.append(i)

        ### THIS SHOULD BE THE SAME WITH
        ### scipy.spatial.distance.pdist(ag.positions)
        ### if not taking into account periodic boundary conditions
        pdist = self_distance_array(ag.positions, box=u.dimensions)
        Z = single(pdist)
        
        f = fcluster(Z, t=args.cutoff, criterion='distance')
        result = Counter(f).most_common(args.ncluster)

        for i in range(args.ncluster):
            if args.percentage:
                data[i].append(result[i][1] / Ntot * 100)
            else:
                data[i].append(result[i][1])

        if args.inertia:
            for i in range(args.ncluster):
                idx = result[i][0]
                bA  = f == idx
                newag = ag[bA]
                newag.positions -= newag.center_of_geometry()

                Ixx = np.sum(newag.positions[:,1]**2 + newag.positions[:,2]**2)
                Iyy = np.sum(newag.positions[:,0]**2 + newag.positions[:,2]**2)
                Izz = np.sum(newag.positions[:,0]**2 + newag.positions[:,1]**2)

                Ixy = -np.sum(newag.positions[:,0] * newag.positions[:,1])
                Ixz = -np.sum(newag.positions[:,0] * newag.positions[:,2])
                Iyz = -np.sum(newag.positions[:,1] * newag.positions[:,2])

                T = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
                w, v = np.linalg.eig(T)
                k = 1.5 * (w[0]**4 + w[1]**4 + w[2]**4) / (w[0]**2 + w[1]**2 + w[2]**2)**2 - 0.5
                iner[i].append(k)

    final = [frames]
    for i in range(args.ncluster):
        final.append(data[i])
    
    if args.inertia:
        for i in range(args.ncluster):
            final.append(iner[i])
    
    final = np.array(final).T
    
    if args.csv:
        fmt = "%12d," + "%12.6f," * args.ncluster
        if args.inertia:
            fmt += "%12.6f," * args.ncluster
        np.savetxt(args.csv, final, fmt=fmt[:-1])

    if args.png:
        fig, ax = plt.subplots()
        for i in range(1, args.ncluster+1):
            ax.plot(final[:,0], final[:,i])

        ax.set_xlabel('frame')
        if args.percentage:
            ax.set_ylabel('nuc%')
        else:
            ax.set_ylabel('N')

        fig.tight_layout()
        fig.savefig(args.png)

if __name__ == '__main__':
    main()


