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
    parser.add_argument("--no-percentage", action="store_false", dest="percentage", help="use percentage")
    parser.set_defaults(percentage=True)

    return parser.parse_args()


def main():
    args = parse_args()
    u    = mda.Universe(*args.top, *args.trj)
    print("TOP: ", *args.top)
    print("TRJ: ", *args.trj)
    ag   = u.select_atoms(args.sel)
    Ntot = ag.residues.n_residues
    
    frames = []
    data   = {}
    for i in range(args.ncluster): data[i] = []

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
    
    final = [frames]
    for i in range(args.ncluster):
        final.append(data[i])
    final = np.array(final).T
    
    if args.csv:
        np.savetxt(args.csv, final)

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





