import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description='Plot Rosetta results')
parser.add_argument("-f",     help="filename e.g. score.fsc", default="score.fsc")
parser.add_argument("-pdbs",  help="path for pdbs", default="./")
parser.add_argument("-r",     help="rmsd w.r.t. r th best structure (1-based)", default = 1)
parser.add_argument("-add",   help="adding selection: protein and name CA and ???", default = '')
parser.add_argument("-num",   help="number of pdb files to be plotted. Default=1/10 of total pdb files", default=None)
parser.add_argument("-o",     help="output score-rmsd figure name", default="score-rmsd.pdf")

def sort_rmsd(f, r, pdbs, add, num, out):
    scores = []
    score_col = None
    descr_col = None
    if add:
        add = ' and ' + add

    for line in open(f):
        if line.startswith('SCORE'):
            sline = line.split()

            try:
                score = float(sline[score_col])
                descr = sline[descr_col]
                scores.append((score, descr))
            
            except:
                for i in range(len(sline)):
                    if sline[i] == 'total_score' or 'score':
                        score_col = i
                        print("score_col (0-based): ", i)
                    elif sline[i] == 'description':
                        descr_col = i
                        print("descr_col (0-based): ", i)

    scores.sort()

    rs = []
    ss = []

    if not pdbs.endswith('/'):
        pdbs += '/'

    d0 = scores[r-1][1]
    ref = mda.Universe(pdbs + d0 + '.pdb')
    check = 0

    if num == 'None':
        num = int(len(scores)/10)
    else:
        num = int(num)

    for score in scores[0:num]:
        s = score[0]
        d = score[1]
    
        mobile = mda.Universe(pdbs + d + '.pdb')
        old, new = align.alignto(mobile, ref, select="protein and name CA" + add)
    
        rs.append(new)
        ss.append(s)
    
        if check < 10:
            print(s, d, new)
        check += 1
    
    
    fig, ax = plt.subplots()
    ax.plot(rs, ss, linestyle='', marker='o')
    ax.set_xlabel('RMSD (A)')
    ax.set_ylabel('score')
    fig.savefig(out)
    


args = parser.parse_args()
f = args.f
r = int(args.r)
pdbs = args.pdbs
add  = args.add
num  = args.num
out  = args.o

sort_rmsd(f, r, pdbs, add, num, out)




