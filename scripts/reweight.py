import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-f",     help="filename e.g. colvar", default="colvar")
parser.add_argument("-b",     help="time to read (ns)", default=0)
#parser.add_argument("-t",  help="time column", default=1)
parser.add_argument("-r",     help="rbias column", default=False)
parser.add_argument("-cv",    help="cv column", default=False)
parser.add_argument("-T",     help="Temperature (K)", default=310)
parser.add_argument("-min",   help="min in cv", default=-0.5)
parser.add_argument("-max",   help="max in cv", default=+0.5)
parser.add_argument("-nbins", help="nbins in cv", default=100)
parser.add_argument("-sigma", help="sigma in gauss", default=0.1)
parser.add_argument("-mintozero", help="make min(pmf) = 0", default=False)

args = parser.parse_args()

def reweight(f = args.f, b = float(args.b), r = int(args.r), cv = int(args.cv), T=float(args.T),
             minbin = float(args.min), maxbin = float(args.max), nbins = int(args.nbins), 
             sigma = float(args.sigma), mintozero = (str(args.mintozero) == 'True')):

    ## column number from 1-based to 0-based
    print("Temperature(K): ", T)
    beta = 1/0.008311/T
    r  -= 1
    cv -= 1
    
    #arr_time  = []
    arr_rbias = []
    arr_cv    = []

    for line in open(f):
        if line.startswith("#"):
            continue

        sline = line.split()
        if float(sline[0]) < b * 1000:
            continue
        
        #if float(sline[0]) == 0:
        #    print(sline)

        arr_rbias.append(float(sline[r]))
        arr_cv.append(float(sline[cv]))

    arr_rbias = np.array(arr_rbias)
    arr_cv    = np.array(arr_cv)
    
    #print(arr_cv[0:10])
    #print(arr_rbias[0:30])

    arr_erbias = np.exp(beta * arr_rbias)

    prob = np.zeros(nbins + 1)
    bins = np.linspace(minbin, maxbin, nbins+1)

    for cv, erbias in zip(arr_cv, arr_erbias):
        gauss = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( -(bins - cv)**2 / (2 * sigma**2))
        prob += gauss * erbias
    
    prob /= np.sum(prob)
    F = -1/beta * np.log(prob)

    if mintozero:
        F -= F[np.argmin(F)]

    return np.transpose([bins, prob, F])


x = reweight()
np.savetxt('reweight', x)

