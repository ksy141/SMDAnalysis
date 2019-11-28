import numpy as np
import os
import glob
import re
import pandas as pd
import shutil
import subprocess

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

### HOW TO USE?
#import SMDAnalysis as smda
#a = smda.REUS_Analysis()
#a.get_data()
#a.get_prob(b=30)
#a.get_pmfs(b=30)


class REUS_Analysis:
    def __init__(self):
        pass

    def get_data(self, directory = "./", plumed = "plumed.*", colvar="colvar.", nopbc = True):
        data = []
        os.chdir(directory)
        plumed_files = glob.glob(plumed)

        for ifile in plumed_files:
            index = int(ifile.split('.')[1])

            with open(ifile) as f:
                for line in f:
                    if re.search("RESTRAINT", line):
                        kappa = self._get_info('KAPPA=', line)
                        at = self._get_info('AT=', line)
                        at = round(at, 3)

            ### pd.read_csv() is much faster than np.loadtxt()
            df = pd.read_csv(colvar + str(index), 
                             delim_whitespace = True, 
                             header = None, 
                             comment='#', 
                             usecols=[0,1])

            if nopbc:
                cols = df.values
                dy = cols[:,1] - at
                dy -= np.around(dy)
                cols[:,1] = at + dy

            else:
                cols = df.values

            data.append({'at': at, 'kappa': kappa, 'colvar': cols})

        self.data = sorted(data, key=lambda k: k['at'])
        self.at_min = self.data[0]['at']
        self.at_max = self.data[-1]['at']
        self.t_min  = self.data[0]['colvar'][:,0][0]/1000
        self.t_max  = self.data[0]['colvar'][:,0][-1]/1000

    
    def get_prob(self, b = 0, e = 100000, 
                 at_min = -0.6, at_max = 0.6, nbins = 100):

        assert b < self.t_max, "b > t_max, no data to analyze!"

        if e > self.t_max:
            e = self.t_max
        
        count_settings = {'bins': nbins, 'range': (at_min, at_max)}

        b *= 1000 # ns to ps
        e *= 1000 # ns to ps

        _, edges = np.histogram([-100000], **count_settings)
        bins = 0.5 * (edges[1:] + edges[:-1])
        print(bins)
        print(len(bins))
        fig, ax = plt.subplots()

        for w in self.data:
            cols = w['colvar']

            c = cols[(cols[:,0] >= b) & (cols[:,0] <= e)]
            count, _ = np.histogram(c[:,1], **count_settings, density=True)
            ax.plot(bins, count, label=w['at'])
        
        ax.legend()
        fig.tight_layout()
        fig.savefig('prob.pdf')


    def get_pmfs(self, b = 0, e = 100000, nblocks = 5,
                 T = 310, nbins = 100, tol = 1e-10,
                 idx_zero = 0):


        assert b < self.t_max, "b > t_max, no data to analyze!"

        if e > self.t_max:
            e = self.t_max

        b *= 1000
        e *= 1000
        ts = np.linspace(b, e, nblocks + 1)
        
        ### MAKE PMFS folder
        shutil.rmtree('pmfs', ignore_errors=True)
        os.makedirs('pmfs')
        

        ### OUTPUT col
        for w in self.data:
            cols = w['colvar']

            for i in range(nblocks):
                c = cols[(cols[:,0] >= ts[i]) & (cols[:,0] < ts[i+1])]
                np.savetxt('pmfs/col{}_{}'.format(w['at'], i), c)
        

        ### OUTPUT METADATA
        runfile = open('pmfs/run.sh', 'w')
        for i in range(nblocks):
            runfile.write('wham {} {} {} {} {} 0 metadata{}.dat pmf{}.dat\n'.format(\
                    self.at_min, self.at_max, nbins, tol, T, i, i))

            mfile = open('pmfs/metadata{}.dat'.format(i), 'w')
            for w in self.data:
                mfile.write('col{}_{} {} {}\n'.format(w['at'], i, w['at'], w['kappa']))
            mfile.close()
        runfile.close()
        
        ### RUN WHAM
        os.chdir('pmfs')
        subprocess.call(['bash', 'run.sh'])

        ### BLOCK AVERAGE
        pmfs = []
        for i in range(nblocks):
            p = np.loadtxt('pmf{}.dat'.format(i))
            x = p[:,0]
            y = p[:,1]
            y -= y[idx_zero]
            y *= 0.239 # kj to kcal
            pmfs.append(y)
        pmfs = np.array(pmfs)
        ave = np.average(pmfs, axis=0)
        std = np.std(pmfs, axis = 0)
        np.savetxt('pmf', np.transpose([x, ave, std, ave - std, ave + std]))
        
        ### PLOT
        fig, ax = plt.subplots()
        ax.plot(x, ave, color='C0')
        ax.fill_between(x, ave - std, ave + std, color='C0', alpha=0.4)
        fig.tight_layout()
        fig.savefig('pmf.pdf')
        
        os.chdir('../')



    def _get_info(self, string, line):
        _, __, interest = line.partition(string)
        return float(interest.split()[0])
