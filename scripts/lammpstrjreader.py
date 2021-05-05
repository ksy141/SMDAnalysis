import numpy as np
import pandas as pd
import MDAnalysis as mda

class LAMMPSTRJReader:
    def __init__(self):
        """
        u = smda.LAMMPSTRJReader.read('step3.lammpstrj')
        """
        pass

    def read(self, filename):
        
        def logic(x, N, N2):
            if x % N < N2: return True
            return False

        nframes = 0
        pbc = []
        N = 0
        cols = None
        nheader = None


        f = open(filename)
        W = f.readlines()

        for i, line in enumerate(W):
            if line.startswith('ITEM: TIMESTEP'):
                nframes += 1
                continue
        
            if line.startswith('ITEM: NUMBER OF ATOMS'):
                N = int(W[i+1])
                continue
        
            if line.startswith('ITEM: BOX BOUNDS'):
                xsplit = W[i+1].split()
                ysplit = W[i+2].split()
                zsplit = W[i+3].split()
        
                pbcx   = float(xsplit[1]) - float(xsplit[0])
                pbcy   = float(ysplit[1]) - float(ysplit[0])
                pbcz   = float(zsplit[1]) - float(zsplit[0])
        
                pbc.append([pbcx, pbcy, pbcz, 90, 90, 90])
                continue
        
            if line.startswith('ITEM: ATOMS'):
                sl   = line.split()
                cols = sl[2:]
        
                if nheader is None:
                    nheader = i + 1
        
                continue
        
            sl = line.split()
            if len(line.split()) < 3: continue
        
        f.close()
        
        df = pd.read_csv(filename,
                skiprows = lambda x: logic(x, N + nheader, nheader),
                names  = cols,
                header = None,
                delim_whitespace = True)
        
        
        ### unit conversion
        # F: kcal/mol/A to kJ/mol/A
        # v: A/fs to A/ps
        
        bAfce = False; bAvel = False
        if 'fx' in cols: bAfce = True
        if 'vx' in cols: bAvel = True

        if bAfce: df[['fx', 'fy', 'fz']] *= 4.184
        if bAvel: df[['vx', 'vy', 'vz']] *= 1e3
        fac = df[['x', 'y', 'z']].to_numpy().reshape(nframes, N, 3)
        
        print('%6d atoms' %N)
        print('%6d timesteps' %nframes)
        print('vel: ', bAvel)
        print('fce: ', bAfce)

        new = mda.Universe.empty(N, n_residues = N,
                atom_resindex = np.arange(N),
                residue_segindex = [0] * N,
                trajectory = True,
                velocities = bAvel,
                forces     = bAfce)
        
        
        new.load_new(fac, forces = fac, order='fac')
        
        if bAvel: vel = df[['vx', 'vy', 'vz']].to_numpy().reshape(nframes, N, 3)
        if bAfce: fce = df[['fx', 'fy', 'fz']].to_numpy().reshape(nframes, N, 3)
        
        for i, ts in enumerate(new.trajectory):
            new.dimensions = pbc[i]
            if bAvel: new.atoms.velocities = vel[i]
            if bAfce: new.atoms.forces = fce[i]

        return new
        
        
        
        
