from MDAnalysis.core.groups import AtomGroup
from MDAnalysis import Universe
import numpy as np

class LAMMPSTRJWriter:
    def __init__(self):
        """ Write lammpstrj file
        smda.LAMMPSTRJWriter().write(
            'CG.lammpstrj', u or u.select_atoms('name SOL SOD CLA'), 
            types = {} or {'SOL':1, 'SOD':2, 'CLA':3}, b=0, e=1e10, skip=1)
        """
        pass


    def write(self, filename, obj, types={}, pbc='pp pp pp',
              b=0, e=None, skip=1, kj2kcal=True):
        """
        Parameters
        ----------
        filename [str]
        obj = AtomGroup or Universe
        types = {'SOL': 1, 'SOD': 2, 'CLA': 3} or {}
        pbc = 'pp pp pp'
        b = 0    frame begins at b [int]
        e = None frame ends   at e [int]
        skip = 1 [int]
        kj2kcal = True [bool]
        
        When you save only some group of atoms, 
        provide AtomGroup of interest instead of Universe.

        When types is not given:
        It will assign the number based on atom name.
        Use for - if - else loop.
        When types is given for some atoms:
        Problematic!
        """

        if isinstance(obj, AtomGroup):
            u = obj.universe
            ag = obj

        elif isinstance(obj, Universe):
            u = obj
            ag = u.atoms

        else:
            assert 1==0, "input either obj = Universe or AtomGroup"

        try:
            u.atoms.forces
            force_b = True
            w = "{:6d} {:5s} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}\n" 
            w = "%6d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f"
            item = "ITEM: ATOMS id type x y z fx fy fz"
        except:
            force_b = False
            w = "{:6d} {:5s} {:8.3f} {:8.3f} {:8.3f}\n" 
            w = "%6d %5d %8.3f %8.3f %8.3f" 
            item = "ITEM: ATOMS id type x y z"

        ### name2type
        t = 1
        for atom in ag.atoms:
            if atom.name in types.keys():
                atom.type = types[atom.name]
            else:
                atom.type = str(t)
                types[atom.name] = str(t)
                t += 1
        print(types)
        

        ### Write
        f = open(filename, 'w')
        f.close()
        f = open(filename, 'ab')
        for i, ts in enumerate(u.trajectory[b:e:skip]):
            if i%100 == 0:
                print("%d/%d processing..." %(i, u.trajectory.n_frames))
            dim = u.dimensions
            info = """ITEM: TIMESTEP
{:d}
ITEM: NUMBER OF ATOMS
{:d}
ITEM: BOX BOUNDS {:s}
0.000000 {:11.6f}
0.000000 {:11.6f}
0.000000 {:11.6f}
""".format(i, ag.n_atoms, pbc, dim[0], dim[1], dim[2])
            
            f.write(str.encode(info + item + '\n'))
            
            pos = ag.positions
            tys = np.transpose([ag.types]).astype(np.int64)
            ind = np.transpose([np.arange(1, ag.n_atoms + 1)])

            if force_b:
                fce = ag.forces
                if kj2kcal: fce *= 0.239

                data = np.concatenate((
                    ind, tys, pos, fce), axis=-1)

                np.savetxt(f, data, fmt=w)

            else:
                data = np.concatenate((
                    ind, tys, pos), axis=-1)
                
                np.savetxt(f, data, fmt=w)

        f.close()


