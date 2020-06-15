from ..common.frame import Frame
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis import Universe

class LAMMPSTRJWriter:
    def __init__(self):
        """ Write lammpstrj file
        smda.LAMMPSTRJWriter().write(
            'CG.lammpstrj', u or u.select_atoms('name SOL SOD CLA'), 
            types = {} or {'SOL':1, 'SOD':2, 'CLA':3}, b=0, e=1e10, skip=1)
        """
        pass


    def write(self, filename, obj, types={}, pbc='pp pp pp',
              b=0, e=1e10, skip=1):
        """
        Parameters
        ----------
        filename [str]
        obj = AtomGroup or Universe
        types = {'SOL': 1, 'SOD': 2, 'CLA': 3} or {}
        pbc = 'pp pp pp'
        b = 0    [float]
        e = 1e10 [float]
        skip = 1 [int]
        
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

        
        ### name2type
        t = 1
        for atom in ag.atoms:
            if atom.name in types.keys():
                atom.type = types[atom.name]
            else:
                atom.type = t
                types[atom.name] = t
                t += 1
        print(types)
        

        ### Write
        bframe, eframe = Frame().frame(u, b, e)
        nframes = int((eframe - bframe)/skip)
        f = open(filename, 'w')
        for i, ts in enumerate(u.trajectory[bframe:eframe+1:skip]):
            if i%100 == 0:
                print("%d/%d processing..." %(i, nframes))
            dim = u.dimensions
            f.write("""ITEM: TIMESTEP
{:d}
ITEM: NUMBER OF ATOMS
{:d}
ITEM: BOX BOUNDS {:s}
0.000000 {:7.3f}
0.000000 {:7.3f}
0.000000 {:7.3f}
ITEM: ATOMS id type xu yu zu
""".format(
    i, ag.n_atoms, pbc, dim[0], dim[1], dim[2]))

            for n, atom in enumerate(ag, 1):
                p = atom.position
                f.write(" {:d} {:d} {:.3f} {:.3f} {:.3f}\n".format(
                    n, atom.type, p[0], p[1], p[2]))
        
        f.close()


