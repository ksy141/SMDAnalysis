from __future__ import absolute_import, division, print_function
import numpy as np
from ..common.frame import Frame

class Membrane:
    def __init__(self):
        pass

    def dimensions(self, u, b=0, e=10000, skip=1):
        bframe, eframe = Frame().frame(u, b, e)
        t = []; x = []; y = []; z = []
        for ts in u.trajectory[bframe:eframe:skip]:
            pbcx, pbcy, pbcz = u.dimensions[0:3]
            t.append(ts.time/1000)
            x.append(pbcx/10)
            y.append(pbcy/10)
            z.append(pbcz/10)

        t = np.array(t); x = np.array(x); y = np.array(y); z = np.array(z)

        return t, x, y, z


    def apl(self, u, b=0, e=10000, skip=1):
        N = len(u.select_atoms("name P"))/2
        t, x, y, z = self.dimensions(u, b, e, skip)
        return t, x*y/N





