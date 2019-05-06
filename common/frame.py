import numpy as np

class Frame:

    '''
    input:
    b = time (ns)
    e = time (ns)
    
    output:
    bframe = (int) frame corresponding to b (ns)
    eframe = (int) frame corresponding to e (ns)
    
    use:
    bframe, eframe = smda.Frame().frame(u, b, e)
    for ts in u.trajectory[bframe, eframe]:
        DO!
    '''

    def __init__(self):
        pass

    def frame(self, u, b, e):
        sb = u.trajectory[0].time
        se = u.trajectory[-1].time
        sf = len(u.trajectory)
        if se < e*1000: e = se/1000
        if sb > b*1000: b = sb/1000
        assert b  <= e,  "b > e"
        assert b  <= se, "b > se"
        
        bframe = None
        eframe = None

        for ts in u.trajectory:
            if ts.time >= b*1000:
                bframe = ts.frame
                break

        for ts in u.trajectory:
            if ts.time >= e*1000:
                eframe = ts.frame
                break
        
        #bframe = int((b*1000 - sb)/(se-sb) * sf)
        #eframe = int((e*1000 - sb)/(se-sb) * sf)
        return bframe, eframe


 
