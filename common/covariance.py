import numpy as np
from decimal import Decimal

class Covariance:
    def __init__(self):
        pass

    def calculate(self, lag=0, features=None):
        print("==================================")
        print("start covariance calculation")
        nfeatures = len(features)
        for feature in features:
            assert feature.shape == features[0].shape, "different size"
        features = np.array(features)
        
        dt = features[0][2][0] - features[0][1][0]
        nf = len(features[0][:,0])
        print("nfeatures: ", nfeatures)
        print("dt: %.2f ns" %dt)
        print("lag: %.2f ns" %lag)
        print("nframes: ", nf)

        if lag != 0:
            assert dt < lag, "dt > lag"

        if Decimal('%f' %lag) % Decimal('%f' %dt) < 0.001:
            print("good lag time")
        else:
            print("bad lag time")
            lag = (Decimal('%f' %lag) // Decimal('%f' %dt)) * dt
            print("adjusting lag time to %.2f ns" %lag)

        df = int(Decimal('%f' %lag)//Decimal('%f' %dt))
        
        n = 0
        covar = np.zeros((nfeatures, nfeatures))
        aver = np.average(features[:, :, 1], axis=1)
        print("<X> = ", aver)
        while (n + df) < nf:
            X = features[:, n, 1]
            Y = features[:, n + df, 1]
            covar += np.outer(X, Y)
            n += 1
        covar /= n
        covar -= np.outer(aver, aver)

        return covar


