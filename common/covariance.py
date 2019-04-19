import numpy as np

class Covariance:
    def __init__(self):
        pass

    def calculate(self, lag=0, features=None):
        nfeatures = len(features)
        for feature in features:
            assert feature.shape == features[0].shape, "different size"
        features = np.array(features)
        
        dt = features[0][1][0] - features[0][0][0]
        nf = len(features[0][:,0])
        print("nfeatures: ", nfeatures)
        print("dt: %.2f ns" %dt)
        print("nframes: ", nf)

        if lag != 0:
            assert dt < lag, "dt > lag"

        if lag % dt < 0.0001:
            print("good lag time")
        else:
            print("bad lag time")
            lag = (lag // dt) * dt
            print("adjusting lag time to %.2f ns" %lag)

        df = int(lag//dt)
        
        n = 0
        covar = np.zeros((nfeatures, nfeatures))
        print("==================================")
        print("start covariance calculation")
        while (n + df) < nf:
            X = features[:, n, 1]
            Y = features[:, n + df, 1]
            X1 = X - np.average(X)
            Y1 = Y - np.average(Y)
            covar += np.outer(X1, Y1)
            n += 1
        covar /= n

        return covar


