import sys
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# load vX and vY from a textfile
def loadTextFile(filename):
    vX, vY = [], []
    with open(filename, "r") as fh:
        assert (next(fh).strip("\r\n") == "BSPLINE")
        for line in fh:
            if line[0] != "#":
                ndim, numv, ndegree = [int(tmp) for tmp in line.split()]
                break
        assert (ndim == 2)
        knots = [float(tmp) for tmp in next(fh).split()]
        for line in fh:
            x, y = [float(tmp) for tmp in line.split()]
            vX.append(x)
            vY.append(y)

    assert (numv == len(vX))

    # plt.plot(vX, vY, color="magenta")
    # plt.plot(vX, vY, 'ob', picker=True, pickradius=5)
    k = ndegree + 1

    pX, pY = np.array([]), np.array([])
    for j in range(k - 1, numv):
        q = np.zeros((k, 50, 2))
        u = np.linspace(0, 1, 50) * (knots[j + 1] - knots[j]) + knots[j]
        for i in range(j - k + 1, j + 1):
            q[i - (j - k + 1), :, 0] = vX[i]
            q[i - (j - k + 1), :, 1] = vY[i]
        for r in range(1, k):
            h = k - r
            # construct control points for order h spline
            for i in range(j, j - h, -1):
                alpha = (u - knots[i]) / (knots[i + h] - knots[i])
                q[i - (j - k + 1)] = (1.0 - alpha)[:, np.newaxis] * q[i - 1 - (j - k + 1)] + alpha[:, np.newaxis] * q[i - (j - k + 1)]
        pX = np.concatenate((pX, q[k - 1, :, 0]))
        pY = np.concatenate((pY, q[k - 1, :, 1]))
    # plt.plot(pX, pY, 'b-')
    # plt.draw()
    # plt.show()
    return pX, pY

def sweepbspline(args):
    # log parameters
    print(args)

    pathY, pathZ = loadTextFile(args.pathfile)
    pathX = np.zeros_like(pathY)
    curveX, curveY = loadTextFile(args.curvefile)
    curveZ = np.zeros_like(curveX)

    return

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("--pathfile", required=True, type=str,
                        help="the name of a file containing the the path or trajectory of the sweep curve")
parser.add_argument("--curvefile", required=True, type=str,
                        help="the name of a file containing the definitiion of an open B-spline curve that is “swept” along the pathline")
parser.add_argument("--outfile", required=True, type=str,
                        help="the name of the output file")

if __name__ == '__main__':
    args = parser.parse_args()
    sweepbspline(args)
