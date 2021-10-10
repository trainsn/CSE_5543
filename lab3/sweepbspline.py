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
    alpha = np.hstack((pathX[:, np.newaxis], pathY[:, np.newaxis], pathZ[:, np.newaxis]))
    num_path_points = alpha.shape[0]
    curveX, curveY = loadTextFile(args.curvefile)
    curveZ = np.zeros_like(curveX)
    beta = np.hstack((curveX[:, np.newaxis], curveY[:, np.newaxis], curveZ[:, np.newaxis]))
    beta = beta - beta[0]
    num_curve_points = beta.shape[0]

    alpha = np.tile(alpha[:, np.newaxis, :], (1, num_curve_points, 1))
    beta = np.tile(beta[np.newaxis, :, :], (num_path_points, 1, 1))

    q = alpha + beta

    f = open(args.outfile, "w")
    f.write("OFF\n{:d} {:d} 0\n".format(num_path_points * num_curve_points, (num_path_points - 1) * (num_curve_points - 1)))

    for i in range(num_path_points):
        for j in range(num_curve_points):
            f.write("{:f} {:f} {:f}\n".format(q[i, j, 0], q[i, j, 1], q[i, j, 2]))

    for i in range(num_path_points - 1):
        for j in range(num_curve_points - 1):
            f.write("4  {:d} {:d} {:d} {:d}\n".format(
                i * num_curve_points + j, (i+1) * num_curve_points + j,
                (i+1) * num_curve_points + j+1, i * num_curve_points + j+1))

    return

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("--pathfile", required=True, type=str,

                        help="the name of a file containing the the path or trajectory of the sweep curve")
parser.add_argument("--curvefile", required=True, type=str,
                        help="the name of a file containing the definitiion of an open B-spline curve that is “swept” along the pathline")
parser.add_argument("--outfile", default="out.off", type=str,
                        help="the name of the output file")

if __name__ == '__main__':
    args = parser.parse_args()
    sweepbspline(args)
