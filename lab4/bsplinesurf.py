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

    infile1Y, infile1Z = loadTextFile(args.infile1)
    infile1X = np.zeros_like(infile1Y)
    alpha = np.hstack((infile1X[:, np.newaxis], infile1Y[:, np.newaxis], infile1Z[:, np.newaxis]))
    num_alpha_points = alpha.shape[0]
    infile2X, infile2Y = loadTextFile(args.infile2)
    infile2Z = np.zeros_like(infile2X)
    beta = np.hstack((infile2X[:, np.newaxis], infile2Y[:, np.newaxis], infile2Z[:, np.newaxis]))
    beta = beta - beta[0]
    num_beta_points = beta.shape[0]

    alpha = np.tile(alpha[:, np.newaxis, :], (1, num_beta_points, 1))
    beta = np.tile(beta[np.newaxis, :, :], (num_alpha_points, 1, 1))

    q = alpha + beta

    f = open(args.outfile, "w")
    f.write("OFF\n{:d} {:d} 0\n".format(num_alpha_points * num_beta_points, (num_alpha_points - 1) * (num_beta_points - 1)))

    for i in range(num_alpha_points):
        for j in range(num_beta_points):
            f.write("{:f} {:f} {:f}\n".format(q[i, j, 0], q[i, j, 1], q[i, j, 2]))

    for i in range(num_alpha_points - 1):
        for j in range(num_beta_points - 1):
            f.write("4  {:d} {:d} {:d} {:d}\n".format(
                i * num_beta_points + j, (i+1) * num_beta_points + j,
                (i+1) * num_beta_points + j+1, i * num_beta_points + j+1))

    return

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("--infile1", required=True, type=str,
                        help="the first B-spline curve")
parser.add_argument("--infile2", required=True, type=str,
                        help="the second B-spline curve")
parser.add_argument("--outfile", default="out.off", type=str,
                        help="the name of the output file")

if __name__ == '__main__':
    args = parser.parse_args()
    sweepbspline(args)
