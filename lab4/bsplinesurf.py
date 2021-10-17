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

    return ndegree, knots, np.array(vX), np.array(vY)


def sweepbspline(args):
    # log parameters
    print(args)

    d1, knots1, y1, z1 = loadTextFile(args.infile1)
    x1 = np.zeros_like(y1)
    q1 = np.hstack((x1[:, np.newaxis], y1[:, np.newaxis], z1[:, np.newaxis]))
    numv1 = q1.shape[0]
    k1 = d1 + 1

    d2, knots2, x2, y2 = loadTextFile(args.infile2)
    z2 = np.zeros_like(x2)
    q2 = np.hstack((x2[:, np.newaxis], y2[:, np.newaxis], z2[:, np.newaxis]))
    q2 = q2 - q2[0]
    numv2 = q2.shape[0]
    k2 = d2 + 1

    q1 = np.tile(q1[:, np.newaxis, :], (1, numv2, 1))
    q2 = np.tile(q2[np.newaxis, :, :], (numv1, 1, 1))
    q = q1 + q2  # [numv1, numv2]

    for j in range(k2 - 1, numv2):
        q_tmp = np.zeros((numv1, k2, 50, 3))
        u = np.linspace(0, 1, 50) * (knots2[j + 1] - knots2[j]) + knots2[j]
        for i in range(j - k2 + 1, j + 1):
            q_tmp[:, i - (j - k2 + 1), :, :] = q[:, i:i+1, :]
        for r in range(1, k2):
            h = k2 - r
            for i in range(j, j - h, -1):
                alpha = (u - knots2[i]) / (knots2[i + h] - knots2[i])
                q_tmp[:, i - (j - k2 + 1)] = (1.0 - alpha)[np.newaxis, :, np.newaxis] * q_tmp[:, i - 1 - (j - k2 + 1)] + \
                                             alpha[np.newaxis, :, np.newaxis] * q_tmp[:, i - (j - k2 + 1)]
        if j == k2 - 1:
            q_hat = q_tmp[:, k2 - 1, :, :]
        else:
            q_hat = np.concatenate((q_hat, q_tmp[:, k2 - 1, :, :]), axis=1) # [numv1, nump2, 3]

    for j in range(k1 - 1, numv1):
        p_tmp = np.zeros((q_hat.shape[1], k1, 50, 3))
        u = np.linspace(0, 1, 50) * (knots1[j + 1] - knots1[j]) + knots1[j]
        for i in range(j - k1 + 1, j + 1):
            p_tmp[:, i - (j - k1 + 1), :, :] = q_hat[i, :, :][:, np.newaxis, :]
        for r in range(1, k1):
            h = k1 - r
            for i in range(j, j - h, - 1):
                alpha = (u - knots1[i]) / (knots1[i + h] - knots1[i])
                p_tmp[:, i - (j - k1 + 1)] = (1.0 - alpha)[np.newaxis, :, np.newaxis] * p_tmp[:, i - 1 - (j - k1 + 1)] + \
                                             alpha[np.newaxis, :, np.newaxis] * p_tmp[:, i - (j - k1 + 1)]
        if j == k1 - 1:
            p = p_tmp[:, k1 - 1, :, :]
        else:
            p = np.concatenate((p, p_tmp[:, k1 - 1, :, :]), axis=1)

    p = np.transpose(p, (1, 0, 2))

    f = open(args.outfile, "w")
    f.write("OFF\n{:d} {:d} 0\n".format(p.shape[0] * p.shape[1], (p.shape[0] - 1) * (p.shape[1] - 1)))

    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            f.write("{:f} {:f} {:f}\n".format(p[i, j, 0], p[i, j, 1], p[i, j, 2]))

    for i in range(p.shape[0] - 1):
        for j in range(p.shape[1] - 1):
            f.write("4  {:d} {:d} {:d} {:d}\n".format(
                i * p.shape[1] + j, (i + 1) * p.shape[1] + j,
                (i + 1) * p.shape[1] + j + 1, i * p.shape[1] + j + 1))

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
