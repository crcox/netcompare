#!/usr/bin/env python
"""Test if two groups of brain networks are significantly non-overlapping.

This program implements the non-parametric permutation-based method described
in:

Simpson, S. L., Lyday, R. G., Hayasaka, S., Marsh, A. P., & Laurienti, P. J.
    (2013). A permutation testing framework to compare groups of brain
    networks.  Frontiers in Computational Neuroscience, 7, 171.
    http://doi.org/10.3389/fncom.2013.00171

In brief, the procedure involves computing the Jaccard index between every pair
of networks, within and between groups. The statistic of interest is the mean
Jaccard index within groups, j_w, divided by the mean Jaccard index between
groups, j_b, so r = j_w / j_b.

This statistic is then recomputed after permuting the group labels, so that the
grouping is essentially at random. This is done as many times as you specify
(default 10,000) to estimate an empirical null distribution for the test
statistic. A p-value is determined by the portion of random re-groupings that
yield a larger statistic than the true grouping.

NB!!! The current implementation assumes a paired test. This has a few
consequences. The first is that A and B must contain the same number of
volumes, and A and B are paired by index (e.g., the first volume in each is
subject 1, and so on). To respect that volumes are paired this way, the
permutation procedure is somewhat constrained. Condition reassignments are done
by flipping the labels within pairs at random---A[1] and B[1] will always be in
different groups, the group labels might be flipped on any given permutation.
"""
import os
import numpy as np
import nibabel as nib
import argparse
import itertools
import numpy.random as random

parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('A', type=str, help='Path to nifti file containing sample of networks of type A.')
parser.add_argument('B', type=str, help='Path to nifti file containing sample of networks of type B.')
parser.add_argument('--nperms', type=int, default=10000, help='Number of permutations to establish empirical null.')
parser.add_argument('--alpha', type=float, default=0.05, help='Threshold of statistical significance.')
parser.add_argument('--permlog', type=str, default='permstat.1D', help='Name of logfile that will contain Jaccard stat at each permutation.')
parser.add_argument('--histfig', action='store_true', help='Generate a histogram displaying the permutation distribution, with true statistic overlayed.')
args = parser.parse_args()

if args.histfig:
    import matplotlib.pyplot as plt

Ahandle = nib.load(args.A)
Bhandle = nib.load(args.B)

# Load 4D data x,y,z,subject
A = np.squeeze(Ahandle.get_data())
B = np.squeeze(Bhandle.get_data())

X = np.concatenate((A,B),axis=3)
R = np.zeros((X.shape[3],X.shape[3]))
# Fill in symmetric matrix
for i,j in itertools.permutations(range(X.shape[3]),2):
    p = np.abs(X[:,:,:,i]) > 0
    q = np.abs(X[:,:,:,j]) > 0
    R[i,j] = np.sum(p & q) / float(np.sum(p | q))

m = A.shape[3]
n = m + B.shape[3]

Nw = ((A.shape[3]**2) - A.shape[3]) / 2 + ((B.shape[3]**2) - B.shape[3]) / 2
Nb = A.shape[3] * B.shape[3]

# Compute Jaccard statistic
within = (np.sum(np.tril(R[1:m,1:m],-1)) + np.sum(np.tril(R[m:n,m:n]))) / float(Nw)
between = np.sum(R[1:m,m:n]) / float(Nb)
r = within / between

# Permutations
permlist = np.zeros(args.nperms)
for i in xrange(args.nperms):
    flip_a = list(np.nonzero(random.rand(A.shape[3]) > 0.5)[0])
    flip_b = [x + m for x in flip_a]
    flip = range(X.shape[3])
    for (x,y) in zip(flip_a,flip_b):
        flip[x] = y
        flip[y] = x
    Rperm = np.copy(R)
    Rperm[:,:] = Rperm[flip,:]
    Rperm[:,:] = Rperm[:,flip]
    within = (np.sum(np.tril(Rperm[1:m,1:m],-1)) + np.sum(np.tril(Rperm[m:n,m:n]))) / float(Nw)
    between = np.sum(Rperm[1:m,m:n]) / float(Nb)
    permlist[i] = within / between

with open(args.permlog, 'w') as f:
    for x in permlist:
        f.write("{:.8f}\n".format(x))

# Stats
p = np.sum(permlist > r) / args.nperms
h = p < args.alpha

print "Jaccard Statistic: {:.4f}".format(r)
print "p-value: {:.4f} (alpha: {:.4f})".format(p, args.alpha)

if h:
    print "Reject the null hypothesis (conditions differ)."
else:
    print "Do not reject the null hypotheses (conditions do not differ)."

if args.histfig:
    plt.hist(permlist)
    plt.axvline(r, color='r')
    plt.title("Empirical null distribution, with true value")
    plt.xlabel("Jaccard statistic")
    plt.ylabel("Frequency")
    plt.savefig('permdist.pdf')
