# netcompare

Test if two groups of brain networks are significantly non-overlapping.

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
