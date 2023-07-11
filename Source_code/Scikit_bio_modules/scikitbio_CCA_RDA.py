# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.linalg import svd, lstsq

from ordination_results import OrdinationResults
from utils import corr, svd_rank, scale
#from skbio.util._decorator import experimental

#@experimental(as_of="0.4.0")
def cca(y, x, scaling=1):
    r"""Compute canonical (also known as constrained) correspondence
    analysis.

    Canonical (or constrained) correspondence analysis is a
    multivariate ordination technique. It appeared in community
    ecology [1]_ and relates community composition to the variation in
    the environment (or in other factors). It works from data on
    abundances or counts of samples and constraints variables,
    and outputs ordination axes that maximize sample separation among species.

    It is better suited to extract the niches of taxa than linear
    multivariate methods because it assumes unimodal response curves
    (habitat preferences are often unimodal functions of habitat
    variables [2]_).

    As more environmental variables are added, the result gets more
    similar to unconstrained ordination, so only the variables that
    are deemed explanatory should be included in the analysis.

    Parameters
    ----------
    y : DataFrame
        Samples by features table (n, m)
    x : DataFrame
        Samples by constraints table (n, q)
    scaling : int, {1, 2}, optional
        Scaling type 1 maintains :math:`\chi^2` distances between rows.
        Scaling type 2 preserves :math:`\chi^2` distances between columns.
        For a more detailed explanation of the interpretation, check Legendre &
        Legendre 1998, section 9.4.3.

    Returns
    -------
    OrdinationResults
        Object that stores the cca results.

    Raises
    ------
    ValueError
        If `x` and `y` have different number of rows
        If `y` contains negative values
        If `y` contains a row of only 0's.
    NotImplementedError
        If scaling is not 1 or 2.

    See Also
    --------
    ca
    rda
    OrdinationResults

    Notes
    -----
    The algorithm is based on [3]_, \S 11.2, and is expected to give
    the same results as ``cca(y, x)`` in R's package vegan, except
    that this implementation won't drop constraining variables due to
    perfect collinearity: the user needs to choose which ones to
    input.

    Canonical *correspondence* analysis shouldn't be confused with
    canonical *correlation* analysis (CCorA, but sometimes called
    CCA), a different technique to search for multivariate
    relationships between two datasets. Canonical correlation analysis
    is a statistical tool that, given two vectors of random variables,
    finds linear combinations that have maximum correlation with each
    other. In some sense, it assumes linear responses of "species" to
    "environmental variables" and is not well suited to analyze
    ecological data.

    References
    ----------
    .. [1] Cajo J. F. Ter Braak, "Canonical Correspondence Analysis: A
        New Eigenvector Technique for Multivariate Direct Gradient
        Analysis", Ecology 67.5 (1986), pp. 1167-1179.

    .. [2] Cajo J.F. Braak and Piet F.M. Verdonschot, "Canonical
        correspondence analysis and related multivariate methods in
        aquatic ecology", Aquatic Sciences 57.3 (1995), pp. 255-289.

    .. [3] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.

    """
    Y = y.values
    X = x.values

    # Perform parameter sanity checks
    if X.shape[0] != Y.shape[0]:
        raise ValueError("The samples by features table 'y' and the samples by"
                         " constraints table 'x' must have the same number of "
                         " rows. 'y': {0} 'x': {1}".format(X.shape[0],
                                                           Y.shape[0]))
    if Y.min() < 0:
        raise ValueError(
            "The samples by features table 'y' must be nonnegative")
    row_max = Y.max(axis=1)
    if np.any(row_max <= 0):
        # Or else the lstsq call to compute Y_hat breaks
        raise ValueError("The samples by features table 'y' cannot contain a "
                         "row with only 0's")
    if scaling not in {1, 2}:
        raise NotImplementedError(
            "Scaling {0} not implemented.".format(scaling))

    # Step 1 (similar to Pearson chi-square statistic)
    grand_total = Y.sum()
    Q = Y / grand_total  # Relative frequencies of Y (contingency table)

    # Features and sample weights (marginal totals)
    column_marginals = Q.sum(axis=0)
    row_marginals = Q.sum(axis=1)

    # Formula 9.32 in Lagrange & Lagrange (1998). Notice that it's an
    # scaled version of the contribution of each cell towards Pearson
    # chi-square statistic.
    expected = np.outer(row_marginals, column_marginals)
    Q_bar = (Q - expected) / np.sqrt(expected)

    # Step 2. Standardize columns of X with respect to sample weights,
    # using the maximum likelihood variance estimator (Legendre &
    # Legendre 1998, p. 595)
    X = scale(X, weights=row_marginals, ddof=0)

    # Step 3. Weighted multiple regression.
    X_weighted = row_marginals[:, None]**0.5 * X
    B, _, rank_lstsq, _ = lstsq(X_weighted, Q_bar)
    Y_hat = X_weighted.dot(B)
    Y_res = Q_bar - Y_hat

    # Step 4. Eigenvalue decomposition
    u, s, vt = svd(Y_hat, full_matrices=False)
    rank = svd_rank(Y_hat.shape, s)
    s = s[:rank]
    u = u[:, :rank]
    vt = vt[:rank]
    U = vt.T

    # Step 5. Eq. 9.38
    U_hat = Q_bar.dot(U) * s**-1

    # Residuals analysis
    u_res, s_res, vt_res = svd(Y_res, full_matrices=False)
    rank = svd_rank(Y_res.shape, s_res)
    s_res = s_res[:rank]
    u_res = u_res[:, :rank]
    vt_res = vt_res[:rank]

    U_res = vt_res.T
    U_hat_res = Y_res.dot(U_res) * s_res**-1

    eigenvalues = np.r_[s, s_res]**2

    # Scalings (p. 596 L&L 1998):
    # feature scores, scaling 1
    V = (column_marginals**-0.5)[:, None] * U

    # sample scores, scaling 2
    V_hat = (row_marginals**-0.5)[:, None] * U_hat

    # sample scores, scaling 1
    F = V_hat * s

    # feature scores, scaling 2
    F_hat = V * s

    # Sample scores which are linear combinations of constraint
    # variables
    Z_scaling1 = ((row_marginals**-0.5)[:, None] *
                  Y_hat.dot(U))
    Z_scaling2 = Z_scaling1 * s**-1

    # Feature residual scores, scaling 1
    V_res = (column_marginals**-0.5)[:, None] * U_res

    # Sample residual scores, scaling 2
    V_hat_res = (row_marginals**-0.5)[:, None] * U_hat_res

    # Sample residual scores, scaling 1
    F_res = V_hat_res * s_res

    # Feature residual scores, scaling 2
    F_hat_res = V_res * s_res

    eigvals = eigenvalues
    if scaling == 1:
        features_scores = np.hstack((V, V_res))
        sample_scores = np.hstack((F, F_res))
        sample_constraints = np.hstack((Z_scaling1, F_res))
    elif scaling == 2:
        features_scores = np.hstack((F_hat, F_hat_res))
        sample_scores = np.hstack((V_hat, V_hat_res))
        sample_constraints = np.hstack((Z_scaling2, V_hat_res))

    biplot_scores = corr(X_weighted, u)

    pc_ids = ['CCA%d' % (i+1) for i in range(len(eigenvalues))]
    sample_ids = y.index
    feature_ids = y.columns
    eigvals = pd.Series(eigenvalues, index=pc_ids)
    samples = pd.DataFrame(sample_scores,
                           columns=pc_ids, index=sample_ids)
    features = pd.DataFrame(features_scores,
                            columns=pc_ids, index=feature_ids)

    biplot_scores = pd.DataFrame(biplot_scores,
                                 index=x.columns,
                                 columns=pc_ids[:biplot_scores.shape[1]])
    sample_constraints = pd.DataFrame(sample_constraints,
                                      index=sample_ids, columns=pc_ids)

    return OrdinationResults(
        "CCA", "Canonical Correspondence Analysis", eigvals, samples,
        features=features, biplot_scores=biplot_scores,
        sample_constraints=sample_constraints,
        proportion_explained=eigvals / eigvals.sum())

def rda(y, x, scale_Y=False, scaling=1):
    r"""Compute redundancy analysis, a type of canonical analysis.
    It is related to PCA and multiple regression because the explained
    variables `y` are fitted to the explanatory variables `x` and PCA
    is then performed on the fitted values. A similar process is
    performed on the residuals.
    RDA should be chosen if the studied gradient is small, and CCA
    when it's large, so that the contingency table is sparse.
    Parameters
    ----------
    y : pd.DataFrame
        :math:`n \times p` response matrix, where :math:`n` is the number
        of samples and :math:`p` is the number of features. Its columns
        need be dimensionally homogeneous (or you can set `scale_Y=True`).
        This matrix is also referred to as the community matrix that
        commonly stores information about species abundances
    x : pd.DataFrame
        :math:`n \times m, n \geq m` matrix of explanatory
        variables, where :math:`n` is the number of samples and
        :math:`m` is the number of metadata variables. Its columns
        need not be standardized, but doing so turns regression
        coefficients into standard regression coefficients.
    scale_Y : bool, optional
        Controls whether the response matrix columns are scaled to
        have unit standard deviation. Defaults to `False`.
    scaling : int
        Scaling type 1 produces a distance biplot. It focuses on
        the ordination of rows (samples) because their transformed
        distances approximate their original euclidean
        distances. Especially interesting when most explanatory
        variables are binary.
        Scaling type 2 produces a correlation biplot. It focuses
        on the relationships among explained variables (`y`). It
        is interpreted like scaling type 1, but taking into
        account that distances between objects don't approximate
        their euclidean distances.
        See more details about distance and correlation biplots in
        [1]_, \S 9.1.4.
    Returns
    -------
    OrdinationResults
        Object that stores the computed eigenvalues, the
        proportion explained by each of them (per unit),
        transformed coordinates for feature and samples, biplot
        scores, sample constraints, etc.
    See Also
    --------
    ca
    cca
    OrdinationResults
    Notes
    -----
    The algorithm is based on [1]_, \S 11.1, and is expected to
    give the same results as ``rda(y, x)`` in R's package vegan.
    The eigenvalues reported in vegan are re-normalized to
    :math:`\sqrt{\frac{s}{n-1}}` `n` is the number of samples,
    and `s` is the original eigenvalues. Here we will only return
    the original eigenvalues, as recommended in [1]_.
    References
    ----------
    .. [1] Legendre P. and Legendre L. 1998. Numerical
       Ecology. Elsevier, Amsterdam.
    """
    Y = y.values
    X = x.values

    n, p = y.shape
    n_, m = x.shape
    if n != n_:
        raise ValueError(
            "Both data matrices must have the same number of rows.")
    if n < m:
        # Mmm actually vegan is able to do this case, too
        raise ValueError(
            "Explanatory variables cannot have less rows than columns.")

    sample_ids = y.index
    feature_ids = y.columns
    # Centre response variables (they must be dimensionally
    # homogeneous)
    Y = scale(Y, with_std=scale_Y)
    # Centre explanatory variables
    X = scale(X, with_std=False)

    # Distribution of variables should be examined and transformed
    # if necessary (see paragraph 4 in p. 580 L&L 1998)

    # Compute Y_hat (fitted values by multivariate linear
    # regression, that is, linear least squares). Formula 11.6 in
    # L&L 1998 involves solving the normal equations, but that fails
    # when cond(X) ~ eps**(-0.5). A more expensive but much more
    # stable solution (fails when cond(X) ~ eps**-1) is computed
    # using the QR decomposition of X = QR:
    # (11.6) Y_hat = X [X' X]^{-1} X' Y
    #              = QR [R'Q' QR]^{-1} R'Q' Y
    #              = QR [R' R]^{-1} R'Q' Y
    #              = QR R^{-1} R'^{-1} R' Q' Y
    #              = Q Q' Y
    # and B (matrix of regression coefficients)
    # (11.4) B = [X' X]^{-1} X' Y
    #          = R^{-1} R'^{-1} R' Q' Y
    #          = R^{-1} Q'
    # Q, R = np.linalg.qr(X)
    # Y_hat = Q.dot(Q.T).dot(Y)
    # B = scipy.linalg.solve_triangular(R, Q.T.dot(Y))
    # This works provided X has full rank. When not, you can still
    # fix it using R's pseudoinverse or partitioning R. To avoid any
    # issues, like the numerical instability when trying to
    # reproduce an example in L&L where X was rank-deficient, we'll
    # just use `np.linalg.lstsq`, which uses the SVD decomposition
    # under the hood and so it's also more expensive.
    B, _, rank_X, _ = lstsq(X, Y)
    Y_hat = X.dot(B)
    # Now let's perform PCA on the fitted values from the multiple
    # regression
    u, s, vt = svd(Y_hat, full_matrices=False)
    # vt are the right eigenvectors, which is what we need to
    # perform PCA. That is, we're changing points in Y_hat from the
    # canonical basis to the orthonormal basis given by the right
    # eigenvectors of Y_hat (or equivalently, the eigenvectors of
    # the covariance matrix Y_hat.T.dot(Y_hat))
    # See 3) in p. 583 in L&L 1998
    rank = svd_rank(Y_hat.shape, s)
    # Theoretically, there're at most min(p, m, n - 1) non-zero eigenvalues

    U = vt[:rank].T  # U as in Fig. 11.2

    # Ordination in the space of response variables. Its columns are
    # sample scores. (Eq. 11.12)
    F = Y.dot(U)
    # Ordination in the space of explanatory variables. Its columns
    # are fitted sample scores. (Eq. 11.13)
    Z = Y_hat.dot(U)

    # Canonical coefficients (formula 11.14)
    # C = B.dot(U)  # Not used

    Y_res = Y - Y_hat
    # PCA on the residuals
    u_res, s_res, vt_res = svd(Y_res, full_matrices=False)
    # See 9) in p. 587 in L&L 1998
    rank_res = svd_rank(Y_res.shape, s_res)
    # Theoretically, there're at most min(p, n - 1) non-zero eigenvalues as

    U_res = vt_res[:rank_res].T
    F_res = Y_res.dot(U_res)  # Ordination in the space of residuals

    eigenvalues = np.r_[s[:rank], s_res[:rank_res]]

    # Compute scores
    if scaling not in {1, 2}:
        raise NotImplementedError("Only scalings 1, 2 available for RDA.")
    # According to the vegan-FAQ.pdf, the scaling factor for scores
    # is (notice that L&L 1998 says in p. 586 that such scaling
    # doesn't affect the interpretation of a biplot):
    eigvals = pd.Series(
        eigenvalues, index=['RDA%d' % (i+1) for i in range(len(eigenvalues))])
    const = np.sum(eigenvalues**2)**0.25
    if scaling == 1:
        scaling_factor = const
    elif scaling == 2:
        scaling_factor = eigenvalues / const
    feature_scores = np.hstack((U, U_res)) * scaling_factor
    sample_scores = np.hstack((F, F_res)) / scaling_factor

    feature_scores = pd.DataFrame(
        feature_scores, index=feature_ids,
        columns=['RDA%d' % (i+1) for i in range(feature_scores.shape[1])])
    sample_scores = pd.DataFrame(
        sample_scores, index=sample_ids,
        columns=['RDA%d' % (i+1) for i in range(sample_scores.shape[1])])
    # TODO not yet used/displayed
    sample_constraints = np.hstack((Z, F_res)) / scaling_factor
    sample_constraints = pd.DataFrame(
        sample_constraints, index=sample_ids,
        columns=['RDA%d' % (i+1) for i in range(sample_constraints.shape[1])])
    # Vegan seems to compute them as corr(X[:, :rank_X],
    # u) but I don't think that's a good idea. In fact, if
    # you take the example shown in Figure 11.3 in L&L 1998 you
    # can see that there's an arrow for each of the 4
    # environmental variables (depth, coral, sand, other) even if
    # other = not(coral or sand)
    biplot_scores = corr(X, u)
    biplot_scores = pd.DataFrame(
        biplot_scores, index=x.columns,
        columns=['RDA%d' % (i+1) for i in range(biplot_scores.shape[1])])
    # The "Correlations of environmental variables with sample
    # scores" from table 11.4 are quite similar to vegan's biplot
    # scores, but they're computed like this:
    # corr(X, F))
    p_explained = pd.Series(
        eigenvalues / eigenvalues.sum(),
        index=['RDA%d' % (i+1) for i in range(len(eigenvalues))])
    return OrdinationResults('RDA', 'Redundancy Analysis',
                             eigvals=eigvals,
                             proportion_explained=p_explained,
                             features=feature_scores,
                             samples=sample_scores,
                             biplot_scores=biplot_scores,
                             sample_constraints=sample_constraints)