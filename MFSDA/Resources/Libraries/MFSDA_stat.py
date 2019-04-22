from __future__ import print_function

import timeit
import numpy as np
from scipy import stats

from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from stat_lpks import lpks
from stat_sif import sif
from stat_wald_ht import wald_ht
from stat_bstrp_pvalue import bstrp_pvalue

from stat_read_x import read_x


def run_stats(y_design, coord_mat, design_data, var_type):
    
    n, l, m = y_design.shape    

    print("+++++++Construct the design matrix: normalization+++++++")
    x_design = read_x(design_data, var_type)
    p = x_design.shape[1]
    print("The dimension of design matrix is ", str(x_design.shape))

    """+++++++++++++++++++++++++++++++++++"""
    """Step 2. Statistical analysis: including (1) smoothing and (2) hypothesis testing"""

    print("+++++++Local linear kernel smoothing+++++++")
    start = timeit.default_timer()
    efit_beta, efity_design, h_opt = lpks(coord_mat, x_design, y_design)
    stop = timeit.default_timer()
    delta_time = str(stop - start)
    # print(h_opt)
    print("Elapsed time is " + delta_time)

    print("+++++++Kernel smoothing (order = 1) for smooth functions (eta)+++++++")
    start = timeit.default_timer()
    resy_design = y_design - efity_design
    print(np.amax(resy_design))
    print(np.amin(resy_design))
    efit_eta, res_eta, esig_eta = sif(coord_mat, resy_design, h_opt)
    print(np.amax(res_eta))
    print(np.amin(res_eta))
    stop = timeit.default_timer()
    delta_time = str(stop - start)
    print("Elapsed time is " + delta_time)

    print("+++++++Hypothesis testing+++++++")
    # hypothesis: beta_pj(d)=0 v.s. beta_pj(d)~=0 for all j and d
    start = timeit.default_timer()
    lpvals = np.zeros((l, p-1))
    lpvals_fdr = np.zeros((l, p-1))
    gpvals = np.zeros((1, p-1))
    clu_pvals = np.zeros((1, p-1))
    areas = np.zeros((1, p-1))
    num_bstrp = 500  # number of bootstrap samples
    thres = 2

    for pp in range(p-1):
        print("Testing whether the covariate " + str(pp+1) + " is zero or not...")
        """ local and global statistics calculation """
        cdesign = np.zeros((1, p))
        cdesign[0, pp+1] = 1
        gstat, lstat = wald_ht(x_design, efit_beta, esig_eta, cdesign)
        lpvals[:, pp] = 1 - np.squeeze(stats.chi2.cdf(lstat, m))
        lpvals_fdr[:, pp] = fdrcorrection0(lpvals[:, pp])[1]
        ind_thres = -np.log10(lpvals[:, pp]) >= thres
        area = np.sum(ind_thres)

        """ Generate random samples and calculate the corresponding statistics and pvalues """
        gpval, clu_pval = bstrp_pvalue(coord_mat, x_design, y_design, cdesign, gstat, num_bstrp, thres, area)

        gpvals[0, pp] = gpval
        areas[0, pp] = area
        clu_pvals[0, pp] = clu_pval
        print("the global p-value for covariate " + str(pp+1) + " is " + str(gpvals[0, pp]) + "...")
        print("the p-value of most significant subregion for covariate " +
              str(pp+1) + " is " + str(clu_pvals[0, pp]) + "...")

    stop = timeit.default_timer()
    delta_time = str(stop - start)
    print("Elapsed time is " + delta_time)

    return gpvals, lpvals_fdr, clu_pvals, efit_beta, efity_design, efit_eta