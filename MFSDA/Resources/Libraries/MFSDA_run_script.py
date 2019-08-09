"""
Run script: multivariate functional shape data analysis (MFSDA).
Usage: python ./MFSDA_run_script.py ./data/ ./result/

Author: Chao Huang (chaohuang.stat@gmail.com)
Last update: 2017-08-14
"""
from __future__ import print_function

import sys
import numpy as np
from scipy import stats
from scipy.io import loadmat
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from stat_read_x import read_x
from stat_lpks import lpks
from stat_sif import sif
from stat_wald_ht import wald_ht
from stat_bstrp_pvalue import bstrp_pvalue
import timeit
import MFSDA_stat as mfsda

"""
installed all the libraries above
"""

def run_script(input_dir, output_dir):
    """
    Run the commandline script for MFSDA.

    Args:
        input_dir (str): full path to the data folder
        output_dir (str): full path to the output folder
    """

    """+++++++++++++++++++++++++++++++++++"""
    """Step 1. load dataset """
    print("loading data ......")
    print("+++++++Read the surface shape data+++++++")
    shape_file_name = input_dir + "aligned_shapes.mat"
    mat = loadmat(shape_file_name)
    y_design = mat['aligned_shape']
    n, l, m = y_design.shape
    print("The dimension of shape matrix is " + str(y_design.shape))
    print("+++++++Read the sphere coordinate data+++++++")
    template_file_name = input_dir + "template.mat"
    mat = loadmat(template_file_name)
    coord_mat = mat['template']
    # d = coord_mat.shape[1]
    print("+++++++Read the design matrix+++++++")
    design_data_file_name = input_dir + "design_data.txt"
    design_data_tmp = np.loadtxt(design_data_file_name, delimiter=delimiter)
    if len(design_data_tmp.shape) == 1:
        design_data = np.reshape(design_data_tmp, (design_data_tmp.shape[0], 1))
    else:
        design_data = design_data_tmp

    # read the covariate type
    var_type_file_name = input_dir + "var_type.txt"
    var_type = np.loadtxt(var_type_file_name)
    print("+++++++Construct the design matrix: normalization+++++++")
    x_design = read_x(design_data, var_type)
    p = x_design.shape[1]
    print("The dimension of design matrix is " + str(x_design.shape))

    """+++++++++++++++++++++++++++++++++++"""
    """Step 2. Statistical analysis: including (1) smoothing and (2) hypothesis testing"""
    gpvals, lpvals_fdr, clu_pvals, efit_beta, efity_design, efit_eta = mfsda.run_stats(y_design, coord_mat, design_data, var_type)

    """+++++++++++++++++++++++++++++++++++"""
    """Step3. Save all the results"""
    gpvals_file_name = output_dir + "global_pvalue.txt"
    np.savetxt(gpvals_file_name, gpvals)
    lpvals_fdr_file_name = output_dir + "local_pvalue_fdr.txt"
    np.savetxt(lpvals_fdr_file_name, lpvals_fdr)
    clu_pvals_file_name = output_dir + "cluster_pvalue.txt"
    np.savetxt(clu_pvals_file_name, clu_pvals)


if __name__ == '__main__':
    input_dir0 = sys.argv[1]
    output_dir0 = sys.argv[2]

    start_all = timeit.default_timer()
    run_script(input_dir0, output_dir0)
    stop_all = timeit.default_timer()
    delta_time_all = str(stop_all - start_all)
    print("The total elapsed time is " + delta_time_all)
