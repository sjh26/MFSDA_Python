"""
Run script: multivariate functional shape data analysis (MFSDA).
Usage: python ./MFSDA_run_script.py ./data/ ./result/

Author: Chao Huang (chaohuang.stat@gmail.com)
Last update: 2017-08-14
"""

import sys
import numpy as np
from scipy import stats
#from scipy.io import loadmat
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from stat_read_x import read_x
from stat_lpks import lpks
from stat_sif import sif
from stat_wald_ht import wald_ht
from stat_bstrp_pvalue import bstrp_pvalue
# from sklearn.cluster import KMeans
# from scipy.cluster.vq import kmeans2
# from stat_gap import gap
import MFSDA_stat as mfsda
import timeit
import vtk
import argparse
import os
"""
installed all the libraries above
"""

parser = argparse.ArgumentParser(description='Multivariate Functional Shape Data Analysis (MFSDA)')
parser.add_argument('--shapeData', type=str, help='Text file list with vtk filenames, 1 file per line', required=True)
parser.add_argument('--coordData', type=str, help='filename, .vtk shape template', required=True)
parser.add_argument('--covariate', type=str, help='filename, with covariates dim = n x p0 (comma separated or tabulation, without header)', required=True)
#parser.add_argument('--covariateInterest', type=str, help='filename (dim = 1xp0 vector comma separated, 1 or 0 value to indicate/activate covariate of interest)', required=True)
parser.add_argument('--covariateType', help='filename, dim=1xsum(covariateInterest) vector comma separated, 1 or 0 to indicate type of covariate double or int', required=True)
parser.add_argument('--outputDir', help='output directory', default='./out')

def run_script(args):
    """
    Run the commandline script for MFSDA.
    """
    """+++++++++++++++++++++++++++++++++++"""
    """Step 1. load dataset """

    print("loading data ......")
    print("+++++++Read the surface shape data+++++++")    

    fh = open(args.shapeData, 'rU')     

    y_design = []
    nshape = 0
    numpoints = -1

    for vtkfilename in fh.readlines():
        print("Reading", vtkfilename)
        vtkfilename = vtkfilename.rstrip()        
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtkfilename)
        reader.Update()
        shapedata = reader.GetOutput()
        shapedatapoints = shapedata.GetPoints()

        y_design.append([])

        if numpoints == -1:
            numpoints = shapedatapoints.GetNumberOfPoints()

        if numpoints != shapedatapoints.GetNumberOfPoints():
            print("WARNING! The number of points is not the same for the shape:", vtkfilename)

        for i in range(shapedatapoints.GetNumberOfPoints()):
            p = shapedatapoints.GetPoint(i)
            y_design[nshape].append(p)

        nshape += 1

    y_design = np.array(y_design)
    y_design.reshape(nshape, numpoints, 3)

    y_design = np.array(y_design)
    y_design.reshape(nshape, numpoints, 3) 
    print("The dimension of shape matrix is " + str(y_design.shape))

    print("+++++++Read the sphere coordinate data+++++++")
    print("Reading", args.coordData)    
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(args.coordData)
    reader.Update()
    coordData = reader.GetOutput()
    shapedatapoints = coordData.GetPoints()

    if numpoints != shapedatapoints.GetNumberOfPoints():
        print("WARNING! The template does not have the same number of points as the shapes")

    coord_mat = []
    for i in range(shapedatapoints.GetNumberOfPoints()):
        p = shapedatapoints.GetPoint(i)
        coord_mat.append(p)

    coord_mat = np.array(coord_mat)        

    print("+++++++Read the design matrix+++++++")
    design_data_file_name = args.covariate
    _ , file_extension = os.path.splitext(design_data_file_name)

    delimiter = ' '
    if file_extension == '.csv':
        delimiter=','
    
    design_data = np.loadtxt(design_data_file_name, delimiter=delimiter)

    # read the covariate type
    var_type_file_name = args.covariateType
    _ , file_extension = os.path.splitext(var_type_file_name)

    delimiter = ' '
    if file_extension == '.csv':
        delimiter=','

    var_type = np.loadtxt(var_type_file_name, delimiter=delimiter)

    """+++++++++++++++++++++++++++++++++++"""
    """Step 2. Statistical analysis: including (1) smoothing and (2) hypothesis testing"""

    gpvals, lpvals_fdr, clu_pvals = mfsda.run_stats(y_design, coord_mat, design_data, var_type)

    """+++++++++++++++++++++++++++++++++++"""
    """Step3. Save all the results"""
    gpvals_file_name = os.path.join(args.outputDir, "global_pvalue.txt")
    np.savetxt(gpvals_file_name, gpvals)
    lpvals_fdr_file_name = os.path.join(args.outputDir, "local_pvalue_fdr.txt")
    np.savetxt(lpvals_fdr_file_name, lpvals_fdr)
    clu_pvals_file_name = os.path.join(args.outputDir, "cluster_pvalue.txt")
    np.savetxt(clu_pvals_file_name, clu_pvals)

if __name__ == '__main__':
        
    args = parser.parse_args()

    start_all = timeit.default_timer()
    run_script(args)
    stop_all = timeit.default_timer()
    delta_time_all = str(stop_all - start_all)
    print("The total elapsed time is " + delta_time_all)