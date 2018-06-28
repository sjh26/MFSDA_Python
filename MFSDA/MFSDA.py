#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
#import inputData
import sys
import numpy as np
from scipy import stats
import pandas as pd
#from scipy.io import loadmat
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
#from stat_read_x import read_x
#from stat_lpks import lpks
#from stat_sif import sif
#from stat_wald_ht import wald_ht
#from stat_bstrp_pvalue import bstrp_pvalue
import inspect
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
sys.path.append(cmd_subfolder)
from MFSDA_run import run_script
from createShapes import run_script2
# from sklearn.cluster import KMeans
# from scipy.cluster.vq import kmeans2
# from stat_gap import gap
import MFSDA_stat as mfsda
import timeit
import argparse
import json
import os.path
import csv

#
# MFSDA
#

class MFSDA(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
        https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
        """
    
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "MFSDA" # TODO make this more human readable by adding spaces
        self.parent.categories = ["Pvalues"]
        self.parent.dependencies = []
        self.parent.contributors = ["Loic Michoud"] # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
            This is an example of scripted loadable module bundled in an extension.
            It performs a simple thresholding on the input volume and optionally captures a screenshot.
            """
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
            This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
            and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
            """ # replace with organization, grant and thanks.

#
# MFSDAWidget
#

class MFSDAWidget(ScriptedLoadableModuleWidget):
    """Uses ScriptedLoadableModuleWidget base class, available at:
        https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
        """
    
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = MFSDALogic(self)
        self.moduleName = 'MFSDA'
        scriptedModulesPath = eval('slicer.modules.%s.path' % self.moduleName.lower())
        scriptedModulesPath = os.path.dirname(scriptedModulesPath)
        path = os.path.join(scriptedModulesPath, 'Resources', 'UI', '%s.ui' % self.moduleName)
        self.widget = slicer.util.loadUI(path)
        self.layout.addWidget(self.widget)
        # Instantiate and connect widgets ...
        self.stateCSVMeansShape = False
        self.dictShapeModels = dict()
        
        self.lineEdit_vtk = self.logic.get('lineEdit_vtk')
        self.lineEdit_template = self.logic.get('lineEdit_template')
        self.lineEdit_csv = self.logic.get('lineEdit_csv')
        self.lineEdit_output = self.logic.get('lineEdit_output')
        self.lineEdit_ShapePvalue=self.logic.get('lineEdit_ShapePvalue')
        self.lineEdit_csv.connect('directoryChanged(const QString)', self.onCSVFile)
        self.lineEdit_template.connect('currentPathChanged(const QString)', self.onCSVFile)
        self.lineEdit_vtk.connect('currentPathChanged(const QString)', self.onCSVFile)
        self.lineEdit_output.connect('directoryChanged (const QString)', self.onCSVFile)
        self.lineEdit_ShapePvalue.connect('directoryChanged (const QString)', self.onCSVFile)
    
    
    
    
    def onCSVFile(self):
        self.dictShapeModels = dict()
        if not os.path.exists(self.lineEdit_csv.currentPath):
            self.stateCSVMeansShape = False
            return
        if not os.path.exists(self.lineEdit_template.currentPath):
            self.stateCSVMeansShape = False
            return
        if not os.path.exists(self.lineEdit_vtk.currentPath):
            self.stateCSVMeansShape = False
            return
        if not os.path.exists(self.lineEdit_output.directory):
            self.stateCSVMeansShape = False
            return
        if not os.path.exists(self.lineEdit_ShapePvalue.currentPath):
            self.stateCSVMeansShape = False
            return
        condition1 = self.logic.checkExtension(self.lineEdit_vtk.currentPath, ".txt")
        condition2 = self.logic.checkExtension(self.lineEdit_template.currentPath, ".vtk")
        condition4 = self.logic.checkExtension(self.lineEdit_csv.currentPath, ".csv")
        condition5= self.lineEdit_output.directory!='.'
        condition6 = self.logic.checkExtension(self.lineEdit_ShapePvalue.currentPath, ".vtk")
        if not condition1:
            self.lineEdit_vtk.setCurrentPath(" ")
            self.stateCSVDataset = False
            return
        if not condition2:
            self.lineEdit_template.setCurrentPath(" ")
            self.stateCSVDataset = False
            return
        if not condition4:
            self.lineEdit_csv.setCurrentPath(" ")
            self.stateCSVDataset = False
            return
        if not condition5:
            self.stateCSVDataset = False
            return
        if not condition6:
            self.lineEdit_ShapePvalue.setCurrentPath(" ")
            self.stateCSVDataset = False
            return
        PathOutput=os.path.dirname(self.lineEdit_vtk.currentPath)+'/'


        covariates_sa = pd.read_csv(self.lineEdit_csv.currentPath)
        covariates = np.array(covariates_sa.axes[1])
        if type(covariates[0]) == str:
            covariates_sa = pd.read_csv(self.lineEdit_csv.currentPath, usecols=covariates)
        covariates_sa = np.array(covariates_sa)
        #######################################################################
        # Reconstruct cvs file without the name of the covariates
        CSVFile=open(self.lineEdit_output.directory+'/valuesCsv.csv', 'wb')
        with CSVFile:
            writer = csv.writer(CSVFile)
            for i in range(0,len(covariates_sa)):
                writer.writerow(covariates_sa[i])
        CSVFile.close()
        #######################################################################
        # Create the covariateType file (1 or 0 to indicate type of covariate double or int)
        covariateType = np.ones((1,len(covariates)))
        for j in range(0,len(covariates)-1):
            if isinstance(covariates_sa[0][j], int):
                covariateType[j] = 0
        if os.path.exists(self.lineEdit_output.directory+'/CovariateType.csv'):
            os.remove(self.lineEdit_output.directory+'/CovariateType.csv')
        CSVFile=open(self.lineEdit_output.directory+'/CovariateType.csv', 'wb')
        with CSVFile:
            writer = csv.writer(CSVFile)
            writer.writerow(covariateType[0])
        CSVFile.close()
        ######################################################################
        # Rewrite txt file with full path to the vtk files
        tempPath = self.lineEdit_output.directory+'/tempText.txt'
        if os.path.exists(tempPath):
            os.remove(tempPath)
        TXTFile=open(tempPath, 'a')
        fh = open(self.lineEdit_vtk.currentPath, 'rU')
        for vtkfilename in fh.readlines():
            vtkfilename = os.path.dirname(self.lineEdit_vtk.currentPath) + '/' + vtkfilename.rstrip()
            TXTFile.write(vtkfilename)
            TXTFile.write("\n")
        fh.close()        
        TXTFile.close()   
        args=arguments(coordData=self.lineEdit_template.currentPath, covariate=self.lineEdit_output.directory+'/valuesCsv.csv', covariateType=self.lineEdit_output.directory+'/CovariateNames.csv', outputDir=self.lineEdit_output.directory, shapeData=tempPath, shapePath=PathOutput)
        # run_script(args)
        pvaluesPath=self.lineEdit_output.directory+'/pvalues.json'
        efitPath=self.lineEdit_output.directory+'/efit.json'
        ShapePvalue=self.lineEdit_ShapePvalue.currentPath
        outputPath=self.lineEdit_output.directory+'/out.vtk'
        covariates = covariates.tolist()
        print("covariates dans MFSDA", covariates)
        argShapes=arguments(pvalues=pvaluesPath, covariates=covariates, shape=ShapePvalue, efit=efitPath, output=outputPath)
        run_script2(argShapes)
        CSVFile=open(self.lineEdit_output.directory+'/output.csv', 'wb')
        Data=['VTK Files']
        with CSVFile:
            writer = csv.writer(CSVFile)
            writer.writerow(Data)
            writer.writerow([outputPath])
        CSVFile.close()
        parameters = {}
        parameters["CSVFile"] = self.lineEdit_output.directory+'/output.csv'
        module = slicer.modules.shapepopulationviewer
        slicer.cli.run(module, None, parameters, wait_for_completion=True)

#
# MFSDALogic
#

class MFSDALogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
        computation done by your module.  The interface
        should be such that other python code can import
        this class and make use of the functionality without
        requiring an instance of the Widget.
        Uses ScriptedLoadableModuleLogic base class, available at:
        https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
        """
    def __init__(self, interface):
        
        self.interface = interface
        self.table = vtk.vtkTable
        self.colorBar = {'Point1': [0, 0, 1, 0], 'Point2': [0.5, 1, 1, 0], 'Point3': [1, 1, 0, 0]}
#        self.input_Data = inputData.inputData()

    
    def get(self, objectName):
        
        """ Functions to recovery the widget in the .ui file
            """
        return slicer.util.findChild(self.interface.widget, objectName)
    
    def checkExtension(self, filename, extension):
        """ Check if the path given has the right extension"""
        if os.path.splitext(os.path.basename(filename))[1] == extension :
            return True
        elif os.path.basename(filename) == "" or os.path.basename(filename) == " " :
            return False
        slicer.util.errorDisplay('Wrong extension file, a ' + extension + ' file is needed!' +filename)
        return False
    
    def onAddGroupForCreationCSVFile(self):
        """Function to add a group of the dictionary
            - Add the paths of all the vtk files found in the directory given
            of a dictionary which will be used to create the CSV file"""
        directory = self.directoryButton_creationCSVFile.directory.encode('utf-8')
        if directory in self.directoryList:
            index = self.directoryList.index(directory) + 1
            slicer.util.errorDisplay('Path of directory already used for the group ' + str(index))
            return
        
        # Add the paths of vtk files of the dictionary
        self.logic.addGroupToDictionary(self.dictCSVFile, directory, self.directoryList, self.spinBox_group.value)
        condition = self.logic.checkSeveralMeshInDict(self.dictCSVFile)
        
        if not condition:
            # Remove the paths of vtk files of the dictionary
            self.logic.removeGroupToDictionary(self.dictCSVFile, self.directoryList, self.spinBox_group.value)
            return
        
        # Increment of the number of the group in the spinbox
        self.spinBox_group.blockSignals(True)
        self.spinBox_group.setMaximum(self.spinBox_group.value + 1)
        self.spinBox_group.setValue(self.spinBox_group.value + 1)
        self.spinBox_group.blockSignals(False)
        
        # Message for the user
    slicer.util.delayDisplay("Group Added")


class arguments:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class MFSDATest(ScriptedLoadableModuleTest):
    
    
    # def runTest(self):
    #     cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"Testing")))
    #     print("in the Test")
    #     vtk_currentPath=cmd_subfolder+'/shapeData.txt'
    #     template_currentPath=cmd_subfolder+'/ALLM_aligned.vtk'
    #     output_directory=cmd_subfolder+'/output'
    #     try:
    #         os.stat(output_directory)
    #     except:
    #         os.mkdir(output_directory)
    #     covariate_currentPath=cmd_subfolder+'/TEST/covariate.csv'
    #     covariateType_currentPath=cmd_subfolder+'/TEST/covariateType.csv'
    #     ShapePvalue_currentPath=cmd_subfolder+'/01_Left_aligned.vtk'
    #     CovariateNames_currentPath=cmd_subfolder+'/TEST/covariate'
    #     PathOutput=os.path.dirname(vtk_currentPath)+'/'
    #     fshape = open(vtk_currentPath, 'rU')
    #     copy=PathOutput+'/temp.txt'
    #     ftemp= open(copy,'w')
    #     for vtkfilename in fshape.readlines():
    #         ftemp.write(PathOutput+vtkfilename)
    #     fshape.close()
    #     ftemp.close()
    #     args=arguments(coordData=template_currentPath, covariate=covariate_currentPath, covariateType=covariateType_currentPath, outputDir=output_directory, shapeData=copy, shapePath=PathOutput)
    #     run_script(args)
    #     pvaluesPath=output_directory+'/pvalues.json'
    #     efitPath=output_directory+'/efit.json'
    #     ShapePvalue=ShapePvalue_currentPath
    #     outputPath=output_directory+'/out.vtk'
    #     argShapes=arguments(pvalues=pvaluesPath, covariates=CovariateNames_currentPath, shape=ShapePvalue_currentPath, efit=efitPath, output=outputPath)
    #     run_script2(argShapes)
    #     CSVFile=open(output_directory+'/output.csv', 'wb')
    #     Data=['VTK Files']
    #     with CSVFile:
    #         writer = csv.writer(CSVFile)
    #         writer.writerow(Data)
    #         writer.writerow([outputPath])
    #     CSVFile.close()




###################################################################################################
###################################################################################################
        ###########################################################################################
    #        parameters = {}
    #        parameters["CSVFile"] = output_directory+'/output.csv'
    #        module = slicer.modules.shapepopulationviewer
    #        slicer.cli.run(module, None, parameters, wait_for_completion=True)
    
    
    
    
    
    """
        This is the test case for your scripted module.
        Uses ScriptedLoadableModuleTest base class, available at:
        https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
        """


