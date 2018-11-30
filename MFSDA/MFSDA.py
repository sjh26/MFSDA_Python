import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
#import inputData
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),os.path.join('Resources','Libraries')))
import numpy as np
from scipy import stats
#from scipy.io import loadmat
#from statsmodels.sandbox.stats.multicomp import fdrcorrection0
#from stat_read_x import read_x
#from stat_lpks import lpks
#from stat_sif import sif
#from stat_wald_ht import wald_ht
#from stat_bstrp_pvalue import bstrp_pvalue
# from sklearn.cluster import KMeans
# from scipy.cluster.vq import kmeans2
# from stat_gap import gap
import MFSDA_stat as mfsda
import timeit
import argparse
import json
import os.path
import csv 
import time

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
		# test=False
		self.logic = MFSDALogic(self)
		# self.test = MFSDATest(self)
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
		self.lineEdit_covariate = self.logic.get('lineEdit_covariate')
		#self.lineEdit_covariateType = self.logic.get('lineEdit_covariateType')
		self.lineEdit_output = self.logic.get('lineEdit_output')

		self.pushButton_run = self.logic.get('pushButton_run')
		# self.apply_button=self.logic.get('Apply_Button')
		self.lineEdit_CovariateNames = self.logic.get('lineEdit_CovariateNames')
		self.lineEdit_ShapePvalue=self.logic.get('lineEdit_ShapePvalue')

		'''
		self.lineEdit_covariateType.connect('currentPathChanged(const QString)', self.onCSVFile)
		self.lineEdit_covariate.connect('currentPathChanged(const QString)', self.onCSVFile)
		self.lineEdit_template.connect('currentPathChanged(const QString)', self.onCSVFile)
		self.lineEdit_vtk.connect('currentPathChanged(const QString)', self.onCSVFile)
		self.lineEdit_output.connect('directoryChanged (const QString)', self.onCSVFile)
		self.lineEdit_ShapePvalue.connect('directoryChanged (const QString)', self.onCSVFile)
		self.lineEdit_CovariateNames.connect('directoryChanged (const QString)', self.onCSVFile)'''		
		self.pushButton_run.connect('clicked(bool)', self.onCSVFile)


	def onShapeState(self):
	    state=self.MFSDAShapeThread.GetStatusString()
	    if state=='Running'or state=='Scheduled':
	        seconds = time.time()-self.starting_time
	        m, s = divmod(seconds, 60)
	        h, m = divmod(m, 60)
	        if h==0 and m==0:
	            t = "00:%02d" % (s)
	        elif h==0 :
	            t = "%02d:%02d" % (m, s)
	        else:
	            t = "%d:%02d:%02d" % (h, m, s)
	        if int(s) ==0:
	            print("MFSDA shape creation "+self.MFSDAShapeThread.GetStatusString()+"  "+t)

	        self.pushButton_run.setText("Abort shape creation ("+t+")")
	    else:
	        print('MFSDA shape creation done')
	        self.checkThreadTimer.stop()
	        self.checkThreadTimer.disconnect('timeout()', self.onShapeState)

	        CSVFile=open(self.lineEdit_output.directory+'/output.csv', 'wb')
	        Data=['VTK Files']
	        with CSVFile:
	        	writer = csv.writer(CSVFile)
	        	writer.writerow(Data)
	        	writer.writerow([self.lineEdit_ShapePvalue.currentPath])
	        CSVFile.close()
	        self.pushButton_run.setText("Run")
	        self.pushButton_run.connect('clicked(bool)', self.onCSVFile)
	        self.pushButton_run.disconnect('clicked(bool)', self.onKillComputationShape)
	        parameters = {}
	        parameters["CSVFile"] = self.lineEdit_output.directory+'/output.csv'
	        module = slicer.modules.shapepopulationviewer
	        slicer.cli.run(module, None, parameters, wait_for_completion=True)


	    return

	def onKillComputationShape(self):
        
	    self.pushButton_run.clicked.disconnect()
	    self.pushButton_run.connect('clicked()',self.onCSVFile)
	    self.pushButton_run.setText("Run")

	    self.checkThreadTimer.stop()
	    self.checkThreadTimer.timeout.disconnect()

	    self.MFSDAShapeThread.Cancel()
	    minutes = int((time.time()-self.starting_time)/60)
	    print('Computation Stopped after '+str(minutes)+' min')


	def onComputationState(self):
	    state=self.MFSDAThread.GetStatusString()
	    if state=='Running'or state=='Scheduled':
	        seconds = time.time()-self.starting_time
	        m, s = divmod(seconds, 60)
	        h, m = divmod(m, 60)
	        if h==0 and m==0:
	            t = "00:%02d" % (s)
	        elif h==0 :
	            t = "%02d:%02d" % (m, s)
	        else:
	            t = "%d:%02d:%02d" % (h, m, s)
	        if int(s) ==0:
	            print("MFSDA computation "+self.MFSDAThread.GetStatusString()+"  "+t)

	        self.pushButton_run.setText("Abort computation ("+t+")")
	    else:
	        print('MFSDA Computation done')
	        self.checkThreadTimer.stop()

	        self.param = {}
	        self.param["shape"] = self.lineEdit_ShapePvalue.currentPath
	        self.param["pvalues"] = self.lineEdit_output.directory+'/pvalues.json'
	        self.param["efit"] = self.lineEdit_output.directory+'/efit.json'
	        self.param["covariates"] = self.lineEdit_CovariateNames.currentPath
	        self.param["output"] = self.lineEdit_output.directory+'/out.vtk'

	        MFSDAShapemodule = slicer.modules.mfsda_createshapes
	        self.MFSDAShapeThread=slicer.cli.run(MFSDAShapemodule, None, self.param, wait_for_completion=False)

		    
	        self.checkThreadTimer.disconnect('timeout()', self.onComputationState)
	        self.checkThreadTimer.connect('timeout()', self.onShapeState)
	        self.checkThreadTimer.start(1000)

	        self.pushButton_run.connect('clicked(bool)', self.onKillComputationShape)
	        self.pushButton_run.disconnect('clicked(bool)', self.onKillComputation)


	def onKillComputation(self):
	    
	    self.pushButton_run.clicked.disconnect()
	    self.pushButton_run.connect('clicked()',self.onCSVFile)
	    self.pushButton_run.setText("Run")

	    self.checkThreadTimer.stop()
	    self.checkThreadTimer.timeout.disconnect()

	    self.MFSDAThread.Cancel()
	    minutes = int((time.time()-self.starting_time)/60)
	    print('Computation Stopped after '+str(minutes)+' min')

	def onCSVFile(self):
		print('coucou')
		self.dictShapeModels = dict()
		'''if not os.path.exists(self.lineEdit_covariateType.currentPath):
			self.stateCSVMeansShape = False
			return'''
		if not os.path.exists(self.lineEdit_covariate.currentPath):
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
		"""if not os.path.exists(self.lineEdit_CovariateNames.currentPath):
			self.stateCSVMeansShape = False
			return"""
		condition1 = self.logic.checkExtension(self.lineEdit_vtk.currentPath, ".txt")
		condition2 = self.logic.checkExtension(self.lineEdit_template.currentPath, ".vtk")
		condition3 = self.logic.checkExtension(self.lineEdit_covariate.currentPath, ".csv")
		#condition4 = self.logic.checkExtension(self.lineEdit_covariateType.currentPath, ".csv")
		condition5= self.lineEdit_output.directory!='.'
		condition6 = self.logic.checkExtension(self.lineEdit_ShapePvalue.currentPath, ".vtk")
		condition7= self.lineEdit_CovariateNames.currentPath!='.'

		if not condition1:
			self.lineEdit_vtk.setCurrentPath(" ")
			self.stateCSVDataset = False
			return
		if not condition2:
			self.lineEdit_template.setCurrentPath(" ")
			self.stateCSVDataset = False
			return
		if not condition3:
			self.lineEdit_covariate.setCurrentPath(" ")
			self.stateCSVDataset = False
			return
		'''if not condition4:
			self.lineEdit_covariateType.setCurrentPath(" ")
			self.stateCSVDataset = False
			return'''
		if not condition5:
			#self.lineEdit_output.setDirectory(" ")
			self.stateCSVDataset = False
			return
		if not condition6:
			self.lineEdit_ShapePvalue.setCurrentPath(" ")
			self.stateCSVDataset = False
			return
		'''if not condition7:
			self.lineEdit_CovariateNames.setCurrentPath(" ")
			self.stateCSVDataset = False
			return'''
		PathOutput=os.path.dirname(self.lineEdit_vtk.currentPath)+'/'

		# print('bouton !!!!',self.apply_button.clicked())
		#print(self.lineEdit_covariateType.currentPath)
		print(self.lineEdit_covariate.currentPath)
		print(self.lineEdit_template.currentPath)
		print(self.lineEdit_vtk.currentPath)
		print(self.lineEdit_output.directory)
		print(self.lineEdit_ShapePvalue.currentPath)
		print(self.lineEdit_CovariateNames.currentPath)		
		print(PathOutput)



		self.param = {}
		self.param["shapeData"] = self.lineEdit_vtk.currentPath
		self.param["coordData"] = self.lineEdit_template.currentPath
		self.param["covariate"] = self.lineEdit_covariate.currentPath
		#self.param["covariateType"] = self.lineEdit_covariateType.currentPath
		self.param["outputDir"] = self.lineEdit_output.directory


		self.starting_time=time.time()
		MFSDAmodule = slicer.modules.mfsda_run
		self.MFSDAThread=slicer.cli.run(MFSDAmodule, None, self.param, wait_for_completion=False)


		self.checkThreadTimer=qt.QTimer()
		self.checkThreadTimer.connect('timeout()', self.onComputationState)
		self.checkThreadTimer.start(1000)

		self.pushButton_run.disconnect('clicked(bool)', self.onCSVFile)
		self.pushButton_run.connect('clicked(bool)', self.onKillComputation)

		return
 
        '''		self.param = {}
		self.param["shape"] = self.lineEdit_ShapePvalue.currentPath
		self.param["pvalues"] = self.lineEdit_output.directory+'/pvalues.json'
		self.param["efit"] = self.lineEdit_output.directory+'/efit.json'
		self.param["covariates"] = self.lineEdit_CovariateNames.currentPath
		self.param["output"] = self.lineEdit_output.directory+'/out.vtk'

		MFSDAShapemodule = slicer.modules.MFSDA_createShapes
		self.MFSDAShapeThread=slicer.cli.run(MFSDAShapemodule, None, self.param, wait_for_completion=False)

		return

		CSVFile=open(self.lineEdit_output.directory+'/output.csv', 'wb')
		Data=['VTK Files']
		with CSVFile:
			writer = csv.writer(CSVFile)
			writer.writerow(Data)
			writer.writerow([self.lineEdit_ShapePvalue.currentPath])
		CSVFile.close()
		parameters = {}
		parameters["CSVFile"] = self.lineEdit_output.directory+'/output.csv'
		module = slicer.modules.shapepopulationviewer
		slicer.cli.run(module, None, parameters, wait_for_completion=True)

		return


		args=arguments(coordData=self.lineEdit_template.currentPath, covariate=self.lineEdit_covariate.currentPath, covariateType=self.lineEdit_covariateType.currentPath, outputDir=self.lineEdit_output.directory, shapeData=self.lineEdit_vtk.currentPath, shapePath=PathOutput)
		self.logic.run_script(args)


		pvaluesPath=self.lineEdit_output.directory+'/pvalues.json'
		efitPath=self.lineEdit_output.directory+'/efit.json'
		ShapePvalue=self.lineEdit_ShapePvalue.currentPath
		outputPath=self.lineEdit_output.directory+'/out.vtk'
		argShapes=arguments(pvalues=pvaluesPath, covariates=self.lineEdit_CovariateNames.currentPath, shape=self.lineEdit_ShapePvalue.currentPath, efit=efitPath, output=outputPath)
		self.logic.run_Shape(argShapes)


		CSVFile=open(self.lineEdit_output.directory+'/output.csv', 'wb')
		Data=['VTK Files']
		with CSVFile:
			writer = csv.writer(CSVFile)
			writer.writerow(Data)
			writer.writerow([self.lineEdit_ShapePvalue.currentPath])
		CSVFile.close()
		parameters = {}
		parameters["CSVFile"] = self.lineEdit_output.directory+'/output.csv'
		module = slicer.modules.shapepopulationviewer
		slicer.cli.run(module, None, parameters, wait_for_completion=True)'''
			



		# self.lineEdit_4.setCurrentPath(" ")
		# Parameters Area
		#
	# 	self.dictShapeModels = dict()

 #        # Check if the path exists:
 #        if not os.path.exists(self.lineEdit_4.currentPath):
 #            self.stateCSVDataset = False
 #            self.enableNetwork()
 #            # return

 #        # print("------ Selection of a Dataset ------")
 #        # Check if it's a CSV file
 #        condition1 = self.logic.checkExtension(self.lineEdit_4.currentPath, ".csv")
 #        if not condition1:
 #            self.lineEdit_4.setCurrentPath(" ")
 #            self.stateCSVDataset = False
 #            self.enableNetwork()
 #            # return
	
	

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
		#self.input_Data = inputData.inputData()
			
			
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

	def run_script(self, args):
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
			vtkfilename=args.shapePath+vtkfilename
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
		#print('y_design', y_design, 'y_design', coord_mat, 'design_data', design_data, 'var_type', var_type)
		gpvals, lpvals_fdr, clu_pvals, efit_beta, efity_design, efit_eta = mfsda.run_stats(y_design, coord_mat, design_data, var_type)
			
		"""+++++++++++++++++++++++++++++++++++"""
		"""Step3. Save all the results"""
			
		if not os.path.exists(args.outputDir):
			os.makedirs(args.outputDir)

		pvalues = {}
		pvalues['Gpvals'] = gpvals.tolist()
		pvalues['clu_pvals'] = clu_pvals.tolist()
		pvalues['Lpvals_fdr'] = lpvals_fdr.tolist()

		with open(os.path.join(args.outputDir,'pvalues.json'), 'w') as outfile:
			json.dump(pvalues, outfile)
			
		efit = {}
		efit['efitBetas'] = efit_beta.tolist()
		efit['efitYdesign'] = efity_design.tolist()
		efit['efitEtas'] = efit_eta.tolist()
		
		with open(os.path.join(args.outputDir,'efit.json'), 'w') as outfile:
			json.dump(efit, outfile)

		# if __name__ == '__main__':
			
	def run_Shape(self, args):

		with open(args.pvalues) as data_file:
			pvalues = json.load(data_file)
			lpvalues = np.array(pvalues['Lpvals_fdr'])

		with open(args.efit) as data_file:
			efit = json.load(data_file)
			efitbetas = np.array(efit['efitBetas'])

		covariates = args.covariates

		reader = vtk.vtkPolyDataReader()
		reader.SetFileName(args.shape)
		reader.Update()
		shapedata = reader.GetOutput()

		pointdata = shapedata.GetPointData()	

		print("Adding pvalues...", lpvalues.shape)

		for j in range(lpvalues.shape[1]):

			arr = vtk.vtkDoubleArray()
			name = 'pvalue_'
			if covariates and j < len(covariates):
				name += covariates[j]
			else:
				name += str(j)

			arr.SetName(name)

			for i in range(lpvalues.shape[0]):			
				arr.InsertNextTuple([lpvalues[i][j]])
			pointdata.AddArray(arr)


		print("Adding betas...", efitbetas.shape)

		for k in range(efitbetas.shape[2]):
			for i in range(efitbetas.shape[0]):		
				arr = vtk.vtkDoubleArray()

				name = 'betavalues_' + str(k) + "_"

				if covariates and i < len(covariates):
					name += covariates[i]
				else:
					name += str(i)

				arr.SetName(name)
				
				for j in range(efitbetas.shape[1]):				
					arr.InsertNextTuple([efitbetas[i][j][k]])
				pointdata.AddArray(arr)


		print("Writing output vtk file...")

		writer = vtk.vtkPolyDataWriter()
		writer.SetFileName(args.output)
		writer.SetInputData(shapedata)
		writer.Update()
		# 	args = parser.parse_args()
			
		# 	start_all = timeit.default_timer()
		# 	run_script(args)
		# 	stop_all = timeit.default_timer()
		# 	delta_time_all = str(stop_all - start_all)
		# 	print("The total elapsed time is " + delta_time_all)

class arguments:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)				
	
class MFSDATest(ScriptedLoadableModuleTest):

	def SetUp(self):
		slicer.mrmlScene.Clear(0)

	def runTest(self):
		self.SetUp()
		vtk_currentPath='/Users/loicmichoud/Desktop/ficher_test/shapeData.txt'
		template_currentPath='/Users/loicmichoud/Desktop/ficher_test/ALLM_aligned.vtk'
		output_directory='/Users/loicmichoud/Desktop/MySlicerExtensions/output'
		covariate_currentPath='/Users/loicmichoud/Desktop/MFSDA_Python/TEST/covariate.csv'
		covariateType_currentPath='/Users/loicmichoud/Desktop/MFSDA_Python/TEST/covariateType.csv'
		ShapePvalue_currentPath='/Users/loicmichoud/Desktop/ficher_test/01_Left_aligned.vtk'
		CovariateNames_currentPath='/Users/loicmichoud/Desktop/MySlicerExtensions/MFSDA/MFSDA/covariate'
		PathOutput=os.path.dirname(vtk_currentPath)+'/'
		args=arguments(coordData=template_currentPath, covariate=covariate_currentPath, covariateType=covariateType_currentPath, outputDir=output_directory, shapeData=vtk_currentPath, shapePath=PathOutput)
		self.run_script(args)
		pvaluesPath=output_directory+'/pvalues.json'
		efitPath=output_directory+'/efit.json'
		ShapePvalue=ShapePvalue_currentPath
		outputPath=output_directory+'/out.vtk'
		argShapes=arguments(pvalues=pvaluesPath, covariates=CovariateNames_currentPath, shape=ShapePvalue_currentPath, efit=efitPath, output=outputPath)
		self.run_Shape(argShapes)
		CSVFile=open(output_directory+'/output.csv', 'wb')
		Data=['VTK Files']
		with CSVFile:
			writer = csv.writer(CSVFile)
			writer.writerow(Data)
			writer.writerow([ShapePvalue_currentPath])
		CSVFile.close()
		parameters = {}
		parameters["CSVFile"] = output_directory+'/output.csv'
		module = slicer.modules.shapepopulationviewer
		slicer.cli.run(module, None, parameters, wait_for_completion=True)

	def run_script(self, args):
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
			vtkfilename=args.shapePath+vtkfilename
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
#        print('y_design', y_design, 'y_design', coord_mat, 'design_data', design_data, 'var_type', var_type)
		gpvals, lpvals_fdr, clu_pvals, efit_beta, efity_design, efit_eta = mfsda.run_stats(y_design, coord_mat, design_data, var_type)
			
		"""+++++++++++++++++++++++++++++++++++"""
		"""Step3. Save all the results"""
			
		if not os.path.exists(args.outputDir):
			os.makedirs(args.outputDir)

		pvalues = {}
		pvalues['Gpvals'] = gpvals.tolist()
		pvalues['clu_pvals'] = clu_pvals.tolist()
		pvalues['Lpvals_fdr'] = lpvals_fdr.tolist()

		with open(os.path.join(args.outputDir,'pvalues.json'), 'w') as outfile:
			json.dump(pvalues, outfile)
			
		efit = {}
		efit['efitBetas'] = efit_beta.tolist()
		efit['efitYdesign'] = efity_design.tolist()
		efit['efitEtas'] = efit_eta.tolist()
		
		with open(os.path.join(args.outputDir,'efit.json'), 'w') as outfile:
			json.dump(efit, outfile)

		# if __name__ == '__main__':
			
	def run_Shape(self, args):

		with open(args.pvalues) as data_file:
			pvalues = json.load(data_file)
			lpvalues = np.array(pvalues['Lpvals_fdr'])

		with open(args.efit) as data_file:
			efit = json.load(data_file)
			efitbetas = np.array(efit['efitBetas'])

		covariates = args.covariates

		reader = vtk.vtkPolyDataReader()
		reader.SetFileName(args.shape)
		reader.Update()
		shapedata = reader.GetOutput()

		pointdata = shapedata.GetPointData()	

		print("Adding pvalues...", lpvalues.shape)

		for j in range(lpvalues.shape[1]):

			arr = vtk.vtkDoubleArray()
			name = 'pvalue_'
			if covariates and j < len(covariates):
				name += covariates[j]
			else:
				name += str(j)

			arr.SetName(name)

			for i in range(lpvalues.shape[0]):			
				arr.InsertNextTuple([lpvalues[i][j]])
			pointdata.AddArray(arr)


		print("Adding betas...", efitbetas.shape)

		for k in range(efitbetas.shape[2]):
			for i in range(efitbetas.shape[0]):		
				arr = vtk.vtkDoubleArray()

				name = 'betavalues_' + str(k) + "_"

				if covariates and i < len(covariates):
					name += covariates[i]
				else:
					name += str(i)

				arr.SetName(name)
				
				for j in range(efitbetas.shape[1]):				
					arr.InsertNextTuple([efitbetas[i][j][k]])
				pointdata.AddArray(arr)


		print("Writing output vtk file...")

		writer = vtk.vtkPolyDataWriter()
		writer.SetFileName(args.output)
		writer.SetInputData(shapedata)
		writer.Update()


	"""
	This is the test case for your scripted module.
	Uses ScriptedLoadableModuleTest base class, available at:
	https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
	"""
