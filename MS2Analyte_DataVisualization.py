import sys
import os
from PyQt5 import QtCore, QtWidgets
from qtmodern.styles import dark
from qtmodern.windows import ModernWindow
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import pyqtgraph as pg
import numpy as np
import pandas
import random
import time
import pickle
from pyqtgraph.Qt import QtCore, QtGui
from decimal import Decimal, ROUND_DOWN
import cProfile
import pstats
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QLineEdit, QTabWidget, \
				  QGridLayout, QVBoxLayout, QHBoxLayout, QGroupBox, QDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from ms2analyte.file_handling import data_import




# Handle high resolution displays:
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
########### Data Import #################

def fetch_sample_list(input_structure, input_type):
	sample_list = data_import.name_extract(input_structure, input_type)
	
	return sample_list



#######################################################################################
def import_dataframe(input_structure, input_type, sample_name):
	with open((os.path.join(input_structure.output_directory, input_type, sample_name + "_all_replicates_blanked_dataframe.pickle")), 'rb') as f:
		df = pickle.load(f)
	f.close()

	return df



#######################################################################################
def import_experiment_dataframe(input_structure):
	with open((os.path.join(input_structure.output_directory,  input_structure.experiment_name + "_experiment_analyte_overview_tableau_output.csv"))) as f:
		df = pandas.read_csv(f)

	return df
    




#######################################################################################
def import_ms1_dataframe(input_structure, input_type, sample_name):
	with open((os.path.join(input_structure.output_directory, input_type, sample_name + '_replicate_analyte_mass_spectra.pickle')), 'rb') as f:


		inputdata = pickle.load(f)
		

	b = 0
	replicate_analyte_id = []
	average_mass = []
	relative_intensity = []

	while b < len(inputdata):
		i = 0
		while i <len(inputdata[b].replicate_analyte_spectrum):
			replicate_analyte_id.append(inputdata[b].replicate_analyte_id)
			average_mass.append(inputdata[b].replicate_analyte_spectrum[i].average_mass)
			relative_intensity.append(inputdata[b].replicate_analyte_spectrum[i].relative_intensity)
			i+=1


		b+=1
	np.round(average_mass,1)
	relative_intensity_zeros = [0]*len(relative_intensity)
	average_mass_lower = [0]*len(average_mass)
	i=0
	while i < len(average_mass):
		average_mass_lower[i] = average_mass[i] - 0.0001
		i+=1

	average_mass_upper = [0]*len(average_mass)
	i=0
	while i < len(average_mass):
		average_mass_upper[i] = average_mass[i] + 0.0001
		i+=1

	data1 = {'replicate_analyte_id': replicate_analyte_id,
			'average_mass': average_mass,
			'relative_intensity': relative_intensity
			}
	data2 = {'replicate_analyte_id': replicate_analyte_id,
			'average_mass': average_mass_upper,
			'relative_intensity': relative_intensity_zeros
			}

	data3 = {'replicate_analyte_id': replicate_analyte_id,
			'average_mass': average_mass_lower,
			'relative_intensity': relative_intensity_zeros
			}

	df1 = pandas.DataFrame (data1, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

	df2 = pandas.DataFrame (data2, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])

	df3 = pandas.DataFrame (data3, columns = ['replicate_analyte_id', 'average_mass','relative_intensity'])


	df_combine = [df1,df2,df3]

	df_combine = pandas.concat(df_combine)




	with open((os.path.join(input_structure.output_directory, input_type, sample_name + '_R1_analytes.pickle')), 'rb') as f:
		inputdata2 = pickle.load(f)

	b = 0
	analyte_id = []
	max_peak_intensity_mass = []
	max_peak_intensity = []


	while b < len(inputdata2):

		analyte_id.append(inputdata2[b].analyte_id)
		max_peak_intensity_mass.append(inputdata2[b].max_peak_intensity_mass)
		max_peak_intensity.append(inputdata2[b].max_peak_intensity)



		b+=1


	max_peak_intensity_zeros = [0]*len(max_peak_intensity)
	max_peak_intensity_mass_lower = [0]*len(max_peak_intensity_mass)
	i=0
	while i < len(max_peak_intensity_mass):
		max_peak_intensity_mass_lower[i] = max_peak_intensity_mass[i] - 0.0001
		i+=1





	max_peak_intensity_zeros = [0]*len(max_peak_intensity)
	max_peak_intensity_mass_upper = [0]*len(max_peak_intensity_mass)
	i=0
	while i < len(max_peak_intensity_mass):
		max_peak_intensity_mass_upper[i] = float(max_peak_intensity_mass[i] + 0.0001)
		i+=1


	data4 = {'analyte_id': analyte_id,
			'max_peak_intensity_mass': max_peak_intensity_mass,
			'max_peak_intensity': max_peak_intensity
			}
	data5 = {'analyte_id': analyte_id,
			'max_peak_intensity_mass': max_peak_intensity_mass_upper,
			'relative_intensity': max_peak_intensity_zeros
			}

	data6 = {'analyte_id': analyte_id,
			'max_peak_intensity_mass': max_peak_intensity_mass_lower,
			'max_peak_intensity': max_peak_intensity_zeros
			}

	df4 = pandas.DataFrame (data4, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

	df5 = pandas.DataFrame (data5, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

	df6 = pandas.DataFrame (data6, columns = ['analyte_id', 'max_peak_intensity_mass','max_peak_intensity'])

	df_combine2 = [df4,df5,df6]

	df_combine2 = pandas.concat(df_combine2)








	with open((os.path.join(input_structure.output_directory,'Samples',  input_structure.experiment_name + '_experiment_analyte_mass_spectra.pickle')), 'rb') as f:




		inputdata3 = pickle.load(f)
		

	b = 0
	experiment_analyte_id = []
	average_mass = []
	relative_intensity = []

	while b < len(inputdata3):
		i = 0
		while i <len(inputdata3[b].relative_experiment_mass_spectrum):
			experiment_analyte_id.append(inputdata3[b].experiment_analyte_id)
			average_mass.append(inputdata3[b].relative_experiment_mass_spectrum[i].average_mass)
			relative_intensity.append(inputdata3[b].relative_experiment_mass_spectrum[i].relative_intensity)
			i+=1


		b+=1
	np.round(average_mass,1)
	relative_intensity_zeros = [0]*len(relative_intensity)
	average_mass_lower = [0]*len(average_mass)
	i=0
	while i < len(average_mass):
		average_mass_lower[i] = average_mass[i] - 0.0001
		i+=1

	average_mass_upper = [0]*len(average_mass)
	i=0
	while i < len(average_mass):
		average_mass_upper[i] = average_mass[i] + 0.0001
		i+=1

	data7 = {'experiment_analyte_id': experiment_analyte_id,
			'average_mass': average_mass,
			'relative_intensity': relative_intensity
			}
	data8 = {'experiment_analyte_id': experiment_analyte_id,
			'average_mass': average_mass_upper,
			'relative_intensity': relative_intensity_zeros
			}

	data9 = {'experiment_analyte_id': experiment_analyte_id,
			'average_mass': average_mass_lower,
			'relative_intensity': relative_intensity_zeros
			}

	df7 = pandas.DataFrame (data7, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])

	df8 = pandas.DataFrame (data8, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])

	df9 = pandas.DataFrame (data9, columns = ['experiment_analyte_id', 'average_mass','relative_intensity'])


	df_combine3 = [df7,df8,df9]

	df_combine3 = pandas.concat(df_combine3)


	return df_combine, df_combine2, df_combine3

	
#################################################################################################
#### Fetch Input Structure, Input type and Sample name list

input_structure = data_import.input_data_structure()
input_type = "Samples"
sample_names = fetch_sample_list(input_structure, input_type)


#################################### Start GUI  ##################################################

class State:
	def __init__(self, null, blank, massMin, massMax, rtMin, rtMax):

		self.null = null
		self.blank = blank
		self.massMin = massMin
		self.massMax = massMax
		self.rtMin = rtMin
		self.rtMax = rtMax



	def return_null_state(self):

		print("State:", self.null)
		return self.null


	def return_blank_state(self):

		print("State:", self.blank)
		return self.blank


	def return_massMin_state(self):

		print("State:", self.massMin)
		return self.massMin



	def return_massMax_state(self):

		print("State:", self.massMax)
		return self.massMax


	def return_rtMin_state(self):

		print("State:", self.rtMin)
		return self.rtMin


	def return_rtMax_state(self):

		print("State:", self.rtMax)
		return self.rtMax


class Tab_State:
	def __init__(self, sample_state, replicate_state, experiment_state, diversity_state):

		self.sample_state = sample_state
		self.replicate_state = replicate_state
		self.experiment_state = experiment_state
		self.diversity_state = diversity_state

	def return_sample_state(self):

		print("Sample State:", self.sample_state)
		return self.sample_state

	def return_replicate_state(self):

		print("Replicate State:", self.replicate_state)
		return self.replicate_state

	def return_experiment_state(self):

		print("Experiment State:", self.experiment_state)
		return self.experiment_state

	def return_diversity_state(self):

		print("Diversity State:", self.diversity_state)

		return self.diversity_state




class Window (QDialog):



	def __init__(self):
		super(Window, self).__init__()


        # changing the background color to yellow 

		
		app = QtWidgets.QApplication(sys.argv)
		#dark(app) 
		app.setStyle('Fusion')
		# self.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")

		self.setWindowTitle("MS2Analyte")
		self.setGeometry(100,100,1700,900)
		self.setWindowIcon(QIcon('ms2analyte_icon.png'))
		pg.setConfigOption('background', 'w')
		pg.setConfigOption('foreground', 'k')


		self.region_1()
		self.Plot_Sample()
		self.Plot_ms1()
		windowLayout = QVBoxLayout()
		# windowLayout.setStyleSheet("background-color: white;")
		windowLayout.addWidget(self.horizontalGroupBox)
		self.setLayout(windowLayout)
		# extractAction = QAction("&Quit", self)
		# extractAction.setShortcut("Ctrl+Q")
		# extractAction.setStatusTip('Leave The App')
		# extractAction.triggered.connect(self.close_application)
		
		# self.statusBar()
		
		# mainMenu = self.menuBar()
		# mainMenu.resize(100,100)
		
		# fileMenu = mainMenu.addMenu('&File')
		# fileMenu.addAction(extractAction)
		# toolsMenu = mainMenu.addMenu('&Tools')


		 




		# self.setStyleSheet("background-color: white;")


		
######## Hide some plots initially ########
		
		self.mz_r1_plot.hide()
		self.mz_r2_plot.hide()
		self.mz_r3_plot.hide()

		self.rt_r1_plot.hide()
		self.rt_r2_plot.hide()
		self.rt_r3_plot.hide()
		

		self.ex_bOn_plot.hide()
		self.div_plot.hide()


		# self.showMaximized()
		self.show()

	def region_1(self):
		self.horizontalGroupBox = QGroupBox()
		layoutMain = QGridLayout()
		

		self.btn2 = QPushButton("Sample", self)
		self.btn2.clicked.connect(self.Sample_Plot_DF)
		self.btn2.clicked.connect(self.show_sample)
		self.btn2.clicked.connect(self.Plot_Sample)
		self.btn2.clicked.connect(self.hide_replicate)
		self.btn2.clicked.connect(self.hide_experiment)
		self.btn2.clicked.connect(self.hide_diversity)

		# self.btn2.resize(buttonX,buttonY)
		# self.btn2.move(20 + buttonmoveX,100 + buttonmoveY)
		self.btn2.setDisabled(True)

		
		self.btn3 = QPushButton("Replicate", self)
		self.btn3.clicked.connect(self.Replicate_Plot_DF)
		self.btn3.clicked.connect(self.show_replicate)
		self.btn3.clicked.connect(self.Plot_Replicate)
		self.btn3.clicked.connect(self.hide_sample)

		self.btn3.clicked.connect(self.hide_experiment)
		self.btn3.clicked.connect(self.hide_diversity)


		# self.btn3.resize(buttonX,buttonY)
		# self.btn3.move(170 + buttonmoveX,100 + buttonmoveY)

		self.btn4 = QPushButton("Experiment", self)
		self.btn4.clicked.connect(self.Diversity_Plot_DF)
		self.btn4.clicked.connect(self.show_experiment)
		self.btn4.clicked.connect(self.Plot_Experiment)
		self.btn4.clicked.connect(self.hide_sample)
		self.btn4.clicked.connect(self.hide_replicate)
		self.btn4.clicked.connect(self.hide_diversity)
		# self.btn4.resize(buttonX,buttonY)
		# self.btn4.move(320 + buttonmoveX,100 + buttonmoveY)
		
		self.btn5 = QPushButton("Diversity", self)
		self.btn5.clicked.connect(self.Diversity_Plot_DF)
		self.btn5.clicked.connect(self.show_diversity)
		self.btn5.clicked.connect(self.Plot_Diversity)
		self.btn5.clicked.connect(self.hide_sample)
		self.btn5.clicked.connect(self.hide_replicate)
		self.btn5.clicked.connect(self.hide_experiment)



		# self.btn5.resize(buttonX,buttonY)
		# self.btn5.move(470 + buttonmoveX,100 + buttonmoveY)
		

		
		extractAction = QAction(QIcon('zoom.png'), 'Quit', self)
		extractAction.triggered.connect(self.close_application)
		


########## Grid Region 1 #######################

		layoutRegion1 = QGridLayout()
		layoutRegion1.addWidget(self.btn2,0,1)

		layoutRegion1.addWidget(self.btn3,0,2)
		layoutRegion1.addWidget(self.btn4,0,3)
		layoutRegion1.addWidget(self.btn5,0,4)




##################### Right hand side of interface (Region 2) #########################
		self.comboBox = QComboBox(self)
		self.comboAnalyte = QComboBox(self)
		self.comboBox.currentTextChanged.connect(self.reset_all)
		self.comboBox.currentTextChanged.connect(self.Sample_Plot_DF)

		self.comboBox.currentTextChanged.connect(self.Plot_ms1)


		n = 0
		while n < len(sample_names):
			self.comboBox.addItem(sample_names[n])
			
			n+=1	



		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()

		self.chooseSample = QLabel("Choose Sample", self)

		
		# self.chooseAnalyte = QLabel("Choose Analyte", self)

		self.chooseAnalyte = QComboBox(self)
		self.chooseAnalyte.currentTextChanged.connect(self.fill_combobox)
		self.chooseAnalyte.currentTextChanged.connect(self.Sample_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_Sample)
		self.chooseAnalyte.currentTextChanged.connect(self.Replicate_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_Replicate)
		self.chooseAnalyte.currentTextChanged.connect(self.ms1_Plot_DF)
		self.chooseAnalyte.currentTextChanged.connect(self.Plot_ms1)
		self.chooseAnalyte.addItem('Sample Analyte ID')
		self.chooseAnalyte.addItem('Replicate Analyte ID')
		self.chooseAnalyte.addItem('Experiment Analyte ID')
		self.showBlanks = QLabel("Show Blank Data", self)

		
		self.showNull = QLabel("Show Null Data", self)


		comboText=self.comboBox.currentText()
		analyteID=self.chooseAnalyte.currentText()


		total_df = import_dataframe(input_structure, input_type, comboText)

		# sorted_df = total_df['analyte_id'].replace('',np.nan, inplace=True)
		# sorted_df = total_df.dropna(subset=['analyte_id'],inplace=True)
		# print(sorted_df)


		replicate_id_df = total_df.dropna(subset=['replicate_analyte_id'])
		sorted_df = total_df.dropna(subset=['analyte_id'])
		sorted_df=total_df[total_df.analyte_id != None]
		sorted_df = sorted_df.sort_values(by=['analyte_id'])
		# total_df['analyte_id'].replace('',np.nan, inplace=True)
		mz_array=[]
		analyte_array=[]

		for analyte_id in sorted_df.analyte_id.unique():
			mz_array.append([sorted_df[sorted_df["analyte_id"] == analyte_id]["mz"].to_numpy(),
		                     sorted_df[sorted_df["analyte_id"] == analyte_id]["intensity"].to_numpy()])

			analyte_array.append(analyte_id)


		replicate_id_df = replicate_id_df.dropna(subset=['replicate_analyte_id'])
		replicate_id_df = replicate_id_df[replicate_id_df.replicate_analyte_id != None]
		replicate_id_df = replicate_id_df.sort_values(by=['replicate_analyte_id'])

		replicate_analyte_array=[]





		for analyte_id in total_df.analyte_id.unique():
			replicate_analyte_array.append(analyte_id)










		self.comboAnalyte.addItem('Show All')

		# if analyteID == 'Analyte ID':
		n = 0
		while n < len(analyte_array):

			self.comboAnalyte.addItem(str(analyte_array[n]))
			
			n+=1
		# else:
		# 	pass	
		# if analyteID == 'replicate_analyte_id':
		# 	n = 0
		# 	while n < len(replicate_analyte_array):

		# 		self.comboAnalyte.addItem(str(replicate_analyte_array[n]))
				
		# 		n+=1	
		# else:
		# 	pass
		

		self.massRangeLabel = QLabel("m/z Range", self)

  
		self.mMinBox = QLineEdit(self)
		self.mMinBox.setText('0')
		self.mMinBox.setValidator(QIntValidator())
		self.mMinBox.setAlignment(Qt.AlignRight)
		self.mMinBox.setValidator(QRegExpValidator())


		self.mMinBox.textChanged.connect(self.mMinLimit)



		self.mMaxBox = QLineEdit(self)
		self.mMaxBox.setText('1200')
		self.mMaxBox.setValidator(QIntValidator())
		self.mMaxBox.setValidator(QRegExpValidator())
		self.mMaxBox.setAlignment(Qt.AlignRight)
		self.mMaxBox.textChanged.connect(self.mMaxLimit)



		self.RTLabel = QLabel("rt Range", self)


		self.rtMinBox = QLineEdit(self)
		self.rtMinBox.setText('0')
		self.rtMinBox.setValidator(QIntValidator())
		self.rtMinBox.setValidator(QRegExpValidator())
		self.rtMinBox.setAlignment(Qt.AlignRight)
		self.rtMinBox.textChanged.connect(self.rtMinLimit)


		self.rtMaxBox = QLineEdit(self)
		self.rtMaxBox.setText('7')
		self.rtMaxBox.setValidator(QIntValidator())
		self.rtMaxBox.setValidator(QRegExpValidator())
		self.rtMaxBox.setAlignment(Qt.AlignRight)
		self.rtMaxBox.textChanged.connect(self.rtMaxLimit)


		# self.mtoLabel = QLabel("to",self)
		# self.mtoLabel.resize(50,50)
		# self.mtoLabel.move(2235,650 + move_x)

		# self.rttoLabel = QLabel("to",self)
		# self.rttoLabel.resize(50,50)
		# self.rttoLabel.move(2235,800 + move_x)

		self.massCheckbox = QCheckBox('Link m/z and rt',self)



		if self.massCheckbox.isChecked() == False:
			self.mMinBox.setText('0')
			self.mMaxBox.setText('1200')


		
		# self.rtCheckbox = QCheckBox('Link m/z',self)
		# self.rtCheckbox.resize(500,50)
		# self.rtCheckbox.move(2420,800 + move_x)



		self.resetButton = QPushButton("Reset All", self)
		self.resetButton.clicked.connect(self.reset_all)








		self.resetMassButton = QPushButton("Reset", self)
		self.resetMassButton.clicked.connect(self.reset_massRange)







		self.resetRTButton = QPushButton("Reset", self)
		self.resetRTButton.clicked.connect(self.reset_rtRange)




		# if Sample_State == False:
		self.massCheckbox.stateChanged.connect(self.Plot_Sample)
		self.massCheckbox.stateChanged.connect(self.Sample_Plot_DF)
		# self.rtCheckbox.stateChanged.connect(self.Plot_Sample)
		# self.rtCheckbox.stateChanged.connect(self.Sample_Plot_DF)
		self.resetButton.clicked.connect(self.Sample_Plot_DF)
		self.resetButton.clicked.connect(self.Plot_Sample)
		self.resetMassButton.clicked.connect(self.Sample_Plot_DF)
		self.resetMassButton.clicked.connect(self.Plot_Sample)
		self.resetRTButton.clicked.connect(self.Sample_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Sample)
		self.comboAnalyte.currentTextChanged.connect(self.Sample_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Sample)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_ms1)
		# else:
		# 	pass
		#if Sample_State == True and Replicate_State == False:
	# if Replicate_State == False:
		self.massCheckbox.stateChanged.connect(self.Plot_Replicate)
		self.massCheckbox.stateChanged.connect(self.Replicate_Plot_DF)
		# self.rtCheckbox.stateChanged.connect(self.Plot_Replicate)
		# self.rtCheckbox.stateChanged.connect(self.Replicate_Plot_DF)
		self.resetButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetButton.clicked.connect(self.Plot_Replicate)
		self.resetMassButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetMassButton.clicked.connect(self.Plot_Replicate)
		self.resetRTButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Replicate)
		self.resetRTButton.clicked.connect(self.Replicate_Plot_DF)
		self.resetRTButton.clicked.connect(self.Plot_Replicate)
		self.comboAnalyte.currentTextChanged.connect(self.Replicate_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Replicate)
		# self.comboAnalyte.currentTextChanged.connect(self.Experiment_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Diversity_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Experiment)
		self.comboAnalyte.currentTextChanged.connect(self.Diversity_Plot_DF)
		self.comboAnalyte.currentTextChanged.connect(self.Plot_Diversity)
		# else:
		# 	pass

################ Center Left of interface (Region 3) ###################
		


		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()


		
############## Null On Off ####################

		# creating a push button 
		self.Nullbutton = QPushButton("On", self) 

		# setting geometry of button 


		# setting checkable to true 
		self.Nullbutton.setCheckable(True) 

		# setting calling method by button 

		self.Nullbutton.clicked.connect(self.nullButton)


		# setting default color of button to light-grey 
		self.Nullbutton.setStyleSheet("background-color : lightgreen") 

############## Blank On Off ####################

		# creating a push button 
		self.Blankbutton = QPushButton("On", self) 

		# setting geometry of button 


		# setting checkable to true 
		self.Blankbutton.setCheckable(True) 

		# setting calling method by button 

		self.Blankbutton.clicked.connect(self.blankButton)


		# setting default color of button to light-grey 
		self.Blankbutton.setStyleSheet("background-color : lightgreen") 

		self.Blankbutton.clicked.connect(self.Plot_Sample)
		self.Nullbutton.clicked.connect(self.Plot_Sample)

		self.Nullbutton.clicked.connect(self.Plot_Replicate)
		self.Blankbutton.clicked.connect(self.Plot_Replicate)

		self.Blankbutton.clicked.connect(self.Plot_Experiment)

		self.Blankbutton.clicked.connect(self.Plot_Diversity)


########### Grid Region 2 ###########################
		self.Filters = QLabel("Filters")
		self.Filters.setFont(QFont('Times', 12))
		self.Line1 = QLabel("__________________________________________________")
		self.Line2 = QLabel("__________________________________________________")
		self.Line3 = QLabel("__________________________________________________")
		self.Line4 = QLabel("__________________________________________________")
		self.Blank1 = QLabel("   ")
		layoutRegion2 = QGridLayout()

		# layoutRegion2.setRowStretch(0, 1)
		layoutRegion2.addWidget(self.chooseSample,0,0)
		layoutRegion2.addWidget(self.comboBox,0,1)
		layoutRegion2.addWidget(self.chooseAnalyte,1,0)
		layoutRegion2.addWidget(self.comboAnalyte,1,1)

		layoutRegion4 = QGridLayout()
		layoutRegion4.addWidget(self.Blank1,0,0,1,3,QtCore.Qt.AlignCenter)
		layoutRegion4.addWidget(self.Filters,1,0,1,3)
		layoutRegion4.addWidget(self.Line1,2,0,1,3,QtCore.Qt.AlignTop)
		layoutRegion4.addWidget(self.showNull,3,0)
		layoutRegion4.addWidget(self.Nullbutton,3,1)
		layoutRegion4.addWidget(self.showBlanks,4,0)
		layoutRegion4.addWidget(self.Blankbutton,4,1)


		layoutRegion4.addWidget(self.Line2,5,0,1,3)
		layoutRegion4.addWidget(self.massRangeLabel,6,0)
		layoutRegion4.addWidget(self.resetMassButton,6,1)
		layoutRegion4.addWidget(self.mMinBox,7,0)
		layoutRegion4.addWidget(self.mMaxBox,7,1)


		layoutRegion4.addWidget(self.Line3,8,0,1,3)
		layoutRegion4.addWidget(self.RTLabel,9,0)
		layoutRegion4.addWidget(self.resetRTButton,9,1)
		layoutRegion4.addWidget(self.rtMinBox,10,0)
		layoutRegion4.addWidget(self.rtMaxBox,10,1)
		layoutRegion4.addWidget(self.massCheckbox,12,1)

		layoutRegion4.addWidget(self.Line4,11,0,1,3)
		layoutRegion4.addWidget(self.resetButton,12,0)







		comboText = self.comboBox.currentText()


################################################# Sample
		self.mz_plot = pg.PlotWidget(self)
		
		self.mz_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_plot.setTitle(comboText)
		self.mz_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})

		self.mz_plot.setStyleSheet("background-color: white;  border:1px solid black")

		self.rt_plot = pg.PlotWidget(self)
		self.rt_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		self.rt_plot.setTitle(comboText)
		self.ms1_plot = pg.PlotWidget(self)
		self.ms1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.ms1_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'})
		self.ms1_plot.setTitle("MS<sup>1</sup> Spectra")
		self.ms1_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
		self.ms2_plot = pg.PlotWidget(self)
		self.ms2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.ms2_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.ms2_plot.setTitle("MS<sup>2</sup> Spectra")
		self.ms2_plot.setLabel(axis='left',text='Intensity',**{ 'font-size': '10pt'})
################################################# Replicate

############# Replicate #############################


		self.mz_r1_plot = pg.PlotWidget(self)
		self.mz_r1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r1_plot.setTitle("Replicate 1 m/z vs intensity")
		self.mz_r1_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r1_plot.setLabel(axis='left',text='Intensity')
		
		self.mz_r2_plot = pg.PlotWidget(self)
		self.mz_r2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r2_plot.setTitle("Replicate 2 m/z vs intensity")
		self.mz_r2_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r2_plot.setLabel(axis='left',text='Intensity')
		
		self.mz_r3_plot = pg.PlotWidget(self)
		self.mz_r3_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.mz_r3_plot.setTitle("Replicate 3 m/z vs intensity")
		self.mz_r3_plot.setLabel(axis='bottom',text='m/z',**{ 'font-size': '10pt'}, **{'font': 'Italicized'})
		self.mz_r3_plot.setLabel(axis='left',text='Intensity')

		self.rt_r1_plot = pg.PlotWidget(self)
		self.rt_r1_plot.setTitle("Replicate 1 rt vs intesnity")
		self.rt_r1_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_r1_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r1_plot.setLabel(axis='left',text='Intensity')
		
		self.rt_r2_plot = pg.PlotWidget(self)
		self.rt_r2_plot.setTitle("Replicate 2 rt vs intesnity")
		self.rt_r2_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r2_plot.setLabel(axis='left',text='Intensity')
		self.rt_r2_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		self.rt_r3_plot = pg.PlotWidget(self)
		self.rt_r3_plot.setTitle("Replicate 3 rt vs intesnity")
		self.rt_r3_plot.setLabel(axis='bottom',text='rt',**{ 'font-size': '10pt'})
		self.rt_r3_plot.setLabel(axis='left',text='Intensity')
		self.rt_r3_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")



############# Experiment Blank On #############################

		self.ex_bOn_plot = pg.PlotWidget(self)
		self.ex_bOn_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")





		self.div_plot = pg.PlotWidget(self)
		self.div_plot.setStyleSheet("background-color: white; margin:0.001px; border:1px solid black; ")
		# self.ex_bOn_plot.move(20,250 + move_y)

		
		if len(sample_names) < 16:
			self.ex_bOn_plot.resize(2000,100*len(sample_names))
		else:
			self.ex_bOn_plot.resize(2000,1500)
			self.ex_bOn_plot.setRange(yRange=(0,15), padding = 0)
		self.ex_bOn_plot.setLabel(axis='bottom',text='Max Mass')



		if len(sample_names) < 16:
			self.div_plot.resize(2000,100*len(sample_names))
		else:
			self.div_plot.resize(2000,1500)
			self.div_plot.setRange(yRange=(0,15), padding = 0)
		self.div_plot.setLabel(axis='bottom',text='Analyte ID')
##############################################################

		self.ms2analytelabel = QLabel(self)
		pixmap = QPixmap('MS2Analyte_logo.png')

		self.ms2analytelabel.setPixmap(pixmap)
		layoutRegion7 = QGridLayout()

		layoutRegion7.addWidget(self.ms2analytelabel)
##############################################################

############ Grid Region 3 ###################################

		layoutRegion3 = QGridLayout()
		
		layoutRegion3.addWidget(self.mz_plot,0,1)
		layoutRegion3.addWidget(self.rt_plot,1,1)

		layoutRegion5 = QGridLayout()
		layoutRegion5.addWidget(self.ms1_plot,0,1)
		layoutRegion5.addWidget(self.ms2_plot,1,1)

		

		layoutRegion3.addWidget(self.mz_r1_plot,1,0)
		layoutRegion3.addWidget(self.mz_r2_plot,2,0)
		layoutRegion3.addWidget(self.mz_r3_plot,3,0)
		layoutRegion3.addWidget(self.rt_r1_plot,4,0)
		layoutRegion3.addWidget(self.rt_r2_plot,5,0)
		layoutRegion3.addWidget(self.rt_r3_plot,6,0)

		layoutRegion6 = QGridLayout()

		layoutRegion6.addWidget(self.ex_bOn_plot,1,0)

		layoutRegion8 = QGridLayout()

		layoutRegion8.addWidget(self.div_plot,1,0)
############# Create Main Grid ###############################

		# g.addWidget(button, N, 0, 1, N, QtCore.Qt.AlignCenter)
		layoutMain.setColumnStretch(0, 4)
		layoutMain.setColumnStretch(1, 2)
		# layoutMain.setRowStretch(1, 0.5)
		# layoutRegion2.setColumnStretch(0,4)
		layoutMain.addLayout(layoutRegion3,1,0,QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion5,1,1,QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion1,0,0,QtCore.Qt.AlignLeft)

		# layoutMain.setSpacing(0)
		layoutMain.addLayout(layoutRegion2,0,1, 1, 1, QtCore.Qt.AlignRight)
		layoutMain.addLayout(layoutRegion7,0,2, QtCore.Qt.AlignCenter)
		layoutMain.addLayout(layoutRegion4,1,2, 1, 1,QtCore.Qt.AlignTop)
		layoutMain.addLayout(layoutRegion6,1,0, 1,1)
		layoutMain.addLayout(layoutRegion8,1,0, 1,1)
############################################################################



		# layoutMain.addLayout(layoutRegion3_R,1,0)
		# layoutMain.addLayout(layoutRegion3_E,1,0)
		self.horizontalGroupBox.setLayout(layoutMain)







################### Plots #################################


	def fill_combobox(self):

		comboText = self.comboBox.currentText()
		chooseAnalyte = self.chooseAnalyte.currentText()
		total_df = import_dataframe(input_structure, input_type, comboText)

		# sorted_df = total_df['analyte_id'].replace('',np.nan, inplace=True)
		# sorted_df = total_df.dropna(subset=['analyte_id'],inplace=True)
		# print(sorted_df)


		replicate_id_df = total_df.dropna(subset=['replicate_analyte_id'])
		experiment_id_df = total_df.dropna(subset=['experiment_analyte_id'])
		sorted_df = total_df.dropna(subset=['analyte_id'])
		sorted_df=total_df[total_df.analyte_id != None]
		sorted_df = sorted_df.sort_values(by=['analyte_id'])
		# total_df['analyte_id'].replace('',np.nan, inplace=True)
		mz_array=[]
		analyte_array=[]

		for analyte_id in sorted_df.analyte_id.unique():
			mz_array.append([sorted_df[sorted_df["analyte_id"] == analyte_id]["mz"].to_numpy(),
		                     sorted_df[sorted_df["analyte_id"] == analyte_id]["intensity"].to_numpy()])

			analyte_array.append(analyte_id)


		replicate_id_df = replicate_id_df.dropna(subset=['replicate_analyte_id'])
		replicate_id_df = replicate_id_df[replicate_id_df.replicate_analyte_id != None]
		replicate_id_df = replicate_id_df.sort_values(by=['replicate_analyte_id'])

		replicate_analyte_array=[]



		for analyte_id in replicate_id_df.replicate_analyte_id.unique():
			replicate_analyte_array.append(analyte_id)


		experiment_id_df = experiment_id_df.dropna(subset=['experiment_analyte_id'])
		experiment_id_df = experiment_id_df[experiment_id_df.experiment_analyte_id != None]
		experiment_id_df = experiment_id_df.sort_values(by=['experiment_analyte_id'])

		experiment_analyte_array=[]



		for analyte_id in experiment_id_df.experiment_analyte_id.unique():
			experiment_analyte_array.append(analyte_id)
		# self.comboAnalyte.clear()



		if chooseAnalyte == "Sample Analyte ID":


			self.comboAnalyte.blockSignals(True)
			self.comboAnalyte.clear()
			self.comboAnalyte.addItem("Show All")
			n = 0
			while n < len(analyte_array):

				self.comboAnalyte.addItem(str(analyte_array[n]))
				
				n+=1
				self.comboAnalyte.blockSignals(False)
		else:
			pass



		if chooseAnalyte == "Replicate Analyte ID":


			self.comboAnalyte.blockSignals(True)
			self.comboAnalyte.clear()
			self.comboAnalyte.addItem("Show All")

			n = 0
			while n < len(replicate_analyte_array):

				self.comboAnalyte.addItem(str(replicate_analyte_array[n]))
				
				n+=1
				self.comboAnalyte.blockSignals(False)
		
		else:
			pass
		if chooseAnalyte == "Experiment Analyte ID":


			self.comboAnalyte.blockSignals(True)
			self.comboAnalyte.clear()
			self.comboAnalyte.addItem("Show All")

			n = 0
			while n < len(experiment_analyte_array):

				self.comboAnalyte.addItem(str(experiment_analyte_array[n]))
				
				n+=1
				self.comboAnalyte.blockSignals(False)
		
		else:
			pass
	def Plot_Sample(self):
		
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		print('Loading Sample Tab')
		


		if Sample_State == False:
			chooseAnalyte = self.chooseAnalyte.currentText()
			self.mz_plot.disableAutoRange(axis=None)
			self.rt_plot.disableAutoRange(axis=None)
			comboAnalyte=self.comboAnalyte.currentText()
		
			print('Loading Sample Tab')
			
			############### Fetch data for mz_plot and plot sample from comboBox ###############

			sorted_df_mz, sorted_df_rt= self.Sample_Plot_DF()


			mz_array = []
			analyte_mz_array = []
			rt_array = []
			analyte_rt_array = []


			if chooseAnalyte == "Sample Analyte ID":
				for analyte_id in sorted_df_mz.analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz[sorted_df_mz["analyte_id"] == analyte_id]["intensity"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				for analyte_id in sorted_df_rt.analyte_id.unique():
					rt_array.append([sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt[sorted_df_rt["analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array.append(analyte_id)


				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

				analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]



				analyte_colour_mz_array = [0]*len(analyte_mz_array)
				analyte_colour_rt_array = [0]*len(analyte_rt_array)

				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)
					i+=1


				i = 0
				while i < len(analyte_rt_array):

					analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
					i+=1
				print('Loading Sample Tab')
			if chooseAnalyte == "Replicate Analyte ID":

				for analyte_id in sorted_df_mz.replicate_analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz[sorted_df_mz["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				for analyte_id in sorted_df_rt.replicate_analyte_id.unique():
					rt_array.append([sorted_df_rt[sorted_df_rt["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt[sorted_df_rt["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array.append(analyte_id)


				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

				analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]

				analyte_mz_array = [0 if str(i) == 'None' else i for i in analyte_mz_array]

				analyte_rt_array = [0 if str(i) == 'None' else i for i in analyte_rt_array]

				analyte_colour_mz_array = [0]*len(analyte_mz_array)
				analyte_colour_rt_array = [0]*len(analyte_rt_array)


				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)

					i+=1


				i = 0
				while i < len(analyte_rt_array):

					analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
					i+=1

				print('Loading Sample Tab')

			if chooseAnalyte == "Experiment Analyte ID":

				for analyte_id in sorted_df_mz.experiment_analyte_id.unique():

					mz_array.append([sorted_df_mz[sorted_df_mz["experiment_analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz[sorted_df_mz["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])

					analyte_mz_array.append(analyte_id)



				for analyte_id in sorted_df_rt.experiment_analyte_id.unique():
					rt_array.append([sorted_df_rt[sorted_df_rt["experiment_analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt[sorted_df_rt["experiment_analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array.append(analyte_id)


				analyte_mz_array = [0 if str(i) == 'nan' else i for i in analyte_mz_array]

				analyte_rt_array = [0 if str(i) == 'nan' else i for i in analyte_rt_array]

				analyte_mz_array = [0 if str(i) == 'None' else i for i in analyte_mz_array]

				analyte_rt_array = [0 if str(i) == 'None' else i for i in analyte_rt_array]

				analyte_colour_mz_array = [0]*len(analyte_mz_array)
				analyte_colour_rt_array = [0]*len(analyte_rt_array)


				i = 0
				while i < len(analyte_mz_array):
					analyte_colour_mz_array[i] = int((analyte_mz_array[i] + 1) % 299)

					i+=1


				i = 0
				while i < len(analyte_rt_array):

					analyte_colour_rt_array[i] = int((analyte_rt_array[i] + 1) % 299)
					i+=1

				print('Loading Sample Tab')




			self.rt_plot.clear()
			self.mz_plot.clear()

			if comboAnalyte != "Show All":
				alpha=1
				print('Loading Sample Tab..')

				colour=0
				for analyte in mz_array:


					#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					#colour+=1

					if analyte_mz_array[colour] == float(comboAnalyte):


						self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],hues=20,values=5,alpha=255),symbol='o', symbolSize=4)
						
					else:
						pass


					self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					
					colour+=1

				colour=0
				for analyte in rt_array:


					if analyte_rt_array[colour] == float(comboAnalyte):


						self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], hues=20,values=5, alpha=255),symbol='o', symbolSize=4)
						

					else:
						pass

					self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], hues=20,values=5, alpha=alpha),symbol='o', symbolSize=3)
					colour+=1


				print('Loading Sample Tab...')
			else:
				alpha=255
				colour=0
				for analyte in mz_array:


					self.test = self.mz_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_mz_array[colour],hues=20,values=5),symbol='o', symbolSize=3)
					#test.setAlpha(0, False)
					colour+=1
				colour=0
				for analyte in rt_array:





					self.rt_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_rt_array[colour], hues=20,values=5),symbol='o', symbolSize=3)
					colour+=1
				print('Loading Sample Tab...')

			############### Fetch data for rt plot and plot sample from comboBox ###############




			



			#analyte_array[0] = 0
			#print(analyte_array)






			self.mz_plot.enableAutoRange(axis=None)
			self.rt_plot.enableAutoRange(axis=None)
			state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())

			self.mz_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
			self.rt_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
		else:
			pass

		print('Sample Tab Loaded Succesfully')


##################################################################################################################	

		return 
		

	def Plot_Replicate(self):
		print('Loading Replicate Tab')
		chooseAnalyte = self.chooseAnalyte.currentText()
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		print(Replicate_State)
		print(Sample_State)
		print(Experiment_State)


		print('Loading Replicate Tab')
		if Replicate_State == False:


			self.mz_r1_plot.disableAutoRange(axis=None)
			self.mz_r2_plot.disableAutoRange(axis=None)
			self.mz_r3_plot.disableAutoRange(axis=None)
			self.rt_r1_plot.disableAutoRange(axis=None)
			self.rt_r2_plot.disableAutoRange(axis=None)
			self.rt_r3_plot.disableAutoRange(axis=None)
			comboAnalyte=self.comboAnalyte.currentText()
		

			
			############### Fetch data for mz_plot and plot sample from comboBox ###############

			sorted_df_mz_r1, sorted_df_rt_r1,sorted_df_mz_r2, sorted_df_rt_r2,sorted_df_mz_r3, sorted_df_rt_r3 = self.Replicate_Plot_DF()


			mz_array_r1 = []
			analyte_mz_array_r1 = []
			rt_array_r1 = []
			analyte_rt_array_r1 = []


			if chooseAnalyte == "Sample Analyte ID":
				print('Loading Replicate Tab')
				for analyte_id in sorted_df_mz_r1.analyte_id.unique():

					mz_array_r1.append([sorted_df_mz_r1[sorted_df_mz_r1["analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r1[sorted_df_mz_r1["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r1.append(analyte_id)






				for analyte_id in sorted_df_rt_r1.analyte_id.unique():
					rt_array_r1.append([sorted_df_rt_r1[sorted_df_rt_r1["analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r1[sorted_df_rt_r1["analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r1.append(analyte_id)


				analyte_mz_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r1]

				analyte_rt_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r1]



				analyte_colour_mz_array_r1 = [0]*len(analyte_mz_array_r1)
				analyte_colour_rt_array_r1 = [0]*len(analyte_rt_array_r1)

				
				for i in range(0, len(analyte_mz_array_r1)):
					analyte_colour_mz_array_r1[i] = int((analyte_mz_array_r1[i] + 1) % 299)
					


				
				for i in range(0, len(analyte_rt_array_r1)):

					analyte_colour_rt_array_r1[i] = int((analyte_rt_array_r1[i] + 1) % 299)
					
				print('Loading Replicate Tab')


			if chooseAnalyte == "Replicate Analyte ID":
				print('Loading Replicate Tab')
				for analyte_id in sorted_df_mz_r1.replicate_analyte_id.unique():

					mz_array_r1.append([sorted_df_mz_r1[sorted_df_mz_r1["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r1[sorted_df_mz_r1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r1.append(analyte_id)




				print('Loading Replicate Tab')

				for analyte_id in sorted_df_rt_r1.replicate_analyte_id.unique():
					rt_array_r1.append([sorted_df_rt_r1[sorted_df_rt_r1["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r1[sorted_df_rt_r1["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r1.append(analyte_id)


				analyte_mz_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r1]

				analyte_rt_array_r1 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r1]

				analyte_mz_array_r1 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r1]

				analyte_rt_array_r1 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r1]

				analyte_colour_mz_array_r1 = [0]*len(analyte_mz_array_r1)
				analyte_colour_rt_array_r1 = [0]*len(analyte_rt_array_r1)

				i = 0
				while i < len(analyte_mz_array_r1):
					analyte_colour_mz_array_r1[i] = int((analyte_mz_array_r1[i] + 1) % 299)

					i+=1


				i = 0
				while i < len(analyte_rt_array_r1):

					analyte_colour_rt_array_r1[i] = int((analyte_rt_array_r1[i] + 1) % 299)
					i+=1




				print('Loading Replicate Tab')




			self.rt_r1_plot.clear()
			self.mz_r1_plot.clear()

			if comboAnalyte != "Show All":
				print('Loading Replicate Tab')
				alpha=1


				colour=0
				for analyte in mz_array_r1:



					if analyte_mz_array_r1[colour] == float(comboAnalyte):


						self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],hues=20,values=5,alpha=255),symbol='o', symbolSize=4)
						
					else:
						pass


					self.mz_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					
					colour+=1

				colour=0
				for analyte in rt_array_r1:


					if analyte_rt_array_r1[colour] == float(comboAnalyte):


						self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], hues=20,values=5, alpha=255),symbol='o', symbolSize=4)
						

					else:
						pass

					self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], hues=20,values=5, alpha=alpha),symbol='o', symbolSize=3)
					colour+=1



			else:
				alpha=255
				colour=0
				for analyte in mz_array_r1:


					self.test = self.mz_r1_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r1[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r1[colour],hues=20,values=5),symbol='o', symbolSize=3)
					#test.setAlpha(0, False)
					colour+=1
				colour=0
				for analyte in rt_array_r1:





					self.rt_r1_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r1[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r1[colour], hues=20,values=5),symbol='o', symbolSize=3)
					colour+=1
				print('Loading Replicate Tab')

			############### Fetch data for rt plot and plot sample from comboBox ###############




			mz_array_r2 = []
			analyte_mz_array_r2 = []
			rt_array_r2 = []
			analyte_rt_array_r2 = []

			if chooseAnalyte == "Sample Analyte ID":
				print('Loading Replicate Tab')
				for analyte_id in sorted_df_mz_r2.analyte_id.unique():

					mz_array_r2.append([sorted_df_mz_r2[sorted_df_mz_r2["analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r2[sorted_df_mz_r2["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r2.append(analyte_id)


					



				for analyte_id in sorted_df_rt_r2.analyte_id.unique():
					rt_array_r2.append([sorted_df_rt_r2[sorted_df_rt_r2["analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r2[sorted_df_rt_r2["analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r2.append(analyte_id)


				analyte_mz_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r2]

				analyte_rt_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r2]



				analyte_colour_mz_array_r2 = [0]*len(analyte_mz_array_r2)
				analyte_colour_rt_array_r2 = [0]*len(analyte_rt_array_r2)
				print('Loading Replicate Tab')
				i = 0
				while i < len(analyte_mz_array_r2):
					analyte_colour_mz_array_r2[i] = int((analyte_mz_array_r2[i] + 1) % 299)
					i+=1


				i = 0
				while i < len(analyte_rt_array_r2):

					analyte_colour_rt_array_r2[i] = int((analyte_rt_array_r2[i] + 1) % 299)
					i+=1


			if chooseAnalyte == "Replicate Analyte ID":
				print('Loading Replicate Tab')
				for analyte_id in sorted_df_mz_r2.replicate_analyte_id.unique():

					mz_array_r2.append([sorted_df_mz_r2[sorted_df_mz_r2["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r2[sorted_df_mz_r2["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r2.append(analyte_id)


					



				for analyte_id in sorted_df_rt_r2.replicate_analyte_id.unique():
					rt_array_r2.append([sorted_df_rt_r2[sorted_df_rt_r2["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r2[sorted_df_rt_r2["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r2.append(analyte_id)


				analyte_mz_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r2]

				analyte_rt_array_r2 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r2]

				analyte_mz_array_r2 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r2]

				analyte_rt_array_r2 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r2]

				analyte_colour_mz_array_r2 = [0]*len(analyte_mz_array_r2)
				analyte_colour_rt_array_r2 = [0]*len(analyte_rt_array_r2)

				i = 0
				while i < len(analyte_mz_array_r2):
					analyte_colour_mz_array_r2[i] = int((analyte_mz_array_r2[i] + 1) % 299)
					i+=1


				i = 0
				while i < len(analyte_rt_array_r2):

					analyte_colour_rt_array_r2[i] = int((analyte_rt_array_r2[i] + 1) % 299)
					i+=1






			self.rt_r2_plot.clear()
			self.mz_r2_plot.clear()
			print('Loading Replicate Tab')
			if comboAnalyte != "Show All":
				print('Loading Replicate Tab')
				alpha=1


				colour=0
				for analyte in mz_array_r2:


					#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					#colour+=1

					if analyte_mz_array_r2[colour] == float(comboAnalyte):


						self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],hues=20,values=5,alpha=255),symbol='o', symbolSize=4)
						
					else:
						pass


					self.mz_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					
					colour+=1

				colour=0
				for analyte in rt_array_r2:


					if analyte_rt_array_r2[colour] == float(comboAnalyte):


						self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], hues=20,values=5, alpha=255),symbol='o', symbolSize=4)
						

					else:
						pass

					self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], hues=20,values=5, alpha=alpha),symbol='o', symbolSize=3)
					colour+=1


				print('Loading Replicate Tab')
			else:
				alpha=255
				colour=0
				for analyte in mz_array_r2:


					self.test = self.mz_r2_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r2[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r2[colour],hues=20,values=5),symbol='o', symbolSize=3)
					#test.setAlpha(0, False)
					colour+=1
				colour=0
				for analyte in rt_array_r2:





					self.rt_r2_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r2[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r2[colour], hues=20,values=5),symbol='o', symbolSize=3)
					colour+=1



			mz_array_r3 = []
			analyte_mz_array_r3 = []
			rt_array_r3 = []
			analyte_rt_array_r3 = []
			print('Loading Replicate Tab')
			if chooseAnalyte == "Sample Analyte ID":
				print('Loading Replicate Tab')
				for analyte_id in sorted_df_mz_r3.analyte_id.unique():

					mz_array_r3.append([sorted_df_mz_r3[sorted_df_mz_r3["analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r3[sorted_df_mz_r3["analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r3.append(analyte_id)


					



				for analyte_id in sorted_df_rt_r3.analyte_id.unique():
					rt_array_r3.append([sorted_df_rt_r3[sorted_df_rt_r3["analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r3[sorted_df_rt_r3["analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r3.append(analyte_id)


				analyte_mz_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r3]

				analyte_rt_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r3]



				analyte_colour_mz_array_r3 = [0]*len(analyte_mz_array_r3)
				analyte_colour_rt_array_r3 = [0]*len(analyte_rt_array_r3)

				i = 0
				while i < len(analyte_mz_array_r3):
					analyte_colour_mz_array_r3[i] = int((analyte_mz_array_r3[i] + 1) % 299)
					i+=1


				i = 0
				while i < len(analyte_rt_array_r3):

					analyte_colour_rt_array_r3[i] = int((analyte_rt_array_r3[i] + 1) % 299)
					i+=1


			if chooseAnalyte == "Replicate Analyte ID":
				for analyte_id in sorted_df_mz_r3.replicate_analyte_id.unique():

					mz_array_r3.append([sorted_df_mz_r3[sorted_df_mz_r3["replicate_analyte_id"] == analyte_id]["mz"].to_numpy(),
				                     sorted_df_mz_r3[sorted_df_mz_r3["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy(),])

					analyte_mz_array_r3.append(analyte_id)


					

				print('Loading Replicate Tab')

				for analyte_id in sorted_df_rt_r3.replicate_analyte_id.unique():
					rt_array_r3.append([sorted_df_rt_r3[sorted_df_rt_r3["replicate_analyte_id"] == analyte_id]["rt"].to_numpy(),
				                     sorted_df_rt_r3[sorted_df_rt_r3["replicate_analyte_id"] == analyte_id]["intensity"].to_numpy()])
					analyte_rt_array_r3.append(analyte_id)


				analyte_mz_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_mz_array_r3]

				analyte_rt_array_r3 = [0 if str(i) == 'nan' else i for i in analyte_rt_array_r3]

				analyte_mz_array_r3 = [0 if str(i) == 'None' else i for i in analyte_mz_array_r3]

				analyte_rt_array_r3 = [0 if str(i) == 'None' else i for i in analyte_rt_array_r3]

				analyte_colour_mz_array_r3 = [0]*len(analyte_mz_array_r3)
				analyte_colour_rt_array_r3 = [0]*len(analyte_rt_array_r3)

				i = 0
				while i < len(analyte_mz_array_r3):
					analyte_colour_mz_array_r3[i] = int((analyte_mz_array_r3[i] + 1) % 299)
					i+=1


				i = 0
				while i < len(analyte_rt_array_r3):

					analyte_colour_rt_array_r3[i] = int((analyte_rt_array_r3[i] + 1) % 299)
					i+=1





			self.rt_r3_plot.clear()
			self.mz_r3_plot.clear()
			print('Loading Replicate Tab')
			if comboAnalyte != "Show All":
				print('Loading Replicate Tab')
				alpha=1


				colour=0
				for analyte in mz_array_r3:


					#self.mz_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_array[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_array[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					#colour+=1

					if analyte_mz_array_r3[colour] == float(comboAnalyte):


						self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],hues=20,values=5,alpha=255),symbol='o', symbolSize=4)
						
					else:
						pass


					self.mz_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],hues=20,values=5,alpha=alpha),symbol='o', symbolSize=3)
					
					colour+=1

				colour=0
				for analyte in rt_array_r3:


					if analyte_rt_array_r3[colour] == float(comboAnalyte):


						self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour], hues=20, values=5, alpha=255), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], hues=20,values=5, alpha=255),symbol='o', symbolSize=4)
						

					else:
						pass

					self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour], hues=20, values=5, alpha=alpha), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], hues=20,values=5, alpha=alpha),symbol='o', symbolSize=3)
					colour+=1



			else:
				alpha=255
				colour=0
				for analyte in mz_array_r3:


					self.test = self.mz_r3_plot.plot( title="test", x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor( analyte_colour_mz_array_r3[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_mz_array_r3[colour],hues=20,values=5),symbol='o', symbolSize=3)
					#test.setAlpha(0, False)
					colour+=1
				colour=0
				for analyte in rt_array_r3:





					self.rt_r3_plot.plot( x=analyte[0], y=analyte[1], pen=None, symbolPen=pg.intColor(analyte_colour_rt_array_r3[colour], hues=20, values=5), symbolBrush=pg.intColor( analyte_colour_rt_array_r3[colour], hues=20,values=5),symbol='o', symbolSize=3)
					colour+=1
			#analyte_array[0] = 0
			#print(analyte_array)

			self.mz_r1_plot.getAxis('left').setStyle(showValues=False)
			self.mz_r2_plot.getAxis('left').setStyle(showValues=False)
			self.mz_r3_plot.getAxis('left').setStyle(showValues=False)
			self.rt_r1_plot.getAxis('left').setStyle(showValues=False)
			self.rt_r2_plot.getAxis('left').setStyle(showValues=False)
			self.rt_r3_plot.getAxis('left').setStyle(showValues=False)




			self.mz_r1_plot.enableAutoRange(axis=None)
			self.mz_r2_plot.enableAutoRange(axis=None)
			self.mz_r3_plot.enableAutoRange(axis=None)
			self.rt_r1_plot.enableAutoRange(axis=None)
			self.rt_r2_plot.enableAutoRange(axis=None)
			self.rt_r3_plot.enableAutoRange(axis=None)
			state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())

			self.mz_r1_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
			self.mz_r2_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
			self.mz_r3_plot.setXRange(float(state.return_massMin_state()),float(state.return_massMax_state()))
			self.rt_r1_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
			self.rt_r2_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
			self.rt_r3_plot.setXRange(float(state.return_rtMin_state()),float(state.return_rtMax_state()))
		else:
			pass
		print('Replicate Tab Loaded Succesfully')
		return

	def Plot_Experiment(self):

		# self.ex_bOn_plot.setYRange(min=-20,max=1)
		# self.ex_bOn_plot.setXRange(min=100,max=1200)
		self.ex_bOn_plot.disableAutoRange(axis=None)
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()



		if Experiment_State == False:
			comboAnalyte=self.comboAnalyte.currentText()

			

			samples, experiment_df = self.Experiment_Plot_DF()

			self.ex_bOn_plot.clear()
			
			for i in range(0, len(samples)):
				analyte_array = []
				experiment_sample_df = experiment_df[experiment_df.sample_name == samples[i]]


				experiment_array_BOn = []


				for experiment_analyte_id in experiment_sample_df.experiment_analyte_id.unique():
					experiment_array_BOn.append([experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["experiment_analyte_max_mass"].to_numpy(), 
									experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["sample_analyte_max_intensity"].to_numpy()])
					analyte_array.append(experiment_analyte_id)
					print(experiment_analyte_id)
				
				analyte_array = [0 if str(z) == 'nan' else z for z in analyte_array]



				analyte_colour_array = [0]*len(analyte_array)


				
				for a in range(0, len(analyte_array)):

					analyte_colour_array[a] = int((analyte_array[a] + 1) % 299)
					

				if comboAnalyte == "Show All":
					
					colour = 0
					for analyte in experiment_array_BOn:





						y=[i]*len(analyte[0])

						x=analyte[0]

						max_int = max(analyte[1])


						test = self.ex_bOn_plot.plot(title= "Test", x=x, y=y, pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5),symbol='o', symbolSize= 10)
						min_int = 50000
						int_x=40
						if max_int < min_int:
							test.setSymbolSize((max_int/230000)*int_x)
						else:
							test.setSymbolSize((max_int/230000)*4)

						colour+=1
				if comboAnalyte != "Show All":
					alpha=30
					colour = 0

					for analyte in experiment_array_BOn:




						max_int = max(analyte[1])
						if analyte_array[colour] == float(comboAnalyte):


							test = self.ex_bOn_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5),symbol='o', symbolSize= 10)
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)
						else:
							pass
						test = self.ex_bOn_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5, alpha = alpha), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5, alpha = alpha),symbol='o', symbolSize= 10)
						min_int = 50000
						int_x=40
						if max_int < min_int:
							test.setSymbolSize((max_int/230000)*int_x)
						else:
							test.setSymbolSize((max_int/230000)*4)

						colour+=1

				


			



			ticks = samples
			ticksdict = dict(enumerate(ticks))
			ax = self.ex_bOn_plot.getAxis('left')
			ax.setTickSpacing(1,1)
			ax.setTicks([ticksdict.items()])
		else:
			pass

		self.ex_bOn_plot.enableAutoRange(axis=None)
	def Plot_Diversity(self):
		print('Loading Diversity Tab')
		# self.div_plot.setYRange(min=-20,max=1)
		# self.div_plot.setXRange(min=0,max=140)
		self.div_plot.disableAutoRange(axis=None)
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		Diversity_State = tab_state.return_diversity_state()


		if Diversity_State == False:
			comboAnalyte=self.comboAnalyte.currentText()

			print('Loading Diversity Tab')

			samples, experiment_df = self.Diversity_Plot_DF()

			self.div_plot.clear()
			
			for i in range(0, len(samples)):
				print('Loading Diversity Tab')
				analyte_array = []
				experiment_sample_df = experiment_df[experiment_df.sample_name == samples[i]]


				experiment_array_BOn = []

				print('Loading Diversity Tab')
				for experiment_analyte_id in experiment_sample_df.experiment_analyte_id.unique():
					experiment_array_BOn.append([experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["experiment_analyte_id"].to_numpy(), 
									experiment_sample_df[experiment_sample_df["experiment_analyte_id"] == experiment_analyte_id]["sample_analyte_max_intensity"].to_numpy()])
					analyte_array.append(experiment_analyte_id)
				
				
				analyte_array = [0 if str(z) == 'nan' else z for z in analyte_array]



				analyte_colour_array = [0]*len(analyte_array)

				print('Loading Diversity Tab')
				
				for a in range(0, len(analyte_array)):

					analyte_colour_array[a] = int((analyte_array[a] + 1) % 299)
					

				if comboAnalyte == "Show All":
					print('Loading Diversity Tab')
					colour = 0
					for analyte in experiment_array_BOn:





						y=[i]*len(analyte[0])

						x=analyte[0]

						max_int = max(analyte[1])


						test = self.div_plot.plot(title= "Test", x=x, y=y, pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5),symbol='o', symbolSize= 10)
						min_int = 50000
						int_x=40
						if max_int < min_int:
							test.setSymbolSize((max_int/230000)*int_x)
						else:
							test.setSymbolSize((max_int/230000)*4)

						colour+=1
				if comboAnalyte != "Show All":
					print('Loading Diversity Tab')
					alpha=30
					colour = 0

					for analyte in experiment_array_BOn:




						max_int = max(analyte[1])
						if analyte_array[colour] == float(comboAnalyte):


							test = self.div_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5),symbol='o', symbolSize= 10)
							min_int = 50000
							int_x=40
							if max_int < min_int:
								test.setSymbolSize((max_int/230000)*int_x)
							else:
								test.setSymbolSize((max_int/230000)*4)
						else:
							pass
						test = self.div_plot.plot(title= "Test", x=analyte[0], y=[i]*len(analyte[0]), pen=None, symbolPen=pg.intColor(analyte_colour_array[colour], hues=20, values=5, alpha = alpha), symbolBrush=pg.intColor(analyte_colour_array[colour],hues=20,values=5, alpha = alpha),symbol='o', symbolSize= 10)
						min_int = 50000
						int_x=40
						if max_int < min_int:
							test.setSymbolSize((max_int/230000)*int_x)
						else:
							test.setSymbolSize((max_int/230000)*4)

						colour+=1

				


			


			print('PLOT SAMPLES', samples)
			
			ticks = samples
			ticksdict = dict(enumerate(ticks))
			ax = self.div_plot.getAxis('left')
			ax.setTickSpacing(1,1)
			ax.setTicks([ticksdict.items()])
		else:
			pass
		self.div_plot.enableAutoRange(axis=None)
		print('Diversity Tab Loaded Succesfully')
##################################################################################################################	

		return 
	def Plot_ms1(self):

		chooseAnalyte = self.chooseAnalyte.currentText()
		comboAnalyte=self.comboAnalyte.currentText()
		comboText=self.comboBox.currentText()
		self.ms1_plot.clear()
		df, df2 , df3= self.ms1_Plot_DF()

		# self.ms1_plot.disableAutoRange(axis=None)
		
		if chooseAnalyte == "Replicate Analyte ID":
			self.ms1_plot.setXRange(min=0,max=1200)
			self.ms1_plot.setYRange(min=5, max=100)
			ms1_data = []
			analyte_id_array = []

			if comboAnalyte != "Show All":
				print("Plot")
				for analyte_id in df.replicate_analyte_id.unique():
					ms1_data.append([df[df["replicate_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
				                     df[df["replicate_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
					analyte_id_array.append(analyte_id)
				analyte_colour_array = [0]*len(analyte_id_array)
				print('analyte_id_array',analyte_id_array)

				i = 0
				while i < len(analyte_id_array):
					analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
					i+=1



				colour=0
				self.ms1_plot.clear()
				for analyte in ms1_data:

					

					if float(comboAnalyte) == analyte_id_array[colour]:
						print("comboAnalyte = colour")
						pen = pg.mkPen(color=(0, 0, 0))
						self.ms1_plot.plot( x=analyte[1], y=analyte[0], pen=pen )


						i = 0
						while i < len(analyte[1]):
							np.round(analyte[1],2)
							np.round(analyte[0],2)
							if analyte[0][i] > 0:
								text = pg.TextItem(str(np.round(analyte[1][i],2)),color=(0,0,0))
								self.ms1_plot.addItem(text)


								text.setPos(analyte[1][i], analyte[0][i])
								i+=1
							else:
								i+=1
					else:

						pass
					colour +=1

			# self.ms1_plot.enableAutoRange(axis=None)




		if chooseAnalyte == "Sample Analyte ID":
			self.ms1_plot.setXRange(min=0,max=1200)
			self.ms1_plot.setYRange(min=50000, max=1000000)
			print('analyte_id')




			








			ms1_data = []
			analyte_id_array = []
			print(df2)
			if comboAnalyte != "Show All":
				print('COMBO ANALYTE',comboAnalyte)
				df=df2[df2.analyte_id == float(comboAnalyte)]
				max_scan = df[df.intensity == df.intensity.max()].scan
				print('MAX_SCAN',max_scan)
				print('COMBO ANALYTE',comboAnalyte)

				max_scan_value = max_scan.to_numpy()

				df=df[df.scan == max_scan_value[0]]
				df.sort_values('scan', inplace=True)





				mz_array = df['mz'].to_numpy()
				intensity_array = df['intensity'].to_numpy()

				print(mz_array)

				intensity_zeros = [0]*len(intensity_array)
				mass_lower = [0]*len(mz_array)


				i=0
				while i < len(mz_array):
				    mass_lower[i] = mz_array[i] - 0.0001
				    i+=1

				mass_upper = [0]*len(mz_array)
				i=0
				while i < len(mz_array):
				    mass_upper[i] = mz_array[i] + 0.0001
				    i+=1

				data1 = {'mz': mz_array,
				        'intensity': intensity_array
				        }
				data2 = {'mz': mass_upper,
				        'intensity': intensity_zeros
				        }

				data3 = {'mz': mass_lower,
				        'intensity': intensity_zeros
				        }

				df_1 = pandas.DataFrame (data1, columns = [ 'mz','intensity'])

				df_2 = pandas.DataFrame (data2, columns = ['mz','intensity'])

				df_3 = pandas.DataFrame (data3, columns = ['mz','intensity'])


				df_combine = [df_1,df_2,df_3]

				df_combine = pandas.concat(df_combine)

				df_combine.sort_values('mz', inplace=True)


				mz_array = df_combine['mz'].to_numpy()
				intensity_array = df_combine['intensity'].to_numpy()
				for analyte_id in df2.analyte_id.unique():
					ms1_data.append([df2[df2["analyte_id"] == analyte_id]["intensity"].to_numpy(),
				                     df2[df2["analyte_id"] == analyte_id]["mz"].to_numpy()])
					analyte_id_array.append(analyte_id)

				

				analyte_colour_array = [0]*len(analyte_id_array)

				analyte_id_array = [0 if str(z) == 'nan' else z for z in analyte_id_array]
				i = 0
				while i < len(analyte_id_array):
					analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
					i+=1



				colour=0
				self.ms1_plot.clear()
				for analyte in ms1_data:
					np.round(mz_array[0],2)
					np.round(intensity_array[0],2)
					if float(comboAnalyte) == analyte_id_array[colour]:
						print("comboAnalyte = colour")
						analyte[1] = [0 if str(i) == 'nan' else i for i in analyte[1]]
						analyte[0] = [0 if str(i) == 'nan' else i for i in analyte[0]]
						pen = pg.mkPen(color=(0, 0, 0))
						self.ms1_plot.plot( x=mz_array, y=intensity_array,pen=pen)
						print('INTENSITY ARRAY',intensity_array[0])
						print('MZ ARRAY',mz_array[0])

						i = 0
						while i < len(mz_array):
							if intensity_array[i] > 0:

								text = pg.TextItem(str(np.round(mz_array[i],2)),color=(0,0,0))
								self.ms1_plot.addItem(text)


								text.setPos(mz_array[i], intensity_array[i])
								i+=1
							else:
								i+=1
					else:

						print("comboAnalyte not = colour")
					colour +=1




		if chooseAnalyte == "Experiment Analyte ID":
			self.ms1_plot.setXRange(min=0,max=1200)
			self.ms1_plot.setYRange(min=5, max=100)
			ms1_data = []
			analyte_id_array = []

			if comboAnalyte != "Show All":
				print("Plot")
				for analyte_id in df3.experiment_analyte_id.unique():
					ms1_data.append([df3[df3["experiment_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
				                     df3[df3["experiment_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
					analyte_id_array.append(analyte_id)
				analyte_colour_array = [0]*len(analyte_id_array)
				print('analyte_id_array',analyte_id_array)

				i = 0
				while i < len(analyte_id_array):
					analyte_colour_array[i] = int((analyte_id_array[i] + 1) % 299)
					i+=1



				colour=0
				self.ms1_plot.clear()
				for analyte in ms1_data:

					if float(comboAnalyte) == analyte_id_array[colour]:
						print("comboAnalyte = colour")
						pen = pg.mkPen(color=(0, 0, 0))
						self.ms1_plot.plot( x=analyte[1], y=analyte[0],pen=pen)

						i = 0
						while i < len(analyte[1]):
							if analyte[0][i] > 0:

								text = pg.TextItem(str(np.round(analyte[1][i],2)),color=(0,0,0))
								self.ms1_plot.addItem(text)


								text.setPos(analyte[1][i], analyte[0][i])
								i+=1
							else:
								i+=1
					else:

						print("comboAnalyte not = colour")
					colour +=1

			# self.ms1_plot.enableAutoRange(axis=None)
		return


	def reset_all(self):
			state = State(self.Nullbutton.setChecked(False), self.Blankbutton.setChecked(False), 0, 1200, 0, 7)
			self.Blankbutton.setStyleSheet("background-color : lightgreen") 
			self.Blankbutton.setText('On')
			self.Nullbutton.setStyleSheet("background-color : lightgreen") 
			self.Nullbutton.setText('On')
			self.massCheckbox.setChecked(False)

			self.mMinBox.setText('0')
			self.mMaxBox.setText('1200')
			self.rtMinBox.setText('0')
			self.rtMaxBox.setText('7')



	def reset_massRange(self):

		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, self.rtMinBox.text(), self.rtMaxBox.text())
		self.mMinBox.setText('0')
		self.mMaxBox.setText('1200')

	def reset_rtRange(self):
		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7)
		self.rtMinBox.setText('0')
		self.rtMaxBox.setText('7')


	def Sample_Plot_DF(self):
		chooseAnalyte = self.chooseAnalyte.currentText()
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		if Sample_State == False:
			if self.massCheckbox.isChecked() == False:
				# if self.rtCheckbox.isChecked() == False:
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, 0, 7)

			#if self.rtCheckbox.isChecked() == False:
			#	state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7)

			if self.massCheckbox.isChecked() == True:
				# if self.rtCheckbox.isChecked() == False:
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())



			# if self.rtCheckbox.isChecked() == True:
			# 	if self.massCheckbox.isChecked() == False:
			# 		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())




			BlankState = state.return_blank_state()
			NullState = state.return_null_state()
			comboText=self.comboBox.currentText()

			############### Fetch data for mz_plot and plot sample from comboBox ###############

			total_df = import_dataframe(input_structure, input_type, comboText)



			if chooseAnalyte == "Sample Analyte ID":
				if NullState == True:
					total_df = total_df.dropna(subset=['analyte_id'])
					total_df=total_df[total_df.analyte_id != None]
					print('df Null Off', total_df.analyte_id)
				else:
					pass
				if BlankState == True:
					total_df=total_df[total_df.blank_analyte != True]
				else:
					pass
				total_df=total_df[total_df.replicate != 2]
				total_df=total_df[total_df.replicate != 3]

				total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
				total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

				sorted_df_mz = total_mz_df.sort_values(by=['mz'])


				total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

				total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
				sorted_df_rt = total_rt_df.sort_values(by=['rt'])
		




			if chooseAnalyte == "Replicate Analyte ID":
				if NullState == True:
					total_df = total_df.dropna(subset=['replicate_analyte_id'])
					total_df=total_df[total_df.replicate_analyte_id != None]

				else:
					pass
				if BlankState == True:
					total_df=total_df[total_df.blank_analyte != True]
				else:
					pass
				total_df=total_df[total_df.replicate != 2]
				total_df=total_df[total_df.replicate != 3]

				total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
				total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

				sorted_df_mz = total_mz_df.sort_values(by=['mz'])


				total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

				total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
				sorted_df_rt = total_rt_df.sort_values(by=['rt'])



			if chooseAnalyte == "Experiment Analyte ID":
				if NullState == True:
					total_df = total_df.dropna(subset=['experiment_analyte_id'])
					total_df=total_df[total_df.experiment_analyte_id != None]

				else:
					pass
				if BlankState == True:
					total_df=total_df[total_df.blank_analyte != True]
				else:
					pass
				total_df=total_df[total_df.replicate != 2]
				total_df=total_df[total_df.replicate != 3]

				total_mz_df = total_df[total_df.rt > float(state.return_rtMin_state())]
				total_mz_df = total_mz_df[total_mz_df.rt < float(state.return_rtMax_state())]

				sorted_df_mz = total_mz_df.sort_values(by=['mz'])


				total_rt_df = total_df[total_df.mz > float(state.return_massMin_state())]

				total_rt_df = total_rt_df[total_rt_df.mz < float(state.return_massMax_state())]
				sorted_df_rt = total_rt_df.sort_values(by=['rt'])
			# else:
			# 	pass
		else:
			sorted_df_mz = []
			sorted_df_rt = []


		return sorted_df_mz, sorted_df_rt



	def Replicate_Plot_DF(self):
		chooseAnalyte = self.chooseAnalyte.currentText()
		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()

		if Replicate_State == False:
			if self.massCheckbox.isChecked() == False:
				# if self.rtCheckbox.isChecked() == False:
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), 0, 1200, 0, 7)

			#if self.rtCheckbox.isChecked() == False:
			#	state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), 0, 7)

			if self.massCheckbox.isChecked() == True:
				# if self.rtCheckbox.isChecked() == False:
				state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())



			# if self.rtCheckbox.isChecked() == True:
			# 	if self.massCheckbox.isChecked() == False:
			# 		state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())




			BlankState = state.return_blank_state()
			NullState = state.return_null_state()
			comboText=self.comboBox.currentText()

			############### Fetch data for mz_plot and plot sample from comboBox ###############

			total_df = import_dataframe(input_structure, input_type, comboText)

			if chooseAnalyte == "Sample Analyte ID":
				if NullState == True:
					total_df = total_df.dropna(subset=['analyte_id'])
					total_df=total_df[total_df.analyte_id != None]
				else:
					pass
				if BlankState == True:
					total_df=total_df[total_df.blank_analyte != True]
				else:
					pass
				total_df_r1=total_df[total_df.replicate != 2]
				total_df_r1=total_df_r1[total_df_r1.replicate != 3]

				total_mz_df_r1 = total_df_r1[total_df_r1.rt > float(state.return_rtMin_state())]
				total_mz_df_r1 = total_mz_df_r1[total_mz_df_r1.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r1 = total_mz_df_r1.sort_values(by=['mz'])


				total_rt_df_r1 = total_df_r1[total_df_r1.mz > float(state.return_massMin_state())]

				total_rt_df_r1 = total_rt_df_r1[total_rt_df_r1.mz < float(state.return_massMax_state())]
				sorted_df_rt_r1 = total_rt_df_r1.sort_values(by=['rt'])

				total_df_r2=total_df[total_df.replicate != 1]
				total_df_r2=total_df_r2[total_df_r2.replicate != 3]

				total_mz_df_r2 = total_df_r2[total_df_r2.rt > float(state.return_rtMin_state())]
				total_mz_df_r2 = total_mz_df_r2[total_mz_df_r2.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r2 = total_mz_df_r2.sort_values(by=['mz'])


				total_rt_df_r2 = total_df_r2[total_df_r2.mz > float(state.return_massMin_state())]

				total_rt_df_r2 = total_rt_df_r2[total_rt_df_r2.mz < float(state.return_massMax_state())]
				sorted_df_rt_r2 = total_rt_df_r2.sort_values(by=['rt'])




				total_df_r3=total_df[total_df.replicate != 1]
				total_df_r3=total_df_r3[total_df_r3.replicate != 2]

				total_mz_df_r3 = total_df_r3[total_df_r3.rt > float(state.return_rtMin_state())]
				total_mz_df_r3 = total_mz_df_r3[total_mz_df_r3.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r3 = total_mz_df_r3.sort_values(by=['mz'])


				total_rt_df_r3 = total_df_r3[total_df_r3.mz > float(state.return_massMin_state())]

				total_rt_df_r3 = total_rt_df_r3[total_rt_df_r3.mz < float(state.return_massMax_state())]
				sorted_df_rt_r3 = total_rt_df_r3.sort_values(by=['rt'])
				#sorted_df_rt = total_df.sort_values(by=['rt'])

			if chooseAnalyte == "Replicate Analyte ID":
				if NullState == True:
					total_df = total_df.dropna(subset=['replicate_analyte_id'])
					total_df=total_df[total_df.replicate_analyte_id != None]
				else:
					pass
				if BlankState == True:
					total_df=total_df[total_df.blank_analyte != True]
				else:
					pass
				total_df_r1=total_df[total_df.replicate != 2]
				total_df_r1=total_df_r1[total_df_r1.replicate != 3]

				total_mz_df_r1 = total_df_r1[total_df_r1.rt > float(state.return_rtMin_state())]
				total_mz_df_r1 = total_mz_df_r1[total_mz_df_r1.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r1 = total_mz_df_r1.sort_values(by=['mz'])


				total_rt_df_r1 = total_df_r1[total_df_r1.mz > float(state.return_massMin_state())]

				total_rt_df_r1 = total_rt_df_r1[total_rt_df_r1.mz < float(state.return_massMax_state())]
				sorted_df_rt_r1 = total_rt_df_r1.sort_values(by=['rt'])

				total_df_r2=total_df[total_df.replicate != 1]
				total_df_r2=total_df_r2[total_df_r2.replicate != 3]

				total_mz_df_r2 = total_df_r2[total_df_r2.rt > float(state.return_rtMin_state())]
				total_mz_df_r2 = total_mz_df_r2[total_mz_df_r2.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r2 = total_mz_df_r2.sort_values(by=['mz'])


				total_rt_df_r2 = total_df_r2[total_df_r2.mz > float(state.return_massMin_state())]

				total_rt_df_r2 = total_rt_df_r2[total_rt_df_r2.mz < float(state.return_massMax_state())]
				sorted_df_rt_r2 = total_rt_df_r2.sort_values(by=['rt'])




				total_df_r3=total_df[total_df.replicate != 1]
				total_df_r3=total_df_r3[total_df_r3.replicate != 2]

				total_mz_df_r3 = total_df_r3[total_df_r3.rt > float(state.return_rtMin_state())]
				total_mz_df_r3 = total_mz_df_r3[total_mz_df_r3.rt < float(state.return_rtMax_state())]



				sorted_df_mz_r3 = total_mz_df_r3.sort_values(by=['mz'])


				total_rt_df_r3 = total_df_r3[total_df_r3.mz > float(state.return_massMin_state())]

				total_rt_df_r3 = total_rt_df_r3[total_rt_df_r3.mz < float(state.return_massMax_state())]
				sorted_df_rt_r3 = total_rt_df_r3.sort_values(by=['rt'])

				#sorted_df_rt = total_df.sort_values(by=['rt'])
		else:
			sorted_df_mz_r1 = []
			sorted_df_mz_r2 = []
			sorted_df_mz_r3 = []
			sorted_df_rt_r1 = []
			sorted_df_rt_r2 = []
			sorted_df_rt_r3 = []

		return sorted_df_mz_r1, sorted_df_rt_r1, sorted_df_mz_r2, sorted_df_rt_r2, sorted_df_mz_r3, sorted_df_rt_r3
################### Connecting Region 1 to Region 3 ###########################
	def Experiment_Plot_DF(self):

		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		print('EXPERIMENT PLOT DF 1')

		if Experiment_State == False:
			print('EXPERIMENT PLOT DF 2')
			state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())

			BlankState = state.return_blank_state()
			NullState = state.return_null_state()

			experiment_df = import_experiment_dataframe(input_structure)
			experiment_df =experiment_df.sort_values(by='sample_processing_order')
			
			if BlankState == True:
				experiment_df=experiment_df[experiment_df.experiment_analyte_is_blank != True]
			else:
				pass
			sample_name_df = experiment_df.sample_name
			sample_name_sorted_df = sample_name_df.drop_duplicates()

			samples = sample_name_sorted_df.values
		else:
			samples=[]
			experiment_df = []

		return samples, experiment_df




	def Diversity_Plot_DF(self):

		tab_state = Tab_State(self.btn2.isEnabled(), self.btn3.isEnabled(), self.btn4.isEnabled(), self.btn5.isEnabled())
		Sample_State = tab_state.return_sample_state()
		Replicate_State = tab_state.return_replicate_state()
		Experiment_State = tab_state.return_experiment_state()
		Diversity_State = tab_state.return_diversity_state()

		if Diversity_State == False:

			state = State(self.Nullbutton.isChecked(), self.Blankbutton.isChecked(), self.mMinBox.text(), self.mMaxBox.text(), self.rtMinBox.text(), self.rtMaxBox.text())

			BlankState = state.return_blank_state()
			NullState = state.return_null_state()

			experiment_df = import_experiment_dataframe(input_structure)
			experiment_df =experiment_df.sort_values(by='sample_processing_order')
			
			if BlankState == True:
				experiment_df=experiment_df[experiment_df.experiment_analyte_is_blank != True]
			else:
				pass
			sample_name_df = experiment_df.sample_name
			sample_name_sorted_df = sample_name_df.drop_duplicates()
			sample_name_sorted_df = sample_name_sorted_df.iloc[::-1]

			samples = sample_name_sorted_df.values

		else:
			samples=[]
			experiment_df = []

		return samples, experiment_df

	def ms1_Plot_DF(self):
		comboText=self.comboBox.currentText()

		df1, df2, df3 = import_ms1_dataframe(input_structure, input_type, comboText)
		df4 = import_dataframe(input_structure, input_type, comboText)
		df1.sort_values('average_mass', inplace=True)
		df2.sort_values('max_peak_intensity_mass', inplace=True)
		df3.sort_values('average_mass', inplace=True)





		df=df4[df4.replicate == 1]

		# df = df.sort_values('mz', inplace= True)
		# sorted_df4.groupby('scan').intensity.max()
		
		
		# print('DF MAX', sorted_df4)





		# ms1_data = []
		# for analyte_id in sorted_df4.analyte_id.unique():
		# 	for max_value in sorted_df4.scan.unique():
		# 		sorted_max_df4 = sorted_df4.loc[sorted_df4['intensity'].idxmax()]
		# 		print(sorted_max_df4)







		# 		for analyte_id in df3.experiment_analyte_id.unique():
		# 			ms1_data.append([df3[df3["experiment_analyte_id"] == analyte_id]["relative_intensity"].to_numpy(),
		# 		                     df3[df3["experiment_analyte_id"] == analyte_id]["average_mass"].to_numpy()])
		# 			analyte_id_array.append(analyte_id)
		return df1, df, df3

	def nullButton(self):
	    # method called by button 

		# if button is checked 
		if self.Nullbutton.isChecked():
			#self.update_tab()
			# setting background color to light-blue 
			self.Nullbutton.setStyleSheet("background-color : lightgrey") 
			self.Nullbutton.setText('Off') 

		# if it is unchecked 
		else: 
			#self.update_tab()
			# set background color back to light-grey 
			self.Nullbutton.setStyleSheet("background-color : lightgreen") 
			self.Nullbutton.setText('On')


		
		



	def blankButton(self):

	    # method called by button 

	  
	    # if button is checked 
		if self.Blankbutton.isChecked():
			#self.update_tab()
	        # setting background color to light-blue 
			self.Blankbutton.setStyleSheet("background-color : lightgrey") 
			self.Blankbutton.setText('Off') 

	    # if it is unchecked 
		else: 
			#self.update_tab()
			# set background color back to light-grey 
			self.Blankbutton.setStyleSheet("background-color : lightgreen") 
			self.Blankbutton.setText('On')


	def mMinLimit(self, text):
		print(self.mMaxBox.text())
		maxvalue=float(self.mMaxBox.text())
		self.mz_plot.setXRange(min=float(text),max=maxvalue)
		self.mz_r1_plot.setXRange(min=float(text),max=maxvalue)
		self.mz_r2_plot.setXRange(min=float(text),max=maxvalue)
		self.mz_r3_plot.setXRange(min=float(text),max=maxvalue)

	def mMaxLimit(self, text):
		print(self.mMinBox.text())
		minvalue=float(self.mMinBox.text())
		self.mz_plot.setXRange(min=minvalue,max=float(text))
		self.mz_r1_plot.setXRange(min=minvalue,max=float(text))
		self.mz_r2_plot.setXRange(min=minvalue,max=float(text))
		self.mz_r3_plot.setXRange(min=minvalue,max=float(text))




	def rtMinLimit(self, text):
		print(self.rtMaxBox.text())
		maxvalue=float(self.rtMaxBox.text())
		self.rt_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r1_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r2_plot.setXRange(min=float(text),max=maxvalue)
		self.rt_r3_plot.setXRange(min=float(text),max=maxvalue)


	def rtMaxLimit(self, text):
		print(self.rtMinBox.text())
		minvalue=float(self.rtMinBox.text())
		self.rt_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r1_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r2_plot.setXRange(min=minvalue,max=float(text))
		self.rt_r3_plot.setXRange(min=minvalue,max=float(text))


	def hide_sample(self):

		#self.btn2.setEnabled(True)
		self.mz_plot.hide()
		self.rt_plot.hide()
		#self.ms1_plot.hide()
		#self.ms2_plot.hide()

	def show_sample(self):

		tab_state = Tab_State(self.btn2.setEnabled(False), self.btn3.setEnabled(True), self.btn4.setEnabled(True), self.btn5.setEnabled(True))
		self.btn2.setEnabled(False)
		self.btn4.setEnabled(True)
		self.btn5.setEnabled(True)
		self.btn3.setEnabled(True)
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)
		#self.comboBox.setEnabled(True)
		self.ms1_plot.show()
		self.ms2_plot.show()


		self.mz_plot.show()
		self.rt_plot.show()
			#self.ms1_plot.show()
			#self.ms2_plot.show()

		
		
	def show_replicate(self):
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(False), self.btn4.setEnabled(True), self.btn5.setEnabled(True))
		self.btn3.setEnabled(False)
		self.btn4.setEnabled(True)
		self.btn2.setEnabled(True)
		self.btn5.setEnabled(True)
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)

		self.mz_r1_plot.show()
		self.mz_r2_plot.show()
		self.mz_r3_plot.show()
		self.rt_r1_plot.show()
		self.rt_r2_plot.show()
		self.rt_r3_plot.show()
		self.ms1_plot.show()
		self.ms2_plot.show()



	def hide_replicate(self):
		self.mz_r1_plot.hide()
		self.mz_r2_plot.hide()
		self.mz_r3_plot.hide()
		self.rt_r1_plot.hide()
		self.rt_r2_plot.hide()
		self.rt_r3_plot.hide()

	def show_experiment(self):
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(True), self.btn4.setEnabled(False), self.btn5.setEnabled(True))
		self.btn4.setEnabled(False)
		self.btn2.setEnabled(True)
		self.btn3.setEnabled(True)
		self.btn5.setEnabled(True)
		# self.Nullbutton.setEnabled(False)
		self.ex_bOn_plot.show()
		# self.chooseSample.setEnabled(False)
		# self.showNull.setEnabled(False)
		self.ms1_plot.show()
		self.ms2_plot.show()
		# if self.Blankbutton.isChecked() == False:
		# 	self.ex_bOn_plot.show()
		# 	#self.ex_bOff_plot.hide()
		# if self.Blankbutton.isChecked()  == True:
		# 	self.ex_bOn_plot.hide()
			#self.ex_bOff_plot.show()


	def hide_experiment(self):

		self.ex_bOn_plot.hide()
		#self.ex_bOff_plot.hide()

	def show_diversity(self):
		tab_state = Tab_State(self.btn2.setEnabled(True), self.btn3.setEnabled(True), self.btn4.setEnabled(True), self.btn5.setEnabled(False))
		self.btn4.setEnabled(True)
		self.btn2.setEnabled(True)
		self.btn3.setEnabled(True)
		self.btn5.setEnabled(False)
		# self.Nullbutton.setEnabled(True)
		# self.chooseSample.setEnabled(True)
		# self.showNull.setEnabled(True)
		self.div_plot.show()
	def hide_diversity(self):
		self.div_plot.hide()
		pass



	def close_application(self):
		choice = QMessageBox.question(self, 'Quit',
											"Are you sure you want to quit?",
											QMessageBox.Yes | QMessageBox.No)
		if choice == QMessageBox.Yes:
			print("Closing Application")
			sys.exit()
		else:
			pass
			
		

		
def run():

	app = QApplication(sys.argv)
	app.setStyle('Fusion')
	GUI = Window()
	sys.exit(app.exec_())
run()
