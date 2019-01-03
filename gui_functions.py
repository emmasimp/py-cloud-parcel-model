import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from scipy.optimize import curve_fit
import math

def getTempVar(self, filename):
	try:
		Temperature = pd.read_csv(filename, skiprows=2, sep='\t')
		Temperature.columns = ['Date','Time','1','2','3','4','5','6','7','8','unknown'] # name columnsT1 (C)	T2 (C)	T3 (C)	T4 (C)	T5 (C)	T6 (C)	T7 (C)	T8 (C)	
		Temperature['datetime'] = Temperature['Date'] + ' ' + Temperature['Time']  
		Temperature['datetime'] = [dt.datetime.strptime(x, '%d/%m/%Y %H:%M:%S') for x in Temperature['datetime']]     # create datetime column
	except:
		Temperature = pd.read_csv(filename, skiprows=1, sep=',')
		Temperature['datetime'] = Temperature['Date'] + ' ' + Temperature[' Time']       # create datetime column
		Temperature['datetime'] = [dt.datetime.strptime(x, '%d/%m/%y %H:%M:%S') for x in Temperature['datetime']]
	
	Temperature.index = Temperature['datetime']

	# find thermocouple with lowest temperature
	try :
		minT = Temperature[['1','2','3','4','5','6','7','8']].sum().idxmin() # minT is an index to the therocouple
	except:
		minT = Temperature[[' T1 (C)',' T2 (C)',' T3 (C)',' T4 (C)',' T5 (C)',' T6 (C)',' T7 (C)',' T8 (C)']].sum().idxmin() # minT is an index to the therocouple

	self.Temperature_data = pd.DataFrame() # dataframes
	self.Temperature_data['Temperature'] = Temperature[minT] + 273.15      # use lowest thermocouple data (minT)

	return self

def getPressVar(self, filename):
	Pressure = pd.read_csv(filename, skiprows=1, sep='\t')
	Pressure['datetime'] = Pressure['Date'] + ' ' + Pressure['Time']       # create datetime column

	try:
		Pressure['datetime'] = [dt.datetime.strptime(x, '%d/%m/%Y %H:%M:%S') for x in Pressure['datetime']]
	except:
		Pressure['datetime'] = [dt.datetime.strptime(x, '%d/%m/%y %H:%M:%S') for x in Pressure['datetime']]

	self.Pressure_data = pd.DataFrame() 
	self.Pressure_data['Pressure'] = Pressure['Pressure (mBar)']
	self.Pressure_data.index = Pressure['datetime']

	return self

def clearlayout(self, layout):
		if layout is not None:
			while layout.count():
				item = layout.takeAt(0)
				widget = item.widget()
				if widget is not None:
					widget.deleteLater()
				else:
					self.clearLayout(item.layout())

def getData(self):
	try:
		f.getTempVar(self, self.temp_filename)
	except:
		self.box = QtGui.QMessageBox()
		self.box.setText('Please select temperature data file')
		self.box.exec_()
	try:
		f.getPressVar(self, self.press_filename)
	except:
		self.box = QtGui.QMessageBox()
		self.box.setText('Please select pressure data file')
		self.box.exec_()

	return self


def getCDPVar(self,filename):
	CDP = pd.DataFrame()
	df = pd.read_csv(filename, skiprows=2, sep=r"\s+",names=range(54)) # more than one delimiter (\s+)
    # label each column with numbers 0 - 54
	df.columns = df.iloc[0]                                            # get column names from first row
	df.drop(df.index[[0]], inplace=True)                               # delete first row
	df['datetime'] = df['Date'] + ' ' + df['Time']                     # create datetime column
	df['datetime'] = [dt.datetime.strptime(x, '%d/%m/%Y %H:%M:%S') for x in df['datetime']]
	df.index = df['datetime']                                          # use datetime as index
	CDP = CDP.append(df)  
	CDP = CDP.sort_index()
	# replace NaNs with zeros
	CDP.fillna(value=0, inplace=True)
	# find indices to columns containing concentrations of each size bin
	bin_first = CDP.columns.get_loc(3.0)
	bin_last = CDP.columns.get_loc(50.0)
	# Create array of binned data, CDP1
	CDP1 = CDP.iloc[:,bin_first:bin_last+1].sum(axis=1)
	# Create variable for sum of sizes above 12microns, ie ice
	bin_12 = CDP.columns.get_loc(12.0)
	bin_5 = CDP.columns.get_loc(5.0)
	CDP12 = CDP.iloc[:,bin_12:bin_last+1].sum(axis=1)
	CDP5 = CDP.iloc[:,bin_5:bin_last+1].sum(axis=1)
	
	#self.CDP_data = pd.DataFrame() 
	#self.CDP_data['Data'] = CDP1
	#self.CDP_data.index = CDP['datetime']

	#self.CDP_data = CDP1
	return CDP1, CDP12