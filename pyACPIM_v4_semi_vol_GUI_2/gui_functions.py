import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from scipy.optimize import curve_fit
import math
from scipy.signal import find_peaks


def getTempVar(self, filename):
	try:
		Temperature = pd.read_csv(filename, skiprows=2, sep='\t')
		Temperature.columns = ['Date','Time','1','2','3','4','5','6','7','8','unknown'] # name columnsT1 (C)    T2 (C)  T3 (C)  T4 (C)  T5 (C)  T6 (C)  T7 (C)  T8 (C)  
		Temperature['datetime'] = Temperature['Date'] + ' ' + Temperature['Time']  
		Temperature['datetime'] = [dt.datetime.strptime(x, '%d/%m/%Y %H:%M:%S') for x in Temperature['datetime']]     # create datetime column
	except:
		Temperature = pd.read_csv(filename, skiprows=1, sep=',')
		Temperature.columns = ['date','time','lab','1','2','3','4','5','6','7','8',' '] # name columns
		Temperature['datetime'] = Temperature['date'] + ' ' + Temperature['time']       # create datetime column
		Temperature['datetime'] = [dt.datetime.strptime(x, '%d/%m/%y %H:%M:%S') for x in Temperature['datetime']]
	
	Temperature.index = Temperature['datetime']

	# find thermocouple with lowest temperature
	try :
		minT = Temperature[['1','2','3','4','5','6','7','8']].sum().idxmin() # minT is an index to the therocouple
	except:
		minT = Temperature[[' T1 (C)',' T2 (C)',' T3 (C)',' T4 (C)',' T5 (C)',' T6 (C)',' T7 (C)',' T8 (C)']].sum().idxmin() # minT is an index to the therocouple

	self.Temperature_data = pd.DataFrame() # dataframes
	self.Temperature_data['Temperature'] = Temperature[minT] + 273.15      # use lowest thermocouple data (minT)
	self.temp_data = Temperature # for calculating least squares fit

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

def getFSSPVar(self,filename):
	FSSP = pd.read_csv(filename,sep='\t')
	return FSSP


def firstNonNan(listfloats):
	""" function that returns first non nan value in list 
	"""
	for item in listfloats:
		if math.isnan(item) == False:
			return item

def Temp_correct(Press, Temp, time_constant):
	""" Function that calculates adiabatic temperaure change from pressure 
		and calculates corrected temperature data using a time consant
		
		Press is a list of pressure values
		
		Temp is a list of temperature values in Kelvin, its helpful to smooth the data before
		passing to this function
	"""
	
	gamma = 287.058/1005.0
	P0 = Press[0]
	T0 = firstNonNan(Temp)
	alpha = (1/time_constant)
	
	adiabatic = np.zeros(len(Press)) # create empty array for calculation of adiabatic T
	Delta = np.zeros_like(adiabatic)
	Delta[0] = T0

	for i in range(0,len(Press)):
		adiabatic[i] = T0*(Press[i]/P0)**(gamma)
		if i == 0:
			Delta[i] = Temp[i]
		else:
			Delta[i] = Temp[i]-Temp[i-1]
			
	T_corr = Delta/alpha + Temp
	T_corr[0] = Temp[0]
	
	return T_corr, adiabatic            

def temp_function(x,T1,T2,T3):
	return T1*np.exp(-T2*x)+T3

def calc_fit_parameters(self,Temperature,Pressure):
	new = pd.DataFrame()
	# start/end times updated with QLineEdits
	start = str(self.start_times[int(self.select_expansion_number.text())-1])[:10]+' '+self.simulation_start_time.text()
	end = str(self.end_times[(int(self.select_expansion_number.text())-1)])[:10]+' '+self.simulation_end_time.text()
	# calculate mean of all 8 thremocouples
	for col in ['1','2','3','4','5','6','7','8']:
		T = Temperature[col][start:end].rolling(2).mean()+ 273.15
		P = Pressure['Pressure'][start:end]
		T = T.resample('1S').pad()
		P = P.resample('1S').pad()
		T0 = firstNonNan(T)
		time_c = 60 # make this an input at some point
		temp_corr,adiabatic = Temp_correct(P,T, time_c)
		new['Temp_corr'+col] = temp_corr
	
	new = new.assign(mean=new.mean(axis=1))
	new = new.assign(max_T=new.max(axis=1))
	new = new.assign(min_T=new.min(axis=1))
	new = new.assign(stdev=new.std(axis=1))
	new = new[2:]
	# calculate fits to mean, mean+stdev, mean-stdev
	x = np.arange(len(new))[:]
	y = new['mean'][:].values
	error = new['stdev'][:]
	y2 = y + error
	y3 = y - error

	y[pd.isnull(y)] = 0.0
	y2[pd.isnull(y2)] = 0.0
	y3[pd.isnull(y3)] = 0.0

	y = pd.to_numeric(y, errors='coerce')
	y2 = pd.to_numeric(y2, errors='coerce')
	y3 = pd.to_numeric(y3, errors='coerce')

	popt, pcov = curve_fit(temp_function, x, y, p0=[3.5,0.017,288],bounds=([1e-5,1e-10,200],[20,1,1e3]))
	popt2, pcov2 = curve_fit(temp_function, x, y2, p0=[3.5,0.017,288],bounds=([1e-5,1e-10,200],[20,1,1e3]))
	popt3, pcov3 = curve_fit(temp_function, x, y3, p0=[3.5,0.017,288],bounds=([1e-10,1e-10,200],[20,1,1e3]))

	self.mean_temp_data = y
	self.min_temp_data = new.min_T
	self.max_temp_data = new.max_T
	self.mean_temp_fit = popt
	self.lower_temp_fit = popt3
	self.upper_temp_fit = popt2

	# pressure data fit
	x = np.arange(len(P))[:]
	popt, pcov = curve_fit(temp_function,x,P,p0=[628,2.6e-3,1000])
	self.press_fit = popt
	self.press_data = P.values

	return self

def find_expansions(self):
	P = self.Pressure_data['Pressure']
	# subsample data to smooth
	P3 = (P[::5]-1000)*-1
	# replace nans with 0
	P3[np.isnan(P3)]=0
	# find the most prominant peaks in data to signal an expanion end
	peaks, _ = find_peaks(P3, prominence=(100, None))
	# save times of end of expansions
	end_times = P3.index[peaks].tolist()
	number_of_expansions = len(end_times)
	# plot data with peaks identified
   # plt.plot(P3.values)
   # plt.plot(peaks,P3[peaks].values,"x")
	# find start times by finding the last identified peak in P[end_time[1]:end_time[2]]
	start_times = []
	for i in np.arange(number_of_expansions):
		if i == 0:
			P5 = P[:end_times[1]]
		elif i < number_of_expansions-1:
			P5 = P[end_times[i]:end_times[i+1]]
		else:
			P5 = P[end_times[i]:]
		peaks, _ = find_peaks(P5)#, prominence=(None,0.01 ),distance=5,width=2) 
		#start_times.append(str(P5.index[peaks[-1]]))
		start_times.append(P5.index[P5.diff(periods=3)==min(P5.diff(periods=3).replace(np.nan,0))].strftime('%Y-%m-%d %H:%M:%S')[0])
	
	start_times = [dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')-dt.timedelta(seconds=5) for x in start_times]
#	start_times = [dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in start_times]

	# make sure start time is < end time
	new_start = []
	new_end = []

	for start,end in zip(start_times[:-1], end_times[1:]):
		if end > start:
			new_start.append(start)
			new_end.append(end)

	self.number_of_expansions = len(new_start)
	self.start_times = new_start
	self.end_times = new_end

	return self