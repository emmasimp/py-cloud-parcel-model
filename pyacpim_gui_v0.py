import sys
from PySide import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import gui_functions as f
import pickle
import numpy as np
from threading import Thread
from netCDF4 import Dataset
import main as model

qt_app = QtGui.QApplication(sys.argv)
finished = QtCore.Signal()

class main_window(QtGui.QWidget):
	

	def __init__(self):
		QtGui.QWidget.__init__(self)

		palette = QtGui.QPalette()
		palette.setColor(QtGui.QPalette.Background,QtGui.QColor('#4c4c4c'))
		self.setPalette(palette)
		# over all layout is three vertical columns (3 horizontal boxes)
		self.layout = QtGui.QHBoxLayout()
		# each column in a series of vertical boxes
		self.column_layout1 = QtGui.QVBoxLayout()
		self.column_layout2 = QtGui.QVBoxLayout()
		self.column_layout3 = QtGui.QVBoxLayout()
		self.simulation_type_layout = QtGui.QVBoxLayout()
		self.parcel_layout = QtGui.QHBoxLayout()
		self.chamber_layout = QtGui.QHBoxLayout()
		self.chamber_layout_temp = QtGui.QHBoxLayout()
		self.chamber_layout_press = QtGui.QHBoxLayout()
		self.grid = QtGui.QGridLayout()
		self.grid_comp = QtGui.QGridLayout()

		self.tab_widget = QtGui.QTabWidget()
		#self.tab_widget.setStyleSheet('background-color:blue')
		self.tab1 = QtGui.QWidget()
		self.tab_widget.setPalette(palette)
		self.tab_widget.addTab(self.tab1,'Initial Conditions')

		self.tab2 = QtGui.QWidget()
		self.tab_widget.addTab(self.tab2,'Aerosol Properties')

		self.tab3 = QtGui.QWidget()
		self.tab_widget.addTab(self.tab3,'Advanced Settings')

		self.setStyleSheet("""QPushButton, QComboBox{
				border: 2px solid gray;
				border-radius: 6px;
				background-color: #4c4c4c;
				color: #d1cfcf;
				min-width: 150px;
				max-width: 200px;
				outline: 0;
				}
				
				QListView{
				outline:0;
				}
				QComboBox { 
				combobox-popup: 0;
				color: #d1cfcf;
				padding: 1px 1px 1px 1px; 
				min-width: 250px
				}
				QComboBox QListView {
				background-color: #4c4c4c;
				color: #d1cfcf;
				min-width: 15em;
				}
				

				QComboBox::drop-down {
    			subcontrol-origin: padding;
    			subcontrol-position: top right;
    			width: 15px;
    			background-color: gray;
    			border-left-width: 1px;
    			border-left-color: #4c4c4c;
    			border-left-style: solid; /* just a single line */
    			border-top-right-radius: 3px; /* same radius as the QComboBox */
    			border-bottom-right-radius: 3px;
    			}

				QPushButton::hover{
				border: #4c4c4c;
				}

				QTabWidget::pane{
				border-top: 2px solid #4c4c4c;
				}
				QTabBar::tab{
				background:#4c4c4c;
				color: #d1cfcf;
                border: 2px solid gray;
                border-bottom-color: #4c4c4c;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
                min-width: 8px;
                padding: 2px;
                }
                
                QTabBar::tab:selected {
                border-color: #4c4c4c;
                border-bottom-color: #4c4c4c;

                }
                QTabBar::tab:!selected {
                margin-top: 3px ; /* make non-selected tabs look smaller */
                }
                QLabel{
                color: #d1cfcf;
                }
                QRadioButton, QCheckBox{
                background-color: #4c4c4c;
                color: #d1cfcf;
                }
                
                QCheckBox::indicator {
                border: 2px solid gray;
                background: none;

                }
                QCheckBox::indicator:checked{
                background: yellow;
                }

                QLineEdit {
                border: 2px solid gray;
                border-radius : 10px;
                padding: 0 8px;
                background: #4c4c4c;
                color: #d1cfcf;
                min-width: 80px;
                max-width: 100px;
                }
                QLineEdit#updraft:{
                max-width: 50px;
                }

                """)

        # Create push buttons
        # column 1
		self.chamber_button = QtGui.QRadioButton('&Chamber', self)
		self.chamber_button.toggled.connect(self.onChamberClicked)
		self.chamber_button.setChecked(True)
		self.parcel_button = QtGui.QRadioButton('&Parcel', self)
		self.parcel_button.toggled.connect(self.onParcelClicked)
		self.sim_type_label = QtGui.QLabel('Simulation type')
		self.add_aerosol_type_button = QtGui.QPushButton('&Add another aerosol type',self)
		self.add_aerosol_type_button.clicked.connect(self.add_aerosol_clicked)
		self.submit_button = QtGui.QPushButton('&Run',self)
		self.submit_button.clicked.connect(self.submit_clicked)
		self.savebutton = QtGui.QPushButton('&Select Output file', self)
		self.savebutton.clicked.connect(self.onSaveButtonClicked)
		self.runtime = QtGui.QLineEdit('270')
		self.runtime_label = QtGui.QLabel('Runtime')
		self.aerosol_composition_button = QtGui.QPushButton('&Aerosol Composition',self)
		self.aerosol_composition_button.clicked.connect(self.aerosol_comp_clicked)
		self.Dlow = QtGui.QLineEdit('10e-9')
		self.Dlow_label = QtGui.QLabel('Dlow')
		self.rkm = QtGui.QLineEdit('1')
		self.rkm_label = QtGui.QLabel('rkm')
		self.time_step = QtGui.QLineEdit('1')
		self.time_step_label = QtGui.QLabel('time step')
		self.alpha_crit = QtGui.QLineEdit('70')
		self.alpha_crit_label = QtGui.QLabel('Alpha crit')
		# column 3
		self.plot_button = QtGui.QPushButton('&Plot')
		self.plot_button.clicked.connect(self.plot)
		




		# input boxes
		self.number_modes = QtGui.QLineEdit('1')
		self.number_modes_enter = QtGui.QPushButton('&Enter',self)
		self.number_modes_enter.setChecked(True)
		self.number_modes_enter.clicked.connect(self.number_modes_enter_clicked)
		self.number_modes_label = QtGui.QLabel()
		self.number_modes_label.setText('number of modes')
		self.num_modes_layout = QtGui.QHBoxLayout()
		self.num_modes_layout.addWidget(self.number_modes_label)
		self.num_modes_layout.addWidget(self.number_modes)
		self.num_modes_layout.addWidget(self.number_modes_enter)

		self.RH_int = QtGui.QLineEdit('0.9')
		self.RH_label = QtGui.QLabel()
		self.P = QtGui.QLineEdit('950')
		self.P_label = QtGui.QLabel('Pressure (hPa)')
		self.RH_label.setText('RH')
		self.T = QtGui.QLineEdit('273')
		self.T_label = QtGui.QLabel('Temperature (K)')
		self.int_con_layout = QtGui.QVBoxLayout()
		self.int_con_layout.addWidget(self.RH_label)
		self.int_con_layout.addWidget(self.RH_int)
		self.int_con_layout.addWidget(self.P_label)
		self.int_con_layout.addWidget(self.P)
		self.int_con_layout.addWidget(self.T_label)
		self.int_con_layout.addWidget(self.T)
		self.freezing_criteria = QtGui.QComboBox(self)
		self.freezing_criteria.addItems(['heterogeneous freezing criteria','RH>1	[default]','Activated Drops', 'Threshold water mass'])
		# column 2
		# canvas for plotting model results
		self.figure = Figure()
		self.canvas = FigureCanvas(self.figure)
		self.toolbar = NavigationToolbar(self.canvas, self)
		# set the layout of figure/canvas
		self.layout_fig = QtGui.QVBoxLayout()
		self.layout_fig.addWidget(self.toolbar)
		self.layout_fig.addWidget(self.canvas)
		
		# column 3
	#	self.dialogbox = QtGui.QComboBox(self)
#		self.dialogbox.addItems(['ice_total','drop_total','temperature','RH'])
#		self.dialogbox2 = QtGui.QComboBox(self)
#		self.dialogbox2.addItems(['ice_total','drop_total','temperature','RH'])
		self.RH = QtGui.QCheckBox('RH',self)
		self.RH.setStyleSheet("background-color: #4c4c4c")
		self.RH.setChecked(True)
		self.temperature = QtGui.QCheckBox('Temperature',self)
		self.temperature.setStyleSheet("background-color: #4c4c4c")
		self.temperature.setChecked(True)
		self.drop_total = QtGui.QCheckBox('Number of drops',self)
		self.drop_total.setStyleSheet("background-color: #4c4c4c")
		self.ice_total = QtGui.QCheckBox('Number of ice crystals',self)
		self.ice_total.setStyleSheet("background-color: #4c4c4c")
		self.pressure= QtGui.QCheckBox('Pressure',self)
		self.pressure.setStyleSheet("background-color: #4c4c4c")
		self.pressure.setChecked(True)
		


		self.column_layout1a = QtGui.QVBoxLayout(self.tab1) 
		self.column_layout1a1 = QtGui.QHBoxLayout()
		self.column_layout1a1.addWidget(self.sim_type_label)
		self.column_layout1a1.addWidget(self.chamber_button)
		self.column_layout1a1.addWidget(self.parcel_button)
		self.column_layout1a.addLayout(self.column_layout1a1)
		self.column_layout1a.addLayout(self.simulation_type_layout)
		self.column_layout1a.addLayout(self.int_con_layout)
		self.column_layout1a.addWidget(self.runtime_label)
		self.column_layout1a.addWidget(self.runtime)
			
		self.aer_props_layout = QtGui.QVBoxLayout(self.tab2)
		self.aer_props_layout.addLayout(self.num_modes_layout)
		self.aer_props_layout.addLayout(self.grid)
		self.aer_props_layout.addLayout(self.grid_comp)
		self.aer_props_layout.addWidget(self.aerosol_composition_button)
		self.aer_props_layout.addWidget(self.add_aerosol_type_button)
		self.aer_props_layout.addWidget(self.freezing_criteria)

		self.advanced_settings_layout = QtGui.QVBoxLayout(self.tab3)
		self.Dlow_layout = QtGui.QHBoxLayout()
		self.Dlow_layout.addWidget(self.Dlow_label)
		self.Dlow_layout.addWidget(self.Dlow)
		self.advanced_settings_layout.addLayout(self.Dlow_layout)
		self.rkm_layout = QtGui.QHBoxLayout()
		self.rkm_layout.addWidget(self.rkm_label)
		self.rkm_layout.addWidget(self.rkm)
		self.advanced_settings_layout.addLayout(self.rkm_layout)
		self.time_step_layout = QtGui.QHBoxLayout()
		self.time_step_layout.addWidget(self.time_step_label)
		self.time_step_layout.addWidget(self.time_step)
		self.alpha_crit_layout = QtGui.QHBoxLayout()
		self.alpha_crit_layout.addWidget(self.alpha_crit_label)
		self.alpha_crit_layout.addWidget(self.alpha_crit)
		self.advanced_settings_layout.addLayout(self.time_step_layout)
		self.advanced_settings_layout.addLayout(self.alpha_crit_layout)


		self.column_layout1.addWidget(self.tab_widget)
		self.column_layout1.addStretch(1)


		self.column_layout1b = QtGui.QHBoxLayout()
		self.column_layout1b.addWidget(self.savebutton)
		self.column_layout1.addLayout(self.column_layout1b)
		self.column_layout1c = QtGui.QHBoxLayout()
		self.column_layout1c.addWidget(self.submit_button)
		self.column_layout1.addLayout(self.column_layout1c)

		self.column_layout2.addLayout(self.layout_fig)
		self.column_layout2.addStretch(1)
		self.column_layout3a = QtGui.QHBoxLayout()
		self.column_layout3a.addWidget(self.RH)
		self.column_layout3a.addWidget(self.temperature)
		self.column_layout3a.addWidget(self.pressure)
		self.column_layout3.addLayout(self.column_layout3a)
		self.column_layout3b = QtGui.QHBoxLayout()
		self.column_layout3b.addWidget(self.ice_total)
		self.column_layout3b.addWidget(self.drop_total)
		self.column_layout3.addLayout(self.column_layout3b)
		self.column_layout3.addWidget(self.plot_button)


		self.layout.addLayout(self.column_layout1)
		self.layout.addLayout(self.column_layout2)
		self.layout.addLayout(self.column_layout3)
		self.plot()
        # Set the VBox layout as the window's main layout
		#self.layout.addStretch(1)
		self.setLayout(self.layout)

		

	def number_modes_enter_clicked(self):
		
		f.clearlayout(self,self.grid)

		rows = ['Number','Diameter','Sigma']
		columns = ['Mode_'+str(x+1) for x in range(int(self.number_modes.text()))]
		defaults = ['100e6','100e-9','0.5']

		for i, r in enumerate(rows):
			self.label = QtGui.QLabel()
			self.label.setText(r)
			self.grid.addWidget(self.label,i+1, 0)
			for j, c in enumerate(columns):
				self.label = QtGui.QLabel()
				self.label.setText(c)
				dummy = r+'_'+c
				setattr(self,dummy,QtGui.QLineEdit(defaults[i]))
				self.grid.addWidget(self.label,0, j+1)
				self.grid.addWidget(getattr(self,dummy),i+1, j+1)
			
		return self

	def aerosol_comp_clicked(self):

		from namelist import Mass_frac

		f.clearlayout(self,self.grid_comp)

		rows = [key for key in Mass_frac.keys()]
		columns = range(int(self.number_modes.text()))
		

		for i, r in enumerate(rows):
			self.label = QtGui.QLabel()
			self.label.setText(r)
			self.grid.addWidget(self.label,i+4, 0)
			for j in columns: # set default to AS for each mode
				dummy = r+'_'+str(j)
				if i == 0:
					setattr(self,dummy,QtGui.QLineEdit('1.0'))
				else:
					setattr(self,dummy,QtGui.QLineEdit('0.0'))
				
				self.grid.addWidget(getattr(self,dummy),i+4, j+1)
			
		return self
		
	def load_data_clicked(self):
		pass

	def advanced_settings_clicked(self):
		pass

	def submit_clicked(self):
		""" this function writes namelist variables inputted in gui
		    to a pickle file which is read in by pyACPIM in variables.py
		"""
		import namelist as n

		try:
			columns = ['Mode_'+str(x+1) for x in range(int(self.number_modes.text()))]
		except:
			print('please enter number of modes')
			return

		number = ['Number_'+ c for c in columns]
		diameter = ['Diameter_'+c for c in columns]
		sigma = ['Sigma_'+c for c in columns]

		number_mode = [float(getattr(self,name).text()) for name in number]
		diameter_mode = [float(getattr(self,name).text()) for name in diameter]
		sigma_mode = [float(getattr(self,name).text()) for name in sigma]

		for key in n.Mass_frac.keys():
			n.Mass_frac[key] = [float((getattr(self,key+'_'+str(i)).text())) for i in range(int(self.number_modes.text()))]
		
		# edit namelist values with values from gui
		nbins = 70
		nmodes = int(self.number_modes.text())
		NAER = number_mode
		D_AER = diameter_mode
		sigma = sigma_mode
		RH = float(self.RH_int.text())
		T = float(self.T.text())
		P = float(self.P.text())*100
		
		if self.chamber_button.isChecked():
			simulation_type = 'Chamber'
		else:
			simulation_type = 'Parcel'
		
		try:
			w = float(self.updraft.text())
		except:
			w = 1

		runtime = int(round(float((self.runtime.text()))))
		Dlow = float(self.Dlow.text())
		rkm = float(self.rkm.text())
		dt = max(int(round(float(self.time_step.text()))),1)
		try:
			PRESS1 = float(self.press1.text())
			PRESS2 = float(self.press2.text())
			Temp1 = float(self.temp1.text())
			Temp2 = float(self.temp2.text())
		except:
			PRESS1 = 1
			PRESS2 = 1
			Temp1 = 1
			Temp2 = 1

		alpha = float(self.alpha_crit.text())

		if self.freezing_criteria.currentText() == 'heterogeneous freezing criteria':
			freezing_crit = 'RH>1'
		else:
			freezing_crit = self.freezing_criteria.currentText()
		
		# write namelist values to pickle file
		list_vars = [nbins, nmodes, NAER, D_AER, sigma, n.Mass_frac, RH, T, P,
		            self.filename,simulation_type, runtime, Dlow, rkm, dt, 
		            PRESS1, PRESS2, Temp1, Temp2, freezing_crit, alpha, w]
		
		with open('test.pkl','wb') as f:
			pickle.dump(list_vars,f)
		print(list_vars)
		# run model
		model.run()

	
	def add_aerosol_clicked(self):
		pass

	def plot(self):
		params = {'ytick.color' : 'gray',
		'xtick.color': 'gray',
		'axes.labelcolor': 'gray',
		'axes.edgecolor':'gray'}
		plt.rcParams.update(params)

		


		# see how many checkboxes have been selected
		vars= ['RH','pressure','temperature','ice_total','drop_total']
		to_plot = [] # things to be plotted
		for var in vars:
			if (getattr(self,var).isChecked()):
				to_plot.append(var)
			else:
				pass


		# get output to plot
		try:
			nc = Dataset(self.filename)
		except:
			nc = Dataset('new_output.nc')
			to_plot = ['RH','pressure','temperature']

		# clear perviously plotted stuff
		self.figure.clear()
		# plot sublots
		for i, thing in enumerate(to_plot):
			n = 100*(len(to_plot))+10+(i+1)
			ax = self.figure.add_subplot(n, axisbg='#4c4c4c')
			ax.plot(nc[thing][:],'y')
			ax.set_title(thing,color='gray')

		self.figure.set_facecolor('#4c4c4c')	
		# sort out subplot layout
		self.figure.tight_layout()   
		# make plot visible on canvas
		self.canvas.draw()
		# close nc file
		nc.close()

	def onParcelClicked(self):
		
		f.clearlayout(self,self.parcel_layout)
		
		self.parcel_layout = QtGui.QHBoxLayout()
			
		f.clearlayout(self,self.chamber_layout_temp)
		f.clearlayout(self,self.chamber_layout_press)

		self.updraft = QtGui.QLineEdit('1')
		self.updraft_label = QtGui.QLabel()
		self.updraft_label.setText('Upraft (ms-1)')
		self.parcel_layout.addWidget(self.updraft_label)
		self.parcel_layout.addWidget(self.updraft)

		self.simulation_type_layout.addLayout(self.parcel_layout)
		
		return self.parcel_layout

	def onChamberClicked(self):

		f.clearlayout(self, self.chamber_layout_temp)
		f.clearlayout(self,self.chamber_layout_press)

		self.chamber_layout_temp = QtGui.QHBoxLayout()
		self.chamber_layout_press = QtGui.QHBoxLayout()
		self.chamber_layout = QtGui.QVBoxLayout()

		f.clearlayout(self,self.parcel_layout)

		params = ['temp1','temp2','press1','press2']
		defaults = ['6.9','0.03','628.12','2.6e-3']

		for i, param in enumerate(params):
			setattr(self,param,QtGui.QLineEdit(defaults[i]))
			self.param_label = QtGui.QLabel()
			self.param_label.setText(param)
			if param[0] == 't':
				self.chamber_layout_temp.addWidget(self.param_label)
				self.chamber_layout_temp.addWidget(getattr(self,param))
			else:
				self.chamber_layout_press.addWidget(self.param_label)
				self.chamber_layout_press.addWidget(getattr(self,param))
		
		self.chamber_layout.addLayout(self.chamber_layout_temp)
		self.chamber_layout.addLayout(self.chamber_layout_press)
		self.simulation_type_layout.addLayout(self.chamber_layout)
						
		return self.chamber_layout

	def onSaveButtonClicked(self):
		self.filename, filter = QtGui.QFileDialog.getSaveFileName(parent=self, caption='Select output file', dir='.', filter='netcdf files (*.nc)')
	
	def run(self):
        # Show the gui
		self.show()
        # Run the qt application
		qt_app.exec_()


 
# Create an instance of the application window and run it
app = main_window()
app.run()