import sys
from PySide2 import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import pandas as pd
import matplotlib.pyplot as plt
import gui_functions as f
import pickle
import numpy as np
from threading import Thread
from netCDF4 import Dataset
import main as model
import constants as c
import time
import threading
import importlib


qt_app = QtWidgets.QApplication(sys.argv)
finished = QtCore.Signal()

class main_window(QtWidgets.QMainWindow):
    
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
        self.resize(5000, 1500) # set size of gui window
        
        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.Background,QtGui.QColor('#4c4c4c'))
        self.setPalette(palette)
        self.setWindowTitle('pyACPIM')
        # set major layouts
        self.central_widget_layout = QtWidgets.QHBoxLayout()
        self.column_layout1 = QtWidgets.QVBoxLayout()
        self.column_layout2 = QtWidgets.QVBoxLayout()
        self.column_layout3 = QtWidgets.QVBoxLayout()
        self.simulation_type_layout = QtWidgets.QVBoxLayout()
        self.parcel_layout = QtWidgets.QHBoxLayout()
        self.chamber_layout = QtWidgets.QHBoxLayout()
        self.description_layout = QtWidgets.QVBoxLayout()
        self.chamber_layout_temp = QtWidgets.QHBoxLayout()
        self.chamber_layout_press = QtWidgets.QHBoxLayout()
        self.gbox_layout = QtWidgets.QVBoxLayout()
        self.grid = QtWidgets.QGridLayout()
        self.grid_comp = QtWidgets.QGridLayout()

        # set-up tabs
        self.tab_widget = QtWidgets.QTabWidget()
        self.tab1 = QtWidgets.QWidget()
    #   self.tab1.resize(2000,1500)
        self.tab_widget.setPalette(palette)
        self.tab_widget.addTab(self.tab1,'Initial Conditions')
        self.tab2 = QtWidgets.QWidget()
        self.tab_widget.addTab(self.tab2,'Aerosol Properties')
        self.tab3 = QtWidgets.QWidget()
        self.tab_widget.addTab(self.tab3,'Advanced Settings')
        # set the style sheet for the GUI
        self.setStyleSheet(self.getStyleSheet('./stylesheet.qss'))

        self.submit_clicked_has_happend = False

        # Create widgets
        # column 1
        # TAB 1
        self.chamber_button = QtWidgets.QRadioButton('&Chamber', self)
        self.chamber_button.toggled.connect(self.onChamberClicked)
        self.chamber_button.setChecked(True)
        self.parcel_button = QtWidgets.QRadioButton('&Parcel', self)
        self.parcel_button.toggled.connect(self.onParcelClicked)
        self.sim_type_label = QtWidgets.QLabel('Simulation type')
        self.RH_int = QtWidgets.QLineEdit('0.9')
        self.RH_label = QtWidgets.QLabel()
        self.P = QtWidgets.QLineEdit('950')
        self.P_label = QtWidgets.QLabel('Pressure (hPa)')
        self.RH_label.setText('RH')
        self.T = QtWidgets.QLineEdit('273')
        self.T_label = QtWidgets.QLabel('Temperature (K)')
        self.runtime = QtWidgets.QLineEdit('270')
        self.runtime_label = QtWidgets.QLabel('Runtime (seconds)')
        # TAB 2
        self.number_modes = QtWidgets.QLineEdit('1')
        self.number_modes_enter = QtWidgets.QPushButton('&Enter',self)
        self.number_modes_enter.setChecked(True)
        self.number_modes_enter.clicked.connect(self.number_modes_enter_clicked)
        self.number_modes_enter.clicked.connect(self.aerosol_comp_clicked)
        self.number_modes_enter.clicked.connect(self.make_add_another_aer_type_button)
        self.number_modes_label = QtWidgets.QLabel()
        self.number_modes_label.setText('number of modes')
        self.num_modes_layout = QtWidgets.QHBoxLayout()
        self.num_modes_layout.addWidget(self.number_modes_label)
        self.num_modes_layout.addWidget(self.number_modes)
        self.num_modes_layout.addWidget(self.number_modes_enter)
        self.freezing_criteria = QtWidgets.QComboBox(self)
        self.freezing_criteria_label = QtWidgets.QLabel()
        self.freezing_criteria_label.setText('Heterogeneous freezing criteria')
        self.freezing_criteria.addItems(['RH>1  [default]','Activated Drops', 'Threshold water mass'])
        self.freezing_criteria.currentIndexChanged.connect(self.make_alpha_crit_lineedit)
        # TAB 3
        self.Dlow = QtWidgets.QLineEdit('10e-9')
        self.Dlow_label = QtWidgets.QLabel('Size of smallest aerosol size bin (m)')
        self.time_step = QtWidgets.QLineEdit('1')
        self.time_step_label = QtWidgets.QLabel('model time step (seconds)')
        #
        self.submit_button = QtWidgets.QPushButton('&Run',self)
        self.submit_button.clicked.connect(self.submit_clicked)
        self.savebutton = QtWidgets.QPushButton('&Select Output file', self)
        self.savebutton.clicked.connect(self.onSaveButtonClicked)
        self.progressBar = QtWidgets.QProgressBar(self)
        self.progressBar.setGeometry(200, 80, 250, 20)
        self.progressBar.setPalette(palette)
        
        # column 2
        # canvas for plotting model results
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        # set the layout of figure/canvas
        self.layout_fig = QtWidgets.QVBoxLayout()
        self.layout_fig.addWidget(self.toolbar)
        self.layout_fig.addWidget(self.canvas)
    
        # column 3
        self.Select_file_to_load_label = QtWidgets.QLabel()
        self.Select_file_to_load_label.setText('Select data file')
        self.Select_file_to_load = QtWidgets.QComboBox(self)
        self.Select_file_to_load.addItems(['Select input file type','Temperature','Pressure', 'CDP','FSSP'])
        self.Select_file_to_load.currentIndexChanged.connect(self.load_file)
    
        self.plot_button = QtWidgets.QPushButton('&Plot')
        self.plot_button.clicked.connect(self.plot)
        self.plot_size_dis_time = QtWidgets.QLineEdit('0')
        self.plot_size_dis_time_label = QtWidgets.QLabel('Time to plot size distribution (s)')
        self.plot_size_dis_button = QtWidgets.QPushButton('&Plot Size Distribution')
        self.plot_size_dis_button.clicked.connect(self.plot_size_dis)
        self.simulation_duration = QtWidgets.QLineEdit('270')
        self.simulation_duration_label = QtWidgets.QLabel('Expansion Duration (seconds)')
        self.simulation_start_time = QtWidgets.QLineEdit('12:59:33') # this needs to be changed to hh:mm:ss
        self.simulation_start_date = QtWidgets.QLineEdit('2014-04-03') # this needs to be changed to yyyy-mm-dd
        self.simulation_start_date_label = QtWidgets.QLabel()
        self.simulation_start_date_label.setText('Expansion start date')
        self.simulation_start_time_label = QtWidgets.QLabel(' and time')
        self.RH          = QtWidgets.QCheckBox('RH',self)
        self.temperature = QtWidgets.QCheckBox('Temperature',self)
        self.drop_total  = QtWidgets.QCheckBox('Number of drops',self)
        self.ice_total   = QtWidgets.QCheckBox('Number of ice crystals',self)
        self.pressure    = QtWidgets.QCheckBox('Pressure',self)
        self.liquid_water_content = QtWidgets.QCheckBox('LWC',self)
        self.RH.setChecked(True)
        self.temperature.setChecked(True)
        self.pressure.setChecked(True)
        self.temp_obs_checkbox = QtWidgets.QCheckBox('Temperature', self)
        self.press_obs_checkbox = QtWidgets.QCheckBox('Pressure',self)
        self.ice_obs_checkbox = QtWidgets.QCheckBox('Ice number', self)
        self.total_obs_checkbox = QtWidgets.QCheckBox('Total number',self)
        self.with_model_radiobutton = QtWidgets.QRadioButton('plot with model output',self)
        self.plot_obs_button = QtWidgets.QPushButton('&Plot')
        self.plot_obs_button.clicked.connect(self.plot_obs)

        # LAYOUTS
        # column 1
        # TAB 1
        self.int_con_layout = QtWidgets.QVBoxLayout()
        self.RH_layout = QtWidgets.QHBoxLayout()
        self.RH_layout.addWidget(self.RH_label)
        self.RH_layout.addWidget(self.RH_int)
        self.P_layout = QtWidgets.QHBoxLayout()
        self.P_layout.addWidget(self.P_label)
        self.P_layout.addWidget(self.P)
        self.T_layout = QtWidgets.QHBoxLayout()
        self.T_layout.addWidget(self.T_label)
        self.T_layout.addWidget(self.T)
        self.int_con_layout.addLayout(self.RH_layout)
        self.int_con_layout.addLayout(self.P_layout)
        self.int_con_layout.addLayout(self.T_layout)
        self.int_con_layout.addStretch()
        self.column_layout1a = QtWidgets.QVBoxLayout(self.tab1) 
        self.column_layout1a.setSpacing(15)
        self.column_layout1a1 = QtWidgets.QHBoxLayout()
        self.column_layout1a1.addWidget(self.chamber_button)
        self.column_layout1a1.addWidget(self.parcel_button)
        self.gbox = QtWidgets.QGroupBox('Simulation type')
        self.gbox.setLayout(self.column_layout1a1)
        self.column_layout1a.addWidget(self.gbox)
        self.column_layout1a.addStretch()
        self.column_layout1a.addLayout(self.simulation_type_layout)
        self.column_layout1a.addStretch()
        self.column_layout1a.addLayout(self.int_con_layout)
        self.column_layout1a2 = QtWidgets.QHBoxLayout()
        self.column_layout1a2.addWidget(self.runtime_label)
        self.column_layout1a2.addWidget(self.runtime)
        self.column_layout1a.addLayout(self.column_layout1a2)
        # TAB 2 
        self.aer_props_layout = QtWidgets.QVBoxLayout(self.tab2)
        self.aer_props_layout.addLayout(self.num_modes_layout)
        self.aer_props_layout.addStretch()
        self.aer_props_layout.addLayout(self.grid)
        self.aer_props_layout.addStretch()
        self.aer_props_layout.addLayout(self.grid_comp)
        self.aer_comp_layout = QtWidgets.QVBoxLayout()
        self.aer_props_layout.addLayout(self.aer_comp_layout)
        self.aer_props_layout.addStretch()
        # TAB 3
        self.advanced_settings_layout = QtWidgets.QVBoxLayout(self.tab3)
        self.Dlow_layout       = QtWidgets.QHBoxLayout()
        self.time_step_layout  = QtWidgets.QHBoxLayout()
        self.Dlow_layout.addWidget(self.Dlow_label)
        self.Dlow_layout.addWidget(self.Dlow)
        self.time_step_layout.addWidget(self.time_step_label)
        self.time_step_layout.addWidget(self.time_step)
        self.advanced_settings_layout.addLayout(self.Dlow_layout)
        self.advanced_settings_layout.addLayout(self.time_step_layout)
        self.het_freezing_layout = QtWidgets.QHBoxLayout()
        self.het_freezing_layout.addWidget(self.freezing_criteria_label)
        self.het_freezing_layout.addWidget(self.freezing_criteria)
        self.advanced_settings_layout.addLayout(self.het_freezing_layout)
        self.advanced_settings_layout.addStretch()

        self.column_layout1.addWidget(self.tab_widget)
        self.column_layout1b = QtWidgets.QHBoxLayout()
        self.column_layout1b.addWidget(self.savebutton)
        self.column_layout1.addLayout(self.column_layout1b)
        self.column_layout1c = QtWidgets.QHBoxLayout()
        self.column_layout1c.addWidget(self.submit_button)
        self.column_layout1c.addWidget(self.progressBar)
        self.column_layout1.addLayout(self.column_layout1c)

        # column 2
        self.column_layout2.addLayout(self.layout_fig)

        # column 3
        self.column_layout3gbox = QtWidgets.QVBoxLayout()
        self.column_layout3a = QtWidgets.QHBoxLayout()
        self.column_layout3a.addWidget(self.RH)
        self.column_layout3a.addWidget(self.temperature)
        self.column_layout3a.addWidget(self.pressure)
        self.column_layout3gbox.addLayout(self.column_layout3a)
        self.column_layout3b = QtWidgets.QHBoxLayout()
        self.column_layout3b.addWidget(self.ice_total)
        self.column_layout3b.addWidget(self.drop_total)
        self.column_layout3b.addWidget(self.liquid_water_content)
        self.column_layout3gbox.addLayout(self.column_layout3b)
        self.column_layout3gbox.addWidget(self.plot_button)
        self.column_layout3e = QtWidgets.QHBoxLayout()
        self.column_layout3e.addWidget(self.plot_size_dis_time_label)
        self.column_layout3e.addWidget(self.plot_size_dis_time)
        self.column_layout3gbox.addLayout(self.column_layout3e)
        self.column_layout3gbox.addWidget(self.plot_size_dis_button)
        self.gbox_model_plotting = QtWidgets.QGroupBox('Plotting model data')
        self.gbox_model_plotting.setLayout(self.column_layout3gbox)
        self.column_layout3.addWidget(self.gbox_model_plotting)
        self.column_layout3.addStretch()
        self.column_layout3e = QtWidgets.QHBoxLayout()
        self.column_layout3e.addWidget(self.Select_file_to_load_label)
        self.column_layout3e.addWidget(self.Select_file_to_load)
        self.column_layout3.addLayout(self.column_layout3e)
        self.column_layout3c = QtWidgets.QHBoxLayout()
        self.column_layout3c.addWidget(self.simulation_start_date_label)
        self.column_layout3c.addWidget(self.simulation_start_date)
        self.column_layout3c.addWidget(self.simulation_start_time_label)
        self.column_layout3c.addWidget(self.simulation_start_time)
        self.column_layout3.addLayout(self.column_layout3c)
        self.column_layout3d = QtWidgets.QHBoxLayout()
        self.column_layout3d.addWidget(self.simulation_duration_label)
        self.column_layout3d.addWidget(self.simulation_duration)
        self.column_layout3.addLayout(self.column_layout3d)
        self.column_layout3.addStretch()
        self.obs_plotting_layout = QtWidgets.QVBoxLayout()
        self.obs_plotting_layout_1 = QtWidgets.QHBoxLayout()
        self.obs_plotting_layout_2 = QtWidgets.QHBoxLayout()
        self.obs_plotting_layout_3 = QtWidgets.QHBoxLayout()
        self.obs_plotting_layout_1.addWidget(self.temp_obs_checkbox)
        self.obs_plotting_layout_1.addWidget(self.press_obs_checkbox)
        self.obs_plotting_layout.addLayout(self.obs_plotting_layout_1)
        self.obs_plotting_layout_2.addWidget(self.ice_obs_checkbox)
        self.obs_plotting_layout_2.addWidget(self.total_obs_checkbox)
        self.obs_plotting_layout.addLayout(self.obs_plotting_layout_2)
        self.obs_plotting_layout_3.addWidget(self.with_model_radiobutton)
        self.obs_plotting_layout_3.addWidget(self.plot_obs_button)
        self.obs_plotting_layout.addLayout(self.obs_plotting_layout_3)
        self.gbox_obs_plotting = QtWidgets.QGroupBox('Plotting observations')
        self.gbox_obs_plotting.setLayout(self.obs_plotting_layout)
        self.column_layout3.addWidget(self.gbox_obs_plotting)
        self.column_layout3.addStretch()
        # add all layouts onto central widget
        self.central_widget_layout.addLayout(self.column_layout1,1)
        self.central_widget_layout.addLayout(self.column_layout2,2)
        self.central_widget_layout.addLayout(self.column_layout3,1)
        self.plot()
        # Set the VBox layout as the window's main layout
        self.central_widget.setLayout(self.central_widget_layout)

    def make_add_another_aer_type_button(self):
        f.clearlayout(self, self.aer_comp_layout)
        self.add_aerosol_type_button = QtWidgets.QPushButton('&Add another aerosol type',self)
        self.add_aerosol_type_button.clicked.connect(self.open_new_dialog)
        self.aer_comp_layout.addWidget(self.add_aerosol_type_button)
        return self

    def make_alpha_crit_lineedit(self):
        if self.freezing_criteria.currentText() == 'Threshold water mass':
            self.alpha_crit, ok = QtWidgets.QInputDialog.getText(None, 'Threshold mass of water', 'Enter a value for A such that Threshold mass of water = A*Volume_INP ',text='70')
        return self

    def load_file(self):
        if self.Select_file_to_load.currentText() == 'Temperature':
            self.onLoadTempButtonClicked()
        if self.Select_file_to_load.currentText() == 'Pressure':
            self.onLoadPressButtonClicked()
        if self.Select_file_to_load.currentText() == 'CDP':
            self.onLoadCDPClicked()
        if self.Select_file_to_load.currentText() == 'FSSP':
            self.onLoadFSSPClicked()

    
    def number_modes_enter_clicked(self):
        
        f.clearlayout(self,self.grid)

        rows = ['Number','Diameter','Sigma']
        columns = ['Mode_'+str(x+1) for x in range(int(self.number_modes.text()))]
        defaults = ['100e6','100e-9','0.5']

        for i, r in enumerate(rows):
            self.label = QtWidgets.QLabel()
            self.label.setText(r)
            self.grid.addWidget(self.label,i+1, 0)
            for j, c in enumerate(columns):
                self.label = QtWidgets.QLabel()
                self.label.setText(c)
                dummy = r+'_'+c
                setattr(self,dummy,QtWidgets.QLineEdit(defaults[i]))
                self.grid.addWidget(self.label,0, j+1)
                self.grid.addWidget(getattr(self,dummy),i+1, j+1)
        
        
        return self

    def aerosol_comp_clicked(self):
        importlib.reload(c)

        f.clearlayout(self,self.grid) # can't delete this here as model will not run
        f.clearlayout(self,self.grid_comp)

        self.number_modes_enter_clicked()

        rows = [key for key in c.aerosol_dict.keys()]
        rows = ['Mass Fraction'] + rows # 
        columns = range(int(self.number_modes.text()))
        
        for i, r in enumerate(rows):
            self.label = QtWidgets.QLabel()
            self.label.setText(r)
            if i == 0:
                self.label.setAlignment(QtCore.Qt.AlignCenter)
                myFont=QtGui.QFont()
                myFont.setBold(True)
                self.label.setFont(myFont)
            self.grid.addWidget(self.label,i+4, 0)
            for j in columns: # set default to AS for each mode
                dummy = r+'_'+str(j)
                if i == 0:
                    continue
                if i == 1:
                    setattr(self,dummy,QtWidgets.QLineEdit('1.0'))
                else:
                    setattr(self,dummy,QtWidgets.QLineEdit('0.0'))
                
                self.grid.addWidget(getattr(self,dummy),i+4, j+1)
            
        return self
        
    def onLoadTempButtonClicked(self):
        self.temp_filename, filter = QtWidgets.QFileDialog.getOpenFileName(parent=self, caption='Select temperature file', dir='.', filter='text files (*.txt *.csv)')
    def onLoadPressButtonClicked(self):
        self.press_filename, filter = QtWidgets.QFileDialog.getOpenFileName(parent=self, caption='Select pressure file', dir='.', filter='text files (*.txt *.S30)')
    def onLoadCDPClicked(self):
        self.CDP_filename, filter = QtWidgets.QFileDialog.getOpenFileName(parent=self, caption='Select CDP file', dir='.', filter='text files (*.min)')
        self.CDP_data_total, self.CDP_data_ice = f.getCDPVar(self,self.CDP_filename)
        
    def onLoadFSSPClicked(self):
        self.FSSP_filename, filter = QtWidgets.QFileDialog.getOpenFileName(parent=self, caption='Select FSSP file', dir='.', filter='text files (*.min)')


    def advanced_settings_clicked(self):
        pass

    def submit_clicked(self):
        """ this function writes namelist variables inputted in gui
            to a pickle file which is read in by pyACPIM in variables.py
        """
        if self.submit_clicked_has_happend:
            self.box = QtWidgets.QMessageBox()
            self.box.setText('Warning '+self.filename+' will be overwritten')
            self.box.setInformativeText('Do you want to select a new output file?')
            self.box.setStandardButtons(self.box.Yes | self.box.No)
            self.box.setDefaultButton(self.box.Yes)
        #   self.box.buttonClicked.connect(self.function_test)

            ret = self.box.exec_()

            if ret == QtWidgets.QMessageBox.Yes:
                self.onSaveButtonClicked()
                self.submit_clicked_has_happend = False
                return

        self.submit_clicked_has_happend = True

        # set a spinning wheel cursor while model is running to stop user interacting with the GUI while
        # it is busy running the model
        kursor = QtGui.QCursor()
        kursor.setShape(QtCore.Qt.WaitCursor)
        QtWidgets.QApplication.setOverrideCursor(kursor)
        ###################################################

        import namelist as n

        try:
            print(self.Number_Mode_1)
        except: 
            self.box = QtWidgets.QMessageBox()
            self.box.setText('please enter number of modes')
            self.box.exec_()
            
        columns = ['Mode_'+str(x+1) for x in range(int(self.number_modes.text()))]          

        number = ['Number_'+ c for c in columns]
        diameter = ['Diameter_'+c for c in columns]
        sigma = ['Sigma_'+c for c in columns]

        number_mode = [float(getattr(self,name).text()) for name in number]
        diameter_mode = [float(getattr(self,name).text()) for name in diameter]
        sigma_mode = [float(getattr(self,name).text()) for name in sigma]
        nmodes = int(self.number_modes.text())

        Mass_frac = {}
        aer_types = c.aerosol_dict.keys()
        try:
            for aer_type in aer_types:
                Mass_frac[aer_type] = np.zeros(nmodes)

            for key in Mass_frac.keys():
                Mass_frac[key] = [float((getattr(self,key+'_'+str(i)).text())) for i in range(int(self.number_modes.text()))]
        except:
            self.box = QtWidgets.QMessageBox()
            self.box.setText('Please enter aerosol type')
            self.box.exec_()

        self.mass_frac_check_sum(nmodes,Mass_frac)
            
        # edit namelist values with values from gui
        nbins = 70
        #nmodes = int(self.number_modes.text())
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
        dt = max(int(round(float(self.time_step.text()))),1)
        
        if self.chamber_button.isChecked():
            try:
                PRESS1 = float(self.press_a_lineedit.text())
                PRESS2 = float(self.press_b_lineedit.text())
                Temp1 = float(self.temp_a_lineedit.text())
                Temp2 = float(self.temp_b_lineedit.text())
            except:
                self.box = QtWidgets.QMessageBox()
                self.box.setText('Error reading pressure and temperature a and b values')
                self.box.setInformativeText('Please check a and b values before trying to run again')
                self.box.exec_()
        else:
            PRESS1 = 1
            PRESS2 = 1
            Temp1 = 1
            Temp2 = 1
        try:
            alpha = float(self.alpha_crit.text())
        except:
            alpha = float('70')

        freezing_crit = self.freezing_criteria.currentText()
        if self.freezing_criteria.currentText() == 'RH>1    [default]':
            freezing_crit = 'RH>1'
        try:
            print(self.filename)
        except:
            self.box = QtWidgets.QMessageBox()
            self.box.setText('Output file name not selected')
            self.box.setInformativeText('Please select output file')
            self.box.exec_()

        # write namelist values to pickle file
        list_vars = [nbins, nmodes, NAER, D_AER, sigma, Mass_frac, RH, T, P,
                    self.filename,simulation_type, runtime, Dlow, dt, 
                    PRESS1, PRESS2, Temp1, Temp2, freezing_crit, alpha, w]
        
        with open('test.pkl','wb') as f:
            pickle.dump(list_vars,f)
        
        # run model
        self.progressBar.setRange(0,runtime)
        self.progressBar.setValue(0)

        t = threading.Thread(target=model.run(self.progressBar))
        t.daemon = True
        t.start()
        # return to normal cursor once model has finished running
        QtWidgets.QApplication.restoreOverrideCursor()
        self.progressBar.setValue(runtime)
    
    def mass_frac_check_sum(self,nmodes,Mass_frac):
        total = np.zeros(nmodes)

        for mode in range(nmodes):
            for key in Mass_frac.keys():
                total[mode] += Mass_frac[key][mode]
            if total[mode] != 1:
                self.box = QtWidgets.QMessageBox()
                self.box.setText('sum of mass fraction does not equal zero')
                self.box.setInformativeText('Please correct mass fraction values')
                self.box.exec_()
                self.submit_clicked_has_happend = False
                return

        
    def getTData(self):
        try:
            f.getTempVar(self, self.temp_filename)
        except:
            self.box = QtWidgets.QMessageBox()
            self.box.setText('Please select input temperature file')
            self.box.exec_()
    def getPData(self):
        try:
            f.getPressVar(self, self.press_filename)
        except:
            self.box = QtWidgets.QMessageBox()
            self.box.setText('Please select input pressure file')
            self.box.exec_()

    def errormessage(self,error_string):
        self.box = QtWidgets.QMessageBox()
        self.box.setText(error_string)
        self.box.exec_()
    
    def plot_obs(self):
        self.figure.clear()
        if self.temp_obs_checkbox.isChecked():
            try:
                f.getTempVar(self, self.temp_filename)
            except:
                self.errormessage('Please select temperature file')
                return
            
        if self.press_obs_checkbox.isChecked():
            try:
                f.getPressVar(self,self.press_filename)
            except:
                self.errormessage('Please select pressure file')
                return
        
        params = {'ytick.color' : '#d1cfcf',
        'xtick.color': '#d1cfcf',
        'axes.labelcolor': '#d1cfcf',
        'axes.edgecolor':'gray'}
        plt.rcParams.update(params)
        
        start = str(self.simulation_start_date.text())+' '+str(self.simulation_start_time.text())
        end = pd.Timestamp(start) + pd.Timedelta(seconds = float(self.simulation_duration.text()))
        
        try:
            x = np.array(range(len(self.Pressure_data[start:end])))
        except:
            try:
                x = np.array(range(len(self.Temperature_data[start:end])))
            except:
                self.errormessage('nothing to plot')


        # see how many checkboxes have been selected
        vars2= ['temp_obs_checkbox','press_obs_checkbox','ice_obs_checkbox','total_obs_checkbox'] # need to add LWC and IWC
        vars= ['temperature','pressure','ice_total','liq_total']#,'liquid_water_content']
        titles = ['Temperature','Pressure','ice crystal number concentration',
                    'total number concentration']
        units = ['K','hPa','g kg$^{-1}$','g kg$^{-1}$']
        unit_conversions = [1,1e-2,1e-6,1e-6]
        to_plot = [] # things to be plotted
        plot_title = []
        ylabels = []
        conversions = []
        for var2, var, title, unit, unit_conversion in zip(vars2,vars,titles, units, unit_conversions):
            if (getattr(self,var2).isChecked()):
                to_plot.append(var)
                plot_title.append(title)
                ylabels.append(unit)
                conversions.append(unit_conversion)
            else:
                pass
        # get output to plot
        try:
            nc = Dataset(self.filename)
        except:
            nc = Dataset('DEFAULT_OUTPUT.nc')   
        # clear previously plotted stuff
        self.figure.clear()
        try:
            xT = len(self.Temperature_data['Temperature'][start:end])
        except:
            pass
        try:    
            xP = len(self.Pressure_data['Pressure'][start:end])
        except:
            pass
        if self.ice_obs_checkbox.isChecked():
            xCDP = len(self.CDP_data_ice[start:end])
        if self.total_obs_checkbox.isChecked():
            try:
                xCDP = len(self.CDP_data_ice[start:end])
            except:
                self.errormessage('Please select CDP file')

        data1 = np.zeros(int(self.runtime.text()))
        data2 = np.zeros(int(self.runtime.text()))
        # plot sublots
        for i, thing in enumerate(to_plot):
            n = 100*(len(to_plot))+10+(i+1) # set-up number of subplot
            ax = self.figure.add_subplot(n, facecolor='#4c4c4c')
            if self.with_model_radiobutton.isChecked() : 
                if thing == 'liq_total':
                    for n in range(2):
                        data1 = data1 + np.sum(np.where(nc['ice_mass'][:,n,:]>3.61e-14,nc['ice_number'][:,n,:],0),axis=1)
                        data2 = data2 + np.sum(np.where(nc['liq_mass'][:,n,:]>3.61e-14,nc['liq_number'][:,n,:],0),axis=1)
                    ax.plot((data1+data2)*conversions[i],'y',label='model')
                else:
                    ax.plot(nc[thing][:]*conversions[i],'y', label='model')
            if thing == 'temperature':
                ax.plot(np.array(range(xT)),self.Temperature_data['Temperature'][start:end],'g', label='observation')
            if thing == 'pressure':
                ax.plot(np.array(range(xP)),self.Pressure_data['Pressure'][start:end],'g',label='observation')
            if thing == 'ice_total':
                ax.plot(np.array(range(xCDP)),self.CDP_data_ice[start:end].rolling(8).mean(),'g',label='observation')
            if thing == 'liq_total':
                ax.plot(np.array(range(xCDP)),self.CDP_data_total[start:end].rolling(8).mean(),'g',label='observation')

            ax.set_title(plot_title[i],color='#d1cfcf')
            ax.set_ylabel(ylabels[i])
            ax.legend()

        self.figure.set_facecolor('#4c4c4c')    
        # sort out subplot layout
        self.figure.tight_layout()   
        # make plot visible on canvas
        self.canvas.draw()
        # close nc file
        nc.close()

        
    def plot(self):
        self.figure.clear()
        params = {'ytick.color' : '#d1cfcf',
        'xtick.color': '#d1cfcf',
        'axes.labelcolor': '#d1cfcf',
        'axes.edgecolor':'gray'}
        plt.rcParams.update(params)

        # see how many checkboxes have been selected
        vars= ['RH','pressure','temperature','ice_total','drop_total','liquid_water_content']
        titles = ['RH','Pressure','Temperature','ice crystal number concentration',
                    'cloud droplet number concentration','LWC']
        units = ['%','hPa','K','cm$^{-3}$','cm$^{-3}$','g m$^{-3}$']
        unit_conversions = [100,1e-2,1,1e-6,1e-6,1e3]
        to_plot = [] # things to be plotted
        plot_title = []
        ylabels = []
        conversions = []
        for var, title, unit, unit_conversion in zip(vars,titles, units, unit_conversions):
            if (getattr(self,var).isChecked()):
                to_plot.append(var)
                plot_title.append(title)
                ylabels.append(unit)
                conversions.append(unit_conversion)
            else:
                pass

        # get output to plot
        try:
            nc = Dataset(self.filename)
        except:
            nc = Dataset('DEFAULT_OUTPUT.nc')
            
        # clear previously plotted stuff
        self.figure.clear()
        # plot sublots
        for i, thing in enumerate(to_plot):
            n = 100*(len(to_plot))+10+(i+1) # set-up number of subplot
            ax = self.figure.add_subplot(n, facecolor='#4c4c4c')
            ax.plot(nc[thing][:]*conversions[i],'#f3982e', label=thing)
            ax.set_title(plot_title[i],color='#d1cfcf')
            ax.set_ylabel(ylabels[i])

        self.figure.set_facecolor('#4c4c4c')    
        # sort out subplot layout
        self.figure.tight_layout()   
        # make plot visible on canvas
        self.canvas.draw()
        # close nc file
        nc.close()

    def plot_size_dis(self):
        # get output to plot
        try:
            nc = Dataset(self.filename)
        except:
            nc = Dataset('DEFAULT_OUTPUT2.nc')

        self.figure.clear()
        
        modes = nc.dimensions['modes'].size
        bins = nc.dimensions['bins'].size
        time = nc.dimensions['time'].size

        Dw = np.zeros([time,modes,bins])
        dlogDw = np.zeros([time,modes,bins-1])
        dNdlogDw = np.zeros([time,modes,bins-1])

        for j in np.arange(modes):
            for i in np.arange(time):
                Dw[i,j,:] = ((nc['liq_mass'][i,j,:]*6)/(np.pi*nc['rho'][i,j,:]))**(1/3)
                dlogDw[i,j,:] = np.diff(np.log(Dw[i,j,:]))
                dNdlogDw[i,j,:] = nc['liq_number'][i,j,1:]/dlogDw[i,j,:]

        t = int(self.plot_size_dis_time.text())
        
        if (t > time - 1):
            t = time - 1

        ax = self.figure.add_subplot(111,facecolor='#4c4c4c')
        for j in np.arange(modes):
           # ax.plot(Dw[1,j,1:]*1e6,dNdlogDw[1,j,:]/1e6,label='Mode '+str(j)+' time = 1')
           # ax.plot(Dw[int(time/2),j,1:]*1e6,dNdlogDw[int(time/2),j,:]/1e6,label='Mode '+str(j)+'time = '+str(int(time/2)))
            ax.plot(Dw[t,j,1:]*1e6,dNdlogDw[t,j,:]/1e6,label='Mode '+str(j+1)+' time = '+str(t))

        ax.set_ylabel('dNdlogD (cm-3)')
        ax.set_xlabel('D (micrometer)')
        ax.set_xscale('log')
        ax.legend()
        self.figure.set_facecolor('#4c4c4c')
        self.figure.tight_layout()
        self.canvas.draw()
    

    def plot_temp_press_fits(self):
        
        self.getTData()
        self.getPData()

        start = str(self.simulation_start_date.text())+' '+str(self.simulation_start_time.text())
        end = pd.Timestamp(start) + pd.Timedelta(seconds = float(self.simulation_duration.text()))

        self.figure.clear()

        xT = np.array(range(len(self.Temperature_data[start:end])))
        Tfit = float(self.temp_a_lineedit.text())*np.exp(-float(self.temp_b_lineedit.text())*xT) + float(self.T.text())

        xP = np.array(range(len(self.Pressure_data[start:end])))
        Pfit = float(self.press_a_lineedit.text())*np.exp(-float(self.press_b_lineedit.text())*xP) + float(self.P.text())


        ax1 = self.figure.add_subplot(211, facecolor='#4c4c4c')
        ax1.plot(xT,self.Temperature_data['Temperature'][start:end],'y')
        ax1.plot(xT,Tfit,'r')
        ax1.set_ylabel('Temperature (K)')
        
        ax2 = self.figure.add_subplot(212, facecolor='#4c4c4c')
        ax2.plot(xP,self.Pressure_data['Pressure'][start:end],'y')
        ax2.plot(xP,Pfit,'r')
        ax2.set_ylabel('Pressure (hPa)')

        self.figure.autofmt_xdate()

        self.figure.set_facecolor('#4c4c4c')
        self.figure.tight_layout()
        self.canvas.draw()

    def onParcelClicked(self):
        
        f.clearlayout(self,self.parcel_layout)
        
        self.parcel_layout = QtWidgets.QHBoxLayout()
        
        f.clearlayout(self,self.gbox_layout)

        self.updraft = QtWidgets.QLineEdit('1')
        self.updraft_label = QtWidgets.QLabel()
        self.updraft_label.setText('Upraft (ms-1)')
        self.parcel_layout.addWidget(self.updraft_label)
        self.parcel_layout.addWidget(self.updraft)

        self.simulation_type_layout.addLayout(self.parcel_layout)
        
        return self.parcel_layout

    def onChamberClicked(self):
        # clear previous chamber layout
        f.clearlayout(self,self.gbox_layout)
        
        self.description_layout = QtWidgets.QVBoxLayout()
        self.chamber_layout_temp = QtWidgets.QHBoxLayout()
        self.chamber_layout_press = QtWidgets.QHBoxLayout()
        self.chamber_layout = QtWidgets.QVBoxLayout()
        self.gbox_layout = QtWidgets.QVBoxLayout()

        # clear previous layouts from pacel button 
        f.clearlayout(self,self.parcel_layout)
        # create labels
        self.description_label1 = QtWidgets.QLabel()
        self.description_label2 = QtWidgets.QLabel()
        self.description_label1.setText('Fit parameters for temperate and pressure profile of the form')
        self.description_label2.setText('y = A*exp(-B*x) + Tint')
        self.temp_label = QtWidgets.QLabel()
        self.temp_label.setText('Temperature')
        self.press_label = QtWidgets.QLabel()
        self.press_label.setText('Pressure   ')
        self.T_a_label = QtWidgets.QLabel()
        self.T_a_label.setText('   A')
        self.T_b_label = QtWidgets.QLabel()
        self.T_b_label.setText('   B')
        self.P_a_label = QtWidgets.QLabel()
        self.P_a_label.setText('      A')
        self.P_b_label = QtWidgets.QLabel()
        self.P_b_label.setText('      B')
        self.plot_fits_button = QtWidgets.QPushButton('Plot fits')
        self.plot_fits_button.clicked.connect(self.plot_temp_press_fits)
        # create line edits
        self.temp_a_lineedit = QtWidgets.QLineEdit('6.9')
        self.temp_b_lineedit = QtWidgets.QLineEdit('0.03')
        self.press_a_lineedit = QtWidgets.QLineEdit('628.12')
        self.press_b_lineedit = QtWidgets.QLineEdit('2.6e-3')
        # add all widgets onto layouts
        self.gbox_fits = QtWidgets.QGroupBox('fit parameters')
        self.description_layout.addWidget(self.description_label1)
        self.description_layout.addWidget(self.description_label2)
        self.chamber_layout_temp.addWidget(self.temp_label)
        self.chamber_layout_temp.addWidget(self.T_a_label)
        self.chamber_layout_temp.addWidget(self.temp_a_lineedit)
        self.chamber_layout_temp.addWidget(self.T_b_label)
        self.chamber_layout_temp.addWidget(self.temp_b_lineedit)
        self.chamber_layout_press.addWidget(self.press_label)
        self.chamber_layout_press.addWidget(self.P_a_label)
        self.chamber_layout_press.addWidget(self.press_a_lineedit)
        self.chamber_layout_press.addWidget(self.P_b_label)
        self.chamber_layout_press.addWidget(self.press_b_lineedit)
        self.chamber_layout.addLayout(self.description_layout)
        self.chamber_layout.addLayout(self.chamber_layout_temp)
        self.chamber_layout.addLayout(self.chamber_layout_press)
    #   self.chamber_layout.addStretch(0)
        self.chamber_layout.addWidget(self.plot_fits_button)
        self.gbox_fits.setLayout(self.chamber_layout)
        self.gbox_layout.addWidget(self.gbox_fits)
        self.simulation_type_layout.addLayout(self.gbox_layout)
        return self.chamber_layout

    def getStyleSheet(self,path):
        f = QtCore.QFile(path)
        f.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text)
        stylesheet = QtCore.QTextStream(f).readAll()
        f.close()
        return stylesheet

    def onSaveButtonClicked(self):
        self.submit_clicked_has_happend = False
        self.filename, filter = QtWidgets.QFileDialog.getSaveFileName(parent=self, caption='Select output file', dir='.', filter='netcdf files (*.nc)')
    
    def open_new_dialog(self):
        self.nd = NewDialog(self)
        self.nd.show()

    def run(self):
        # Show the gui
        self.show()
        # Run the qt application
        qt_app.exec_()

class NewDialog(QtWidgets.QWidget):
    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self)

        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.Background,QtGui.QColor('#4c4c4c'))
        self.setPalette(palette)
        self.resize(100, 200)
        self.setWindowTitle('New Aerosol type')
        self.setStyleSheet(self.getStyleSheet('./stylesheet.qss'))

        self.dialog = QtWidgets.QWidget()

        self.layout = QtWidgets.QVBoxLayout()
        self.row1 = QtWidgets.QHBoxLayout()
        self.row2 = QtWidgets.QHBoxLayout()
        self.row3 = QtWidgets.QHBoxLayout()
        self.row4 = QtWidgets.QHBoxLayout()
        self.INP_description_layout = QtWidgets.QVBoxLayout()
        self.IASSD_layout = QtWidgets.QHBoxLayout()
        self.button_layout = QtWidgets.QHBoxLayout()
        
        self.aerosol_name_label = QtWidgets.QLabel()
        self.aerosol_name_label.setText('Aerosol name')
        self.aerosol_name = QtWidgets.QLineEdit('')
        self.molw_label = QtWidgets.QLabel()
        self.molw_label.setText('Molecular weight (kg mol<sup>-1</sup>)')
        self.molw = QtWidgets.QLineEdit()
        self.density_label = QtWidgets.QLabel()
        self.density_label.setText('Density (kg m<sup>-3</sup>)')
        self.density = QtWidgets.QLineEdit()
        self.kappa_label = QtWidgets.QLabel()
        self.kappa_label.setText('Kappa')
        self.kappa = QtWidgets.QLineEdit()
        self.INP_description = QtWidgets.QLabel()
        self.INP_description.setText('IASSD fit values of the form ns(T) = exp(a*T+b)')
        self.INP_description2 = QtWidgets.QLabel()
        urlLink=" <a href=\"https://journals.ametsoc.org/doi/abs/10.1175/JAS-D-11-0249.1\"> <font size=3 color=#d1cfcf> Niemand et al (2012)</font></a>"
        self.INP_description2.setText('following'+urlLink)
        self.INP_description2.setOpenExternalLinks(True)
        self.IASSD_a_label = QtWidgets.QLabel()
        self.IASSD_a_label.setText('a')
        self.IASSD_a = QtWidgets.QLineEdit()
        self.IASSD_b_label = QtWidgets.QLabel()
        self.IASSD_b_label.setText('b')
        self.IASSD_b = QtWidgets.QLineEdit()

        self.Enter_button = QtWidgets.QPushButton('Enter',self)
        self.Enter_button.clicked.connect(self.onEnterClicked)
        self.Cancel_button = QtWidgets.QPushButton('Cancel',self)
        self.Cancel_button.clicked.connect(self.onCancelClicked)

        self.isINP = QtWidgets.QRadioButton('Add ice nucleating properties', self)
    #   self.isINP.toggled.connect(self.onisINPClicked)
        
        self.row1.addWidget(self.aerosol_name_label)
        self.row1.addWidget(self.aerosol_name)
        self.row2.addWidget(self.molw_label)
        self.row2.addWidget(self.molw)
        self.row3.addWidget(self.density_label)
        self.row3.addWidget(self.density)
        self.row4.addWidget(self.kappa_label)
        self.row4.addWidget(self.kappa)
        self.INP_description_layout.addWidget(self.INP_description)
        self.INP_description_layout.addWidget(self.INP_description2)
        self.IASSD_layout.addWidget(self.IASSD_a_label)
        self.IASSD_layout.addWidget(self.IASSD_a)
        self.IASSD_layout.addWidget(self.IASSD_b_label)
        self.IASSD_layout.addWidget(self.IASSD_b)
        self.button_layout.addWidget(self.Enter_button)
        self.button_layout.addWidget(self.Cancel_button)

        self.layout.addLayout(self.row1)
        self.layout.addLayout(self.row2)
        self.layout.addLayout(self.row3)
        self.layout.addLayout(self.row4)
        self.layout.addWidget(self.isINP)
        self.layout.addLayout(self.INP_description_layout)
        self.layout.addLayout(self.IASSD_layout)
        self.layout.addLayout(self.button_layout)

        self.setLayout(self.layout) 

    def getStyleSheet(self,path):
        f = QtCore.QFile(path)
        f.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text)
        stylesheet = QtCore.QTextStream(f).readAll()
        f.close()
        return stylesheet

    def onCancelClicked(self):
        self.close()

    def ErrorMessage(self,string):
        self.box = QtGui.QMessageBox()
        self.box.setText(string)
        self.box.exec_()

    def onEnterClicked(self):
        # here need to write enter values to constants.py
        new_aerosol_type = list(range(4))
        
        key = self.aerosol_name.text()
        new_aerosol_type[0] = self.molw.text()
        new_aerosol_type[1] = self.density.text()
        new_aerosol_type[3] = self.kappa.text()
        
        #make temporay aerosol dictionary
        try:
            c.aerosol_dict[key] = list(map(float, new_aerosol_type))
        except:
            self.ErrorMessage('Please make sure all aerosol properties are given')
            return

        if self.isINP.isChecked():
            new_INP_type = np.zeros(3)
            new_INP_type[0] = 2
            try:
                new_INP_type[1] = self.IASSD_a.text()
            except:
                self.ErrorMessage('Please enter a value for a')
                return
            try:
                new_INP_type[2] = self.IASSD_b.text()
            except:
                self.ErrorMessage('Please enter a value for b')
                return

            c.ns_dict[key] = list(map(float, new_INP_type))
            
        with open("constants.py", "r") as in_file:
            buf = in_file.readlines()
                     
        with open("constants.py", "w") as out_file:
            for line in buf:
                if line == "aerosol_dict = {'ammonium sulphate': [132.14e-3, 1770, 3, 0.61],\n":
                    line = line + "\'"+key+"\': ["+str(new_aerosol_type[0])+','+str(new_aerosol_type[1])+','+str(new_aerosol_type[2])+','+str(new_aerosol_type[3])+"],"
                if line == "ns_dict = {'Kaolinite':[0, -0.8881, 12.9445, -0.8802, 231.383],\n":
                    if self.isINP.isChecked():
                        line = line + "\'"+key+"\': ["+str(new_INP_type[0])+','+str(new_INP_type[1])+','+str(new_INP_type[2])+"],"
                            
                out_file.write(line)


        # close window
        self.close()
                
 
# Create an instance of the application window and run it
app = main_window()
app.run()
