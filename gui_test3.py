a# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 09:09:54 2017

@author: NTU_Math
"""

import sys
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from functions import *

x_raw = []; t = []; t_todo = []; signal_todo = []


class Window(QMainWindow):
    
    
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        """Menu Bar"""
        bar = self.menuBar()
        fileMenu = bar.addMenu('&File')
        fileMenu.addAction("Open Data")
        fileMenu.triggered[QAction].connect(self.windowaction)        
        runMenu = bar.addMenu('Run')
        runMenu.addAction("Go")
        self.form_widget = FormWidget(self) 
        self.setCentralWidget(self.form_widget) 

        """Window Title"""
        self.setWindowTitle("Time Frequency Analysis")
        
		
    def windowaction(self, q):
        global x_raw
        global t
        global t_todo
        global signal_todo
        if q.text() == "Open Data":
            self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
            
            if self.name[0]:
                data = np.genfromtxt(self.name[0], delimiter=',')
                n = len(data.shape)
                if n == 1:
                    x_raw = data; t = np.arange(0,len(data))
                    t_todo = t; signal_todo = x_raw
                    
                elif n==2:
                    t = data[:,0]; x_raw = data[:,1]
                    t_todo = t; signal_todo = x_raw
                else:
                    print ('Wrong Format')
    
class FormWidget(QWidget):

    def __init__(self, parent):        
        super(FormWidget, self).__init__(parent)
        #Figure1
        self.figure_1 = plt.figure()
        self.canvas_1 = FigureCanvas(self.figure_1)
        self.toolbar_1 = NavigationToolbar(self.canvas_1, self)
        self.ax1 = self.figure_1.add_subplot(111)
        self.ax1.set_title('Signal View')
        self.ax1.set_xlabel('Time')
        
        self.para_layout1 = QHBoxLayout()
        self.button_1 = QPushButton('Update')
        self.button_1.clicked.connect(self.update_1)
        self.detrend_1 = QCheckBox('Detrend')
        self.detrend_1.setChecked(False)
        self.para_layout1.addWidget(self.detrend_1)
        self.para_layout1.addWidget(self.button_1)
        
        #Figure2
        self.figure_2 = plt.figure()
        self.canvas_2 = FigureCanvas(self.figure_2)
        self.toolbar_2 = NavigationToolbar(self.canvas_2, self)
        self.ax2 = self.figure_2.add_subplot(111)
        self.ax2.set_title('STFT')
        self.ax2.set_xlabel('Time')
        self.ax2.set_ylabel('Frequency')

        self.para_layout2 = QHBoxLayout()
        self.zlog_cb2 = QComboBox()
        self.zlog_cb2.addItems(["normal scale","log scale"])
        self.en_cb_2 = QComboBox()
        self.en_cb_2.addItems(["Normal", "Synchrosquezzing", "Reassigment"])
        self.button_2 = QPushButton('Update')
        self.button_2.clicked.connect(self.update_2)
        self.para_layout2.addWidget(self.zlog_cb2)
        self.para_layout2.addWidget(self.en_cb_2)
        self.para_layout2.addWidget(self.button_2)
        
        #Figure 3
        self.figure_3 = plt.figure()
        self.canvas_3 = FigureCanvas(self.figure_3)
        self.toolbar_3 = NavigationToolbar(self.canvas_3, self)
        self.ax3 = self.figure_3.add_subplot(111)
        self.ax3.set_title('Continous Wavelet Transform')
        self.ax3.set_xlabel('Time')
        self.ax3.set_ylabel('Frequency')

        self.para_layout3 = QHBoxLayout()
        self.zlog_cb3 = QComboBox()
        self.zlog_cb3.addItems(["normal scale","log scale"])
        self.en_cb_3 = QComboBox()
        self.en_cb_3.addItems(["Normal", "Synchrosquezzing"])
        self.MW_cb3 = QComboBox()
        self.MW_cb3.addItems(['Cinfc','morse', 'morse-a', 'morse-b', 'morse-c', 'morlet', 'gaussian', 'meyer', 'BL3'])
        self.button_3 = QPushButton('Update')
        self.button_3.clicked.connect(self.update_3)
        self.para_layout3.addWidget(self.zlog_cb3)
        self.para_layout3.addWidget(self.en_cb_3)
        self.para_layout3.addWidget(self.MW_cb3)
        self.para_layout3.addWidget(self.button_3)

        # set the layout
        layout_1 = QVBoxLayout()
        layout_1.addWidget(self.toolbar_1)
        layout_1.addWidget(self.canvas_1)
        layout_1.addLayout(self.para_layout1)
        
        layout_2 = QVBoxLayout()
        layout_2.addWidget(self.toolbar_2)
        layout_2.addWidget(self.canvas_2)
        layout_2.addLayout(self.para_layout2)
        
        layout_3 = QVBoxLayout()
        layout_3.addWidget(self.toolbar_3)
        layout_3.addWidget(self.canvas_3)
        layout_3.addLayout(self.para_layout3)
        
        layout_base = QGridLayout()
        layout_base.addLayout(layout_1,1,1)
        layout_base.addLayout(layout_2,1,2)
        layout_base.addLayout(layout_3,2,1)
        self.setLayout(layout_base)

    def update_1(self):
        global t_todo
        global signal_todo
        if signal_todo == []:
            print ('data not loaded')
            return
        
        t_todo, signal_todo = signal_filter(t,x_raw,self.detrend_1.isChecked())
        self.ax1.hold(False)
        self.ax1.plot(t_todo,signal_todo, '-', color = 'k', linewidth = 0.3)
        self.ax1.set_title('Signal View')
        self.ax1.set_xlabel('Time')
        self.canvas_1.draw()
        
    def update_2(self):
        def plot_2():
            SR = int(SamplingRate.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            hop = int(tDS.text())
            WL = int(WindowLength.text())
            nW_InConceFT = int(NoWindowsInConceFT.text()) 
            WB = int(WindowBandwidth.text()) 
            n_ConceFT = int(NoConceFT.text())
            ifSecond = second_order.isChecked()
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            
            if self.en_cb_2.currentIndex() == 0:
                z,y,_,_,_ = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0, 0, 0,0)
            elif self.en_cb_2.currentIndex() == 1:    
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, n_ConceFT, ifSecond, ifSmooth, ifHemi) 
                else:
                    _, _, z, _, y = ConceFT_sqSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0, ifSecond, 0, 0)
            elif self.en_cb_2.currentIndex() == 2:
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_rsSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, n_ConceFT) 
                else:
                    _, _, z, _, y = ConceFT_rsSTFT_C(signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, hop, WL, nW_InConceFT, WB, 0)     
            
            self.figure_2.clear()
            self.ax2 = self.figure_2.add_subplot(111)
    
            imageSQ(self.figure_2, self.ax2, t_todo, y*SR, np.abs(z), 99.5,self.zlog_cb2.currentIndex(), 'Greys')
            self.ax2.set_title('STFT')
            self.ax2.set_xlabel('Time')
            self.ax2.set_ylabel('Frequency')
            
            self.dialog.accept()
            self.canvas_2.draw()
            
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box2.setEnabled(True)
                else:
                    group_box2.setEnabled(False)
        
        
        
        global t_todo
        global signal_todo
        
        #Confirm Data Loaded
        if signal_todo == []:
            print ('data not loaded')
            return
        
        #Create Dialog for settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit('1')
        SamplingRate.setValidator(QIntValidator())
        WindowLength = QLineEdit('377')
        WindowLength.setValidator(QIntValidator())
        WindowBandwidth = QLineEdit('10') 
        WindowBandwidth.setValidator(QIntValidator())
        tDS = QLineEdit('1')
        tDS.setValidator(QIntValidator())
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        gb1.addRow("Sampling Rate",SamplingRate)
        gb1.addRow("Window Length",WindowLength)
        gb1.addRow("Window Bandwidth",WindowBandwidth)
        gb1.addRow("tDS",tDS)
        gb1.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        group_box1.setLayout(gb1)
        
        #ConceFT Confirm
        ConceFT_cb = QCheckBox("Multitaper")
        ConceFT_cb.stateChanged.connect(lambda:ConceFT_cb_state(ConceFT_cb))
        
            
        #Multitaper box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        NoWindowsInConceFT = QLineEdit('2')
        NoWindowsInConceFT.setValidator(QIntValidator())
        NoConceFT = QLineEdit('20')
        NoConceFT.setValidator(QIntValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Number of windows in ConceFT", NoWindowsInConceFT)
        gb2.addRow("Number of ConceFT", NoConceFT)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        group_box2.setEnabled(False)
        
        #Order Box
        group_box3 = QGroupBox("&Order")
        first_order = QRadioButton("1st")
        second_order = QRadioButton("2nd")
        first_order.setChecked(True)
        gb3 = QVBoxLayout()
        gb3.addWidget(first_order)
        gb3.addWidget(second_order)
        group_box3.setLayout(gb3)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_2)
        
        #Add Widgets and Set Layout
        vbox.addWidget(group_box1)
        if self.en_cb_2.currentIndex() != 0:
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box2)
        if self.en_cb_2.currentIndex() == 1:
            vbox.addWidget(group_box3)
        vbox.addWidget(Run)
        
        #Exececute
        self.setLayout(vbox)
        self.dialog.setWindowTitle("STFT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
        
    def update_3(self):
        
        def plot_3():
            
            #Values for parameters
            SR = int(SamplingRate.text())
            LowFrequencyLimit = float(Low_freq.text())
            HighFrequencyLimit = float(High_freq.text())
            FAR = float(FrequencyAxisResolution.text())
            n_ConceFT = int(NoConceFT.text())
            ifSmooth = Smooth.isChecked()
            ifHemi = Hemi.isChecked()
            
            #Valuese for opts
            opts_in.motherwavelet = MW_type
            opts_in.CENTER = float(CENTER.text())
            opts_in.FWHM = float(FWHM.text())
            opts_in.beta = float(beta.text())
            opts_in.gam = float(gam.text())
            opts_in.dim = int(dim.text())
            if MW_type == 'morse':
                opts_in.k = int(k.text())
            elif MW_type == 'morse-a':
                opts_in.k = k_morse_a.currentIndex()
            
            
            if self.en_cb_3.currentIndex() == 0:
                z, _, _, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
            elif self.en_cb_3.currentIndex() == 1:    
                if ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, n_ConceFT, opts_in, ifSmooth, ifHemi)  
                else:
                    _, _, z, _, y = ConceFT_CWT(t_todo, signal_todo, LowFrequencyLimit, HighFrequencyLimit, FAR, 0, opts_in, 0, 0) 
            
            self.figure_3.clear()
            self.ax3 = self.figure_3.add_subplot(111)
    
            imageSQ(self.figure_3, self.ax3, t_todo, y*SR, np.abs(z), 99.5,self.zlog_cb3.currentIndex(), 'Greys')
            self.ax3.set_title('CWT')
            self.ax3.set_xlabel('Time')
            self.ax3.set_ylabel('Frequency')
            
            self.dialog.accept()
            self.canvas_3.draw()
            
            
        def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box2.setEnabled(True)
                else:
                    group_box2.setEnabled(False)
        
        
        class opts_in:
            motherwavelet = 'Cinfc'
            CENTER = 1
            FWHM = 0.3
            
        global t_todo
        global signal_todo
        
        if signal_todo == []:
            print ('data not loaded')
            return
        
        #Create Dialog for Settings
        self.dialog = QDialog()
        vbox = QVBoxLayout(self.dialog)
        
        #Group Box1 for General Parameeters
        group_box1 = QGroupBox()
        gb1 = QFormLayout()
        SamplingRate = QLineEdit('1')
        SamplingRate.setValidator(QIntValidator())
        FrequencyAxisResolution = QLineEdit('0.001')
        FrequencyAxisResolution.setValidator(QDoubleValidator())
        gb1.addRow("Sampling Rate",SamplingRate)
        gb1.addRow("Frequency Axis Resolution", FrequencyAxisResolution)
        group_box1.setLayout(gb1)
        
        #ConceFT Check Box
        ConceFT_cb = QCheckBox("Multitaper")
        ConceFT_cb.stateChanged.connect(lambda:ConceFT_cb_state(ConceFT_cb))
        
            
        #Multitaper box
        group_box2 = QGroupBox()
        gb2 = QFormLayout()
        Low_freq = QLineEdit('0')
        Low_freq.setValidator(QDoubleValidator())
        High_freq = QLineEdit('0.5')
        High_freq.setValidator(QDoubleValidator())
        NoConceFT = QLineEdit('20')
        NoConceFT.setValidator(QIntValidator())
        Smooth = QCheckBox()
        Hemi = QCheckBox()
        gb2.addRow("Low Frequency Limit",Low_freq)
        gb2.addRow("High_Frequency Limit",High_freq)
        gb2.addRow("Number of ConceFT", NoConceFT)
        gb2.addRow("Smooth",Smooth)
        gb2.addRow("Hemi",Hemi)
        group_box2.setLayout(gb2)
        group_box2.setEnabled(False)
        
        #Mother Wavelet Group Box
        MW_type = self.MW_cb3.currentText()
        MW_group_box = QGroupBox('Mother Wavelet')
        
        MW_gb = QFormLayout()
        l1 = QLabel()
        l1.setText(MW_type)
        CENTER = QLineEdit('1')
        CENTER.setValidator(QDoubleValidator())
        FWHM = QLineEdit('0.3')
        FWHM.setValidator(QDoubleValidator())
        dim = QLineEdit('1')
        dim.setValidator(QIntValidator())
        beta = QLineEdit('1')
        beta.setValidator(QDoubleValidator())
        gam = QLineEdit('1')
        gam.setValidator(QDoubleValidator())
        k = QLineEdit('0')
        k.setValidator(QIntValidator())
        k_morse_a = QComboBox()
        k_morse_a.addItems(['0','1'])
        
        MW_gb.addWidget(l1)
        if (MW_type == 'Cinfc') or (MW_type == 'gaussian'):
            MW_gb.addRow('Center', CENTER)
            MW_gb.addRow('FWHM',FWHM)
        elif MW_type in ['morse', 'morse-a', 'morse-b','morse-c'] :
            MW_gb.addRow('beta',beta)
            MW_gb.addRow('gam',gam)
            if MW_type in ['morse-b', 'morse-c']:
                MW_gb.addRow('dim',dim)
            elif MW_type == 'morse':
                MW_gb.addRow('k',k)
            elif MW_type == 'morse-a':
                MW_gb.addRow('k',k_morse_a)
        MW_group_box.setLayout(MW_gb)
        
        #Run Button
        Run = QPushButton('Run')
        Run.clicked.connect(plot_3)
        
        vbox.addWidget(MW_group_box)
        vbox.addWidget(group_box1)
        if self.en_cb_3.currentIndex() == 1:
            vbox.addWidget(ConceFT_cb)
            vbox.addWidget(group_box2)
        vbox.addWidget(Run)
        
        
        self.setLayout(vbox)
        self.dialog.setWindowTitle("CWT Settings")
        self.dialog.setWindowModality(Qt.ApplicationModal)
        self.dialog.exec_()
        
    
        
        



def main():
   app = QApplication(sys.argv)
   ex = Window()
   ex.show()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()