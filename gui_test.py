"""
Pei-Chun Su
May 2017

"""


import sys
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from functions import *
import matplotlib.pyplot as plt


class Window(QMainWindow):
    
    
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)
        """Menu Bar"""
        bar = self.menuBar()
        fileMenu = bar.addMenu('&File')
        fileMenu.addAction("Open Data")
        fileMenu.addAction("cascade")
        fileMenu.addAction("Tiled")
        fileMenu.triggered[QAction].connect(self.windowaction)        
        runMenu = bar.addMenu('Run')
        runMenu.addAction("Original Signal")
        runMenu.addAction("Time Frequency Analysis")
        runMenu.triggered[QAction].connect(self.run_file)        
        
        """arguments"""
        self.x_raw = []
        self.t = []
        self.count = 0
        """Window Title"""
        self.setWindowTitle("Synchrosquezzing")
        
        
    def windowaction(self, q):
        if q.text() == "Open Data":
            self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
            
            if self.name[0]:
                data = np.genfromtxt(self.name[0], delimiter=',')
                self.t = data[:,0]
                self.x_raw = data[:,1]
        
		
        if q.text() == "cascade":
            self.mdi.cascadeSubWindows()
		
        if q.text() == "Tiled":
            self.mdi.tileSubWindows()
        
        
    def run_file(self, q):
        if q.text() == "Original Signal":
            self.raw_signal_plot()
            
        if q.text() == "Time Frequency Analysis":
            def selectionchange(i):
                if i == 0:
                    group_box1.setEnabled(True)
                    group_box2.setEnabled(False)
                    self.rs.setEnabled(True)
                    self.normal.setChecked(True)
                    self.ConceFT_cb.setEnabled(False)
                    self.ConceFT_cb.setChecked(False)
                if i == 1:
                    group_box1.setEnabled(True)
                    group_box2.setEnabled(False)
                    self.rs.setEnabled(False)
                    self.normal.setChecked(True)
                    self.ConceFT_cb.setEnabled(False)
                    self.ConceFT_cb.setChecked(False)
                    
                if i == 2:
                    group_box1.setEnabled(False)
                    group_box2.setEnabled(False)
                    self.ConceFT_cb.setEnabled(False)
                    self.normal.setChecked(True)
                    self.ConceFT_cb.setChecked(False)
                    
            
            def btnstate(b):
                if self.cb.currentIndex() == 0:    
                    if b.text() == "Normal":
                        if b.isChecked():
                            group_box2.setEnabled(False)
                            self.ConceFT_cb.setEnabled(False)
                            self.ConceFT_cb.setChecked(False)
                            
                    if b.text() == "Synchro Squeezing":
                        if b.isChecked():
                            group_box2.setEnabled(True)
                            self.ConceFT_cb.setEnabled(True)
                    
                    if b.text() == "Reassigment":
                        if b.isChecked():
                            group_box2.setEnabled(False)
                            self.ConceFT_cb.setEnabled(True)
                    
                if self.cb.currentIndex() == 1:
                    if b.text() == "Normal":
                        if b.isChecked():
                            group_box2.setEnabled(False)
                            self.ConceFT_cb.setEnabled(False)
                            self.ConceFT_cb.setChecked(False)
                            
                    if b.text() == "Synchro Squeezing":
                        if b.isChecked():
                            group_box2.setEnabled(False)
                            self.ConceFT_cb.setEnabled(True)
                            self.ConceFT_cb.setChecked(False)
                    
            
            def ConceFT_cb_state(b):
                if b.isChecked() == True:
                    group_box4.setEnabled(True)
                else:
                    group_box4.setEnabled(False)
                
                
                
            self.dialog = QDialog()
            """######################vbox1######################################"""
    
            vbox1 = QVBoxLayout(self.dialog)
            
            """"Combo Box for selection of Time-Frequency analysis method"""
            self.cb = QComboBox()
            self.cb.addItems(["STFT", "CWT", "S-Transform"])
            self.cb.currentIndexChanged.connect(selectionchange)
            
            vbox1.addWidget(self.cb)
            
            """Group Box 1 for selection of Normal/Sq/RS"""
            group_box1 = QGroupBox("&Enhancement")  # the shortcut key is ALT + C
            
            """ radiobuttons"""
            self.normal = QRadioButton('Normal')
            self.sq = QRadioButton('Synchro Squeezing')
            self.rs = QRadioButton('Reassigment')
            self.normal.setChecked(True)
            self.normal.toggled.connect(lambda:btnstate(self.normal))
            self.sq.toggled.connect(lambda:btnstate(self.sq))
            self.rs.toggled.connect(lambda:btnstate(self.rs))
            
            gb1 = QVBoxLayout()
            gb1.addWidget(self.normal)
            gb1.addWidget(self.sq)
            gb1.addWidget(self.rs)
            
            group_box1.setLayout(gb1)
            vbox1.addWidget(group_box1)
            
            """Group Box 2 for selection of 1st/2nd"""
            group_box2 = QGroupBox("&Order")
            self.first_order = QRadioButton("1st")
            self.second_order = QRadioButton("2nd")
            self.first_order.setChecked(True)
            group_box2.setEnabled(False)
            gb2 = QVBoxLayout()
            gb2.addWidget(self.first_order)
            gb2.addWidget(self.second_order)
            group_box2.setLayout(gb2)
            vbox1.addWidget(group_box2)
            
            """Group Box 3 for general parameters"""
            group_box3 = QGroupBox('&Parameters')
            gb3 = QFormLayout()
            self.Low_freq = QLineEdit('0')
            self.Low_freq.setValidator(QDoubleValidator())
            self.High_freq = QLineEdit('0.5')
            self.High_freq.setValidator(QDoubleValidator())
            self.SamplingRate = QLineEdit('100')
            self.SamplingRate.setValidator(QIntValidator())
            self.WindowLength = QLineEdit('377')
            self.WindowLength.setValidator(QIntValidator())
            self.WindowBandwidth = QLineEdit('10') 
            self.WindowBandwidth.setValidator(QIntValidator())
            self.tDS = QLineEdit('1')
            self.tDS.setValidator(QIntValidator())
            self.FrequencyAxisResolution = QLineEdit('0.001')
            self.FrequencyAxisResolution.setValidator(QDoubleValidator())
            gb3.addRow("Low Frequency Limit",self.Low_freq)
            gb3.addRow("High_Frequency Limit",self.High_freq)
            gb3.addRow("Sampling Rate",self.SamplingRate)
            gb3.addRow("Window Length",self.WindowLength)
            gb3.addRow("Window Bandwidth",self.WindowBandwidth)
            gb3.addRow("tDS",self.tDS)
            gb3.addRow("Frequency Axis Resolution", self.FrequencyAxisResolution)
            group_box3.setLayout(gb3)
            vbox1.addWidget(group_box3)
            
            """Check Box for Multitaper"""
            self.ConceFT_cb = QCheckBox("Multitaper")
            self.ConceFT_cb.stateChanged.connect(lambda:ConceFT_cb_state(self.ConceFT_cb))
            self.ConceFT_cb.setEnabled(False)
            vbox1.addWidget(self.ConceFT_cb)
            
            """Group Box 4 for Multitaper parameters"""
            group_box4 = QGroupBox()
            gb4 = QFormLayout()
            self.NoWindowsInConceFT = QLineEdit('2')
            self.NoWindowsInConceFT.setValidator(QIntValidator())
            self.NoConceFT = QLineEdit('20')
            self.NoConceFT.setValidator(QIntValidator())
            self.Smooth = QCheckBox()
            self.Hemi = QCheckBox()
            gb4.addRow("Number of windows in ConceFT", self.NoWindowsInConceFT)
            gb4.addRow("Number of ConceFT", self.NoConceFT)
            gb4.addRow("Smooth",self.Smooth)
            gb4.addRow("Hemi",self.Hemi)
            group_box4.setLayout(gb4)
            group_box4.setEnabled(False)
            vbox1.addWidget(group_box4)
            
            
            button1 = QPushButton('Run')
            button1.clicked.connect(self.run)
            vbox1.addWidget(button1)
            
            

            # set the layout
            
            self.setLayout(vbox1)
            
            self.dialog.setWindowTitle("TF Settings & Run")
            self.dialog.setWindowModality(Qt.ApplicationModal)
            self.dialog.exec_()
            
        
    def run(self, q):
        if self.x_raw == []:
            print("no input data")
            self.dialog.close()
            return
        
        self.count = self.count+1
        widget = QWidget()
        layout = QVBoxLayout()
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        toolbar = NavigationToolbar(self.canvas, widget)
        self.plot()
        layout.addWidget(self.canvas)
        layout.addWidget(toolbar)
        
        widget.setLayout(layout)
        sub = QMdiSubWindow()
        sub.setWidget(widget)
        sub.setWindowTitle("subwindow"+str(self.count))
        sub.setGeometry(0,0,450,450)
        self.mdi.addSubWindow(sub)
        sub.show()
        self.dialog.close()

		
    def raw_signal_plot(self):
        if self.x_raw == []:
            print("no input data")
            return
        
        self.count = self.count+1
        widget = QWidget()
        layout = QVBoxLayout()
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        toolbar = NavigationToolbar(self.canvas, widget)
        toolbar.setGeometry(0,0,1,1)
        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(self.x_raw, '-')
        self.canvas.draw()
        layout.addWidget(self.canvas)
        layout.addWidget(toolbar)
        widget.setLayout(layout)
        sub = QMdiSubWindow()
        sub.setWidget(widget)
        sub.setWindowTitle("subwindow"+str(self.count))
        sub.setGeometry(0,0,450,450)
        self.mdi.addSubWindow(sub)
        sub.show()
        
        
    def plot(self):
        ''' plot some random stuff '''
        # random data
        data = self.x_raw
        time = self.t
        SamplingRate = int(self.SamplingRate.text())
        LowFrequencyLimit = float(self.Low_freq.text())
        HighFrequencyLimit = float(self.High_freq.text())
        FrequencyAxisResolution = float(self.FrequencyAxisResolution.text())
        tDS = int(self.tDS.text())
        WindowLength = int(self.WindowLength.text())
        NoWindowsInConceFT = int(self.NoWindowsInConceFT.text()) 
        WindowBandwidth = int(self.WindowBandwidth.text()) 
        NoConceFT = int(self.NoConceFT.text())
        Second = self.second_order.isChecked()
        Smooth = self.Smooth.isChecked()
        Hemi = self.Hemi.isChecked()
        z = []
        y = []
        
        class init:
            motherwavelet = 'Cinfc'
            CENTER = 1
            FWHM = 0.3
            
        if self.cb.currentIndex() == 0:
            if self.normal.isChecked():
                z,y,_,_,_ = ConceFT_sqSTFT_C(data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, tDS, WindowLength, NoWindowsInConceFT, WindowBandwidth, 0, 0, 0,0)
            elif self.sq.isChecked():
                if self.ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_sqSTFT_C(data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, tDS, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, Second, Smooth, Hemi) 
                else:
                    _, _, z, _, y = ConceFT_sqSTFT_C(data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, tDS, WindowLength, NoWindowsInConceFT, WindowBandwidth, 0, Second, 0, 0)
            elif self.rs.isChecked():
                if self.ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_rsSTFT_C(data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, tDS, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT) 
                else:
                    _, _, z, _, y = ConceFT_rsSTFT_C(data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, tDS, WindowLength, NoWindowsInConceFT, WindowBandwidth, 0) 
        
        elif self.cb.currentIndex() == 1:
            if self.normal.isChecked():
                z, _, _, _, y = ConceFT_CWT(time, data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, 0, init, 0, 0) 
            elif self.sq.isChecked():
                if self.ConceFT_cb.isChecked():
                    _, _, _, z, y = ConceFT_CWT(time, data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, NoConceFT, init, Smooth, Hemi) 
                else:
                    _, _, z, _, y = ConceFT_CWT(time, data, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, 0, init, 0, 0) 
        
        elif self.cb.currentIndex() == 2: 
            y,z = stran(data)
            
            
        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        ax.hold(False)

        # plot data
        imageSQ(ax, time, y*SamplingRate, np.abs(z), 99.5,'Greys')

        # refresh canvas
        self.canvas.draw()
        

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())