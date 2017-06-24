# -*- coding: utf-8 -*-
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

x_raw = []; t = []


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
        if q.text() == "Open Data":
            self.name = QFileDialog.getOpenFileName(self, 'Open Data','' ,"*.txt *.csv")
            
            if self.name[0]:
                data = np.genfromtxt(self.name[0], delimiter=',')
                t = data[:,0]
                x_raw = data[:,1]
    
class FormWidget(QWidget):

    def __init__(self, parent):        
        super(FormWidget, self).__init__(parent)
        self.signal_todo = []
        self.t_todo = []
        self.figure_1 = plt.figure()
        self.figure_2 = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas_1 = FigureCanvas(self.figure_1)
        self.canvas_2 = FigureCanvas(self.figure_2)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar_1 = NavigationToolbar(self.canvas_1, self)
        self.toolbar_2 = NavigationToolbar(self.canvas_2, self)

        # Just some button connected to `plot` method
        self.para_layout1 = QHBoxLayout()
        
        self.button_1 = QPushButton('Update')
        self.button_1.clicked.connect(self.plot_1)
        self.detrend_1 = QCheckBox('Detrend')
        self.detrend_1.setChecked(False)
        
        self.para_layout1.addWidget(self.detrend_1)
        self.para_layout1.addWidget(self.button_1)
        
        self.para_layout2 = QHBoxLayout()
        self.button_2 = QPushButton('Update')
        self.button_2.clicked.connect(self.plot_2)
        self.para_layout2.addWidget(self.button_2)

        # set the layout
        layout_1 = QVBoxLayout()
        layout_1.addWidget(self.toolbar_1)
        layout_1.addWidget(self.canvas_1)
        layout_1.addLayout(self.para_layout1)
        
        layout_2 = QVBoxLayout()
        layout_2.addWidget(self.toolbar_2)
        layout_2.addWidget(self.canvas_2)
        layout_2.addLayout(self.para_layout2)
        
        layout_base = QGridLayout()
        layout_base.addLayout(layout_1,1,1)
        layout_base.addLayout(layout_2,1,2)
        self.setLayout(layout_base)

    def plot_1(self):
        self.t_todo, self.signal_todo = signal_filter(t,x_raw,self.detrend_1.isChecked())
        ax = self.figure_1.add_subplot(111)
        ax.hold(False)
        ax.plot(self.t_todo,self.signal_todo, '-', color = 'k')
        ax.set_title('Signal View')
        ax.set_xlabel('Time')
        self.canvas_1.draw()
    
    def plot_2(self):
        ax = self.figure_2.add_subplot(111)
        ax.hold(False)
        ax.plot(t,x_raw, '-')
        ax.set_title('STFT')
        ax.set_xlabel('Time')
        ax.set_ylabel('Frequency')
        self.canvas_2.draw()
        



def main():
   app = QApplication(sys.argv)
   ex = Window()
   ex.show()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()