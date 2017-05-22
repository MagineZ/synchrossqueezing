import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

import random

class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        hbox = QHBoxLayout()
        vbox1 = QVBoxLayout()
        
        self.cb = QComboBox()
        self.cb.addItems(["STFT", "CWT", "STRAN"])
        self.cb.currentIndexChanged.connect(self.selectionchange)
        
        vbox1.addWidget(self.cb)
        
        self.b1 = QRadioButton("Button1")
        self.b1.setChecked(True)
        self.b1.toggled.connect(lambda:self.btnstate(self.b1))
        vbox1.addWidget(self.b1)
        self.b2 = QRadioButton("Button2")
        self.b2.toggled.connect(lambda:self.btnstate(self.b2))

        vbox1.addWidget(self.b2)

   
        self.button1 = QPushButton('Run')
        vbox1.addWidget(self.button1)
        
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # set the layout
        vbox = QVBoxLayout()
        vbox.addWidget(self.toolbar)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.button)
        
        hbox.addLayout(vbox1)
        hbox.addLayout(vbox)
        self.setLayout(hbox)

    
    def selectionchange(self,i):
        print ("Items in the list are :")
		
        for count in range(self.cb.count()):
            print (self.cb.itemText(count))
        print ("Current index",i,"selection changed ",self.cb.currentText())
    
    def plot(self):
        ''' plot some random stuff '''
        # random data
        data = [random.random() for i in range(10)]

        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(data, '*-')

        # refresh canvas
        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())