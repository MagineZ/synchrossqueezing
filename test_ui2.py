import sys
from PyQt4.QtGui import QApplication, QDialog
from ui_imagedialog import Ui_ImageDialog

app = QApplication(sys.argv)
window = QDialog()
ui = Ui_ImageDialog()
ui.setupUi(window)

window.show()
sys.exit(app.exec_())