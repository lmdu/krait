#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import ctypes
#do not generate pyc file
sys.dont_write_bytecode = True

from PySide.QtCore import *
from PySide.QtGui import *

#create application and set properties
app = QApplication(sys.argv)
app.setOrganizationName('Chengdu University')
app.setOrganizationDomain('http://www.cdu.edu.cn')
app.setApplicationName('Krait')

if os.name == 'nt':
	myappid = 'CDU.Krait.ssr.0.8.3'
	ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

#set style sheet like css
with open('style.qss') as qss:
	app.setStyleSheet(qss.read())

if __name__ == '__main__':
	pixmap = QPixmap("icons/logo.png")
	splash = QSplashScreen(pixmap)
	splash.show()
	app.processEvents()
	from widgets import SSRMainWindow
	#create main windows
	win = SSRMainWindow()
	win.show()
	splash.finish(win)
	#start the main loop
	sys.exit(app.exec_())
