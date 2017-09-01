#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time
#do not generate pyc file
sys.dont_write_bytecode = True

import krait_rc
import config

from PySide.QtCore import *
from PySide.QtGui import *

#create application and set properties
app = QApplication(sys.argv)
app.setOrganizationName('Chengdu University')
app.setOrganizationDomain('http://www.cdu.edu.cn')
app.setApplicationName('Krait')

if os.name == 'nt':
	import ctypes
	myappid = 'CDU.Krait.ssr.%s' % config.VERSION
	ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

#set style sheet like css
app.setStyleSheet(config.STYLE_QSS)

splash_img = QPixmap("icons/splash.png")
splash = QSplashScreen(splash_img)
splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.SplashScreen | Qt.FramelessWindowHint)
splash.setEnabled(False)

probar = QProgressBar(splash)
probar.setStyleSheet('''
QProgressBar{
	border:none;
	width:490px;
	height: 10px;
	text-align: center;
	padding: 0 5px;
}
QProgressBar::chunk{
	background: white;
	border-radius: 7px;
}
''')
probar.setMaximum(10)
probar.setGeometry(0, splash_img.height()-30, splash_img.width(), 10)


splash.show()
for i in range(10):
	time.sleep(1)
	probar.setValue(i)
	splash.showMessage('print for %s' % i, Qt.AlignCenter | Qt.AlignBottom, Qt.white)
	app.processEvents()

from widgets import SSRMainWindow
#create main windows
win = SSRMainWindow()
win.show()

splash.finish(win)
#start the main loop
sys.exit(app.exec_())
