#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
#do not generate pyc file
sys.dont_write_bytecode = True

from PySide.QtGui import QApplication
from PySide.QtSql import QSqlDatabase, QSqlQuery
from db import open_database
import config

#create application and set properties
app = QApplication(sys.argv)
app.setOrganizationName('Chengdu University')
app.setOrganizationDomain('http://www.cdu.edu.cn')
app.setApplicationName('Niblet')

#set style sheet like css
with open('style.qss') as qss:
	app.setStyleSheet(qss.read())

#connect to sqlite database
db = QSqlDatabase.addDatabase('QSQLITE')
db.setDatabaseName(config.SSR_DB)
db.open()
QSqlQuery("PRAGMA synchronous=OFF;")

if __name__ == '__main__':
	from widgets import SSRMainWindow
	#create main windows
	win = SSRMainWindow()
	win.show()

	#start the main loop
	sys.exit(app.exec_())
