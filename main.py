#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from PySide.QtGui import QApplication
from PySide.QtSql import QSqlDatabase, QSqlQuery
from widgets import SSRMainWindow
from db import open_database

#create application and set properties
app = QApplication(sys.argv)
app.setOrganizationName('Chengdu University')
app.setOrganizationDomain('http://www.cdu.edu.cn')
app.setApplicationName('Krait')

#set style sheet like css
with open('style.qss') as qss:
	app.setStyleSheet(qss.read())

#connect to sqlite database
db = QSqlDatabase.addDatabase('QSQLITE')
db.setDatabaseName(':memory:')
db.open()
QSqlQuery("PRAGMA synchronous=OFF;")

#create main windows
win = SSRMainWindow()
win.show()

#start the main loop
sys.exit(app.exec_())
