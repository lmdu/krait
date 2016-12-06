#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PySide.QtSql import *
from ssr import SSR


DB = QSqlDatabase.addDatabase('QSQLITE')
DB.setDatabaseName(':memory:')
DB.open()
DB.transaction()

class TableModel(object):
	table = None
	fields = []
	
	def __init__(self):
		self.query = QSqlQuery()
		self.createTable()
		self.prepareInsert()
		DB.commit()

	def __iter__(self):
		return self.fetchSSRs()

	def createTable(self):
		sql = ",".join(["%s %s" % (f, t) for f,t,_ in self.fields])
		sql = "CREATE TABLE IF NOT EXISTS %s (%s)" % (self.table, sql)
		self.query.exec_(sql)

	def getTotalCounts(self):
		self.query.exec_("SELECT COUNT(1) FROM %s LIMIT 1" % self.table)
		while self.query.next():
			return int(self.query.value(0))

	def fetchAll(self):
		self.query.exec_("SELECT * FROM %s" % self.table)
		while self.query.next():
			yield SSR({field[0]: field[2](self.query.value(idx)) for idx, field in enumerate(self.fields)})

	def prepareInsert(self):
		sql = ",".join([":%s" % f for f,_,_ in self.fields])
		sql = "INSERT INTO %s VALUES (%s)" % (self.table, sql)
		self.query.prepare(sql)

	def insert(self, ssr):
		for field in ssr:
			self.query.bindValue(":%s" % field, ssr[field])
		self.query.exec_()

class PerfectTableModel(TableModel):
	table = 'ssr'
	fields = [
		("ID", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("smotif", "TEXT", str),
		("mlength", "INTEGER", int),
		("repeat", "INTEGER", int),
		("length", "INTEGER", int)
	]
	
	def __init__(self):
		super(PerfectTableModel, self).__init__()

class CompoundTableModel(TableModel):
	table = 'cssr'
	fields = [
		("ID", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("smotif", "TEXT", str),
		("complexity", "INTEGER", int),
		("cssrs", "TEXT", str),
		("compound", "TEXT", str)
	]

	def __init__(self):
		super(CompoundTableModel, self).__init__()

#PerfectTableModel()