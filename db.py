#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PySide.QtCore import QObject
from PySide.QtSql import QSqlDatabase, QSqlQuery
from utils import Data

def open_database(dbname=':memory:'):
	db = QSqlDatabase.database()
	db.setDatabaseName(dbname)
	if not db.open():
		raise Exception("Can not connect to database %s" % dbname)
	query = QSqlQuery()
	query.exec_("PRAGMA synchronous=OFF;")
	return db

class SSRTable(QObject):
	table = None
	fields = []
	
	def __init__(self):
		self.query = QSqlQuery()
		self._createTable()
		self.prepareInsert()

	def __iter__(self):
		return self.fetchSSRs()

	def _createTable(self):
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
			yield Data({field[0]: field[2](self.query.value(idx)) for idx, field in enumerate(self.fields)})

	def prepareInsert(self):
		sql = ",".join([":%s" % f for f,_,_ in self.fields])
		sql = "INSERT INTO %s VALUES (%s)" % (self.table, sql)
		self.query.prepare(sql)

	def insert(self, data):
		for field in data:
			self.query.bindValue(":%s" % field, data[field])
		self.query.exec_()

class PerfectSSRTable(SSRTable):
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
		super(PerfectSSRTable, self).__init__()

class CompoundSSRTable(SSRTable):
	table = 'cssr'
	fields = [
		("ID", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("smotif", "TEXT", str),
		("complexity", "INTEGER", int),
		("length", "INTEGER", int),
		("cssrs", "TEXT", str),
		("compound", "TEXT", str)
	]

	def __init__(self):
		super(CompoundSSRTable, self).__init__()

class FastaSSRTable(SSRTable):
	table = 'fasta'
	fields = [
		("fid", "INTEGER PRIMARY KEY", int),
		("path", "TEXT", str)
	]

class SequenceSSRTable(SSRTable):
	table = 'sequence'
	fields = [
		("sid", "INTEGER PRIMARY KEY", int),
		("name", "TEXT", str),
		("fid", "INTEGER", int)
	]
