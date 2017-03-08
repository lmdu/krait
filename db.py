#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import apsw
from PySide.QtCore import QObject
from PySide.QtSql import QSqlDatabase, QSqlQuery

import config
import utils

sql = """
CREATE TABLE IF NOT EXISTS `ssr`(
	sid INTEGER PRIMARY KEY,
	sequence TEXT,
	start INTEGER,
	stop INTEGER,
	motif TEXT,
	standard TEXT,
	type INTEGER,
	repeat INTEGER,
	length INTEGER
);
CREATE TABLE IF NOT EXISTS `cssr`(
	cid INTEGER PRIMARY KEY
	sequence TEXT,
	start INTEGER,
	stop INTEGER,
	motif TEXT,
	standard TEXT,
	complexity INTEGER,
	length INTEGER,
	component TEXT,
	structure TEXT
);

CREATE TABLE IF NOT EXISTS `ltr`(
	lid INTEGER PRIMARY KEY,
	sequence TEXT,
	start INTEGER,
	stop INTEGER,
	motif TEXT,
	type INTEGER,
	repeat INTEGER,
	length INTEGER
);

CREATE TABLE IF NOT EXISTS `issr`(
	iid INTEGER PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS `fasta`(
	fid INTEGER PRIMARY KEY,
	path TEXT
);

CREATE TABLE IF NOT EXISTS `sequence`(
	sid INTEGER PRIMARY KEY,
	name TEXT,
	fid INTEGER
);

CREATE TABLE IF NOT EXISTS `meta`(
	name, TEXT,
	value, TEXT
);

"""


class DB:
	_conn = None
	def __init__(self):
		if self._conn is None:
			self._conn = apsw.Connection(config.DATABASE)

	def getCursor(self):
		return self._conn.cursor()

	def __del__(self):
		self._conn.close()








__all__ = [
	"open_database",
	"MicrosatelliteTable",
	"CompoundTable",
	"SatelliteTable",
	"FastaTable",
	"SequenceTable",
	"MetaTable"
]

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
		self._createTable()
		self.prepareInsert()

	def __iter__(self):
		return self.fetchSSRs()

	def _createTable(self):
		sql = ",".join(["%s %s" % (f, t) for f,t,_ in self.fields])
		sql = "CREATE TABLE IF NOT EXISTS %s (%s)" % (self.table, sql)
		QSqlQuery(sql)

	def recordCounts(self):
		'''
		get the number of records in the table
		@return int, counts of records
		'''
		return self.get("SELECT COUNT(1) FROM %s LIMIT 1" % self.table)


	def fetchAll(self):
		'''
		get all the records in the table and return one in each time
		'''
		query = QSqlQuery("SELECT * FROM %s" % self.table)
		while query.next():
			yield Data({field[0]: field[2](query.value(idx)) for idx, field in enumerate(self.fields)})

	def query(self, sql):
		query = QSqlQuery(sql)
		rec = query.record()
		fields = {rec.fieldName(i): i for i in range(rec.count())}
		while query.next():
			yield Data({field: query.value(fields[field]) for field in fields})

	def get(self, sql):
		query = QSqlQuery(sql)
		while query.next():
			return query.value(0)

	def prepareInsert(self):
		'''
		prepare insert data into table
		'''
		sql = ",".join([":%s" % f for f,_,_ in self.fields])
		sql = "INSERT INTO %s VALUES (%s)" % (self.table, sql)
		self._query = QSqlQuery()
		self._query.prepare(sql)

	def insert(self, data):
		'''
		insert data into table
		@para data dict, a record with fields and values
		'''
		for field in data:
			self._query.bindValue(":%s" % field, data[field])
		self._query.exec_()


class MicrosatelliteTable(SSRTable):
	table = 'ssr'
	fields = [
		("mid", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("standard", "TEXT", str),
		("type", "INTEGER", int),
		("repeat", "INTEGER", int),
		("length", "INTEGER", int)
	]
	
	def __init__(self):
		super(MicrosatelliteTable, self).__init__()

class CompoundTable(SSRTable):
	table = 'cssr'
	fields = [
		("cid", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("standard", "TEXT", str),
		("complexity", "INTEGER", int),
		("length", "INTEGER", int),
		("component", "TEXT", str),
		("structure", "TEXT", str)
	]

	def __init__(self):
		super(CompoundTable, self).__init__()

class SatelliteTable(SSRTable):
	table = 'satellite'
	fields = [
		("sid", "INTEGER PRIMARY KEY", int),
		("sequence", "TEXT", str),
		("start", "INTEGER", int),
		("stop", "INTEGER", int),
		("motif", "TEXT", str),
		("type", "INTEGER", int),
		("repeat", "INTEGER", int),
		("length", "INTEGER", int)
	]

class FastaTable(SSRTable):
	table = 'fasta'
	fields = [
		("fid", "INTEGER PRIMARY KEY", int),
		("path", "TEXT", str)
	]

class SequenceTable(SSRTable):
	table = 'sequence'
	fields = [
		("sid", "INTEGER PRIMARY KEY", int),
		("name", "TEXT", str),
		("fid", "INTEGER", int)
	]

class MetaTable(SSRTable):
	table = 'meta'
	fields = [
		("name", "TEXT", str),
		("value", "TEXT", str)
	]

	def __init__(self):
		super(MetaTable, self).__init__()

	def getMeta(self, name):
		'''
		get the meta information 
		@para name str, meta name
		@return str
		'''
		self._query.exec_("SELECT value FROM meta WHERE name='%s' LIMIT 1" % name)
		while self._query.next():
			return self._query.value(0)
