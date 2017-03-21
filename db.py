#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import apsw
import random
from PySide.QtCore import QObject
from PySide.QtSql import QSqlDatabase, QSqlQuery

import config
import utils

conn = apsw.Connection(':memory:')

class Database:
	def __init__(self):
		if not self.get_tables():
			self.create_table()
		self.get_cursor().execute("PRAGMA synchronous=OFF;")

	def create_table(self):
		self.get_cursor().execute(config.CREATE_TABLES_SQL)

	def get_fields(self, table):
		cursor = self.get_cursor()
		for row in cursor.execute("SELECT * FROM %s LIMIT 1" % table):
			return [col[0] for col in cursor.getdescription()]

	def drop_tables(self):
		for table in self.get_tables():
			self.get_cursor().execute("DROP TABLE %s" % table)

	def get_tables(self):
		sql = "SELECT name FROM sqlite_master WHERE type='table'"
		return [row[0] for row in self.get_cursor().execute(sql)]

	def get_cursor(self):
		return conn.cursor()

	def get_one(self, sql):
		for row in self.get_cursor().execute(sql):
			if row:
				if len(row) == 1:
					return row[0]
				else:
					return row
			return None


	def get_all(self, sql):
		return self.get_cursor().execute(sql).fetchall()

	def get_columns(self, sql):
		return [row[0] for row in self.get_cursor().execute(sql)]

	def open(self, dbfile):
		source = apsw.Connection(dbfile)
		with conn.backup("main", source, "main") as b:
			b.step()
		self.create_table()

	def save(self, dbfile):
		target = apsw.Connection(dbfile)
		with target.backup("main", conn, "main") as b:
			b.step()

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
