#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import apsw
import random

import config
import utils

class Row:
	def __init__(self, names, values):
		self.names = names
		self.values = values

	def getValues(self):
		return self.values

	def value(self, name):
		idx = self.names.index(name)
		return self.values[idx]

	def getKeys(self):
		return self.names

	def __iter__(self):
		for v in self.values:
			yield v

	def __nonzero__(self):
		return True if self.values else False

	def __len__(self):
		return len(self.values)

	def __getitem__(self, key):
		return self.values[key]

	def __getattr__(self, name):
		idx = self.names.index(name)
		return self.values[idx]

def row_factory(cursor, row):
	fields = [name for name, _ in cursor.getdescription()]
	return Row(fields, row)

conn = apsw.Connection(':memory:')
conn.setrowtrace(row_factory)
conn.cursor().execute("PRAGMA synchronous=0")

class Database:
	def __init__(self):
		if not self.get_tables():
			self.create_table()

	def get_size(self):
		return apsw.memoryused()

	def get_cursor(self):
		return conn.cursor()

	def create_table(self):
		self.query(config.CREATE_TABLES_SQL)

	def get_last_insert_rowid(self):
		return conn.last_insert_rowid()

	def get_fields(self, table):
		'''
		get all column names of table in sqlite, note: cursor
		getdescription must fetch row can get column name
		@para table str, table name
		@return list, column names
		'''
		cursor = self.get_cursor()
		for row in cursor.execute("SELECT * FROM %s LIMIT 1" % table):
			return [col[0] for col in cursor.getdescription()]

	def drop_tables(self):
		for table in self.get_tables():
			self.query("DROP TABLE %s" % table)
		self.query("DROP INDEX IF EXISTS loc")

	def get_tables(self):
		sql = "SELECT name FROM sqlite_master WHERE type='table'"
		return [row[0] for row in self.query(sql)]

	def get_one(self, sql):
		for row in self.query(sql):
			try:
				return row[0]
			except IndexError:
				return None

	def get_all(self, sql):
		'''
		get data from multiple rows
		'''
		return self.query(sql).fetchall()

	def get_row(self, sql):
		'''
		get data from one row
		'''
		for row in self.query(sql):
			return row

	def get_column(self, sql):
		'''
		get all values in a column
		'''
		return [row[0] for row in self.query(sql)]

	def is_empty(self, table):
		'''
		check the table is empty or not
		'''
		if self.get_one("SELECT 1 FROM %s LIMIT 1" % table):
			return False
		return True

	def get_option(self, name):
		opt = self.get_one("SELECT value FROM option WHERE name='%s' LIMIT 1" % name)
		if opt and opt.isdigit():
			return int(opt)
		return opt

	def set_option(self, name, value):
		if self.get_option(name):
			self.get_cursor().execute("UPDATE option SET value=? WHERE name=?", (value, name))
		else:
			self.get_cursor().execute("INSERT INTO option VALUES (?,?,?)", (None, name, value))

	def query(self, sql):
		return self.get_cursor().execute(sql)

	def insert(self, sql, rows):
		self.begin()
		self.get_cursor().executemany(sql, rows)
		self.commit()

	def open(self, dbfile):
		source = apsw.Connection(dbfile)
		with conn.backup("main", source, "main") as b:
			b.step()
		self.create_table()

	def save(self, dbfile):
		target = apsw.Connection(dbfile)
		return target.backup("main", conn, "main")
		#with target.backup("main", conn, "main") as b:
		#	b.step()

	def memory(self):
		return apsw.memoryused()

	def begin(self):
		self.query("BEGIN;")

	def commit(self):
		self.query("COMMIT;")

	def clear(self, table):
		self.query("DELETE FROM %s" % table)

	def changes(self):
		return conn.changes()
