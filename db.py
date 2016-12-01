#!/usr/bin/env python
try:
	import sqlite3
except ImportError:
	from pysqlite2 import dbapi2 as sqlite3

#create a SSR object from dict
class Row(sqlite3.Row):
	def __getattr__(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

#connect to ssr database and create ssr table
class SSRDB:
	def __init__(self, db):
		self.conn = sqlite3.connect(db)
		self.conn.row_factory = Row
		self.cursor = self.conn.cursor()

	def createTable(self, sql):
		c = self.conn.cursor()
		try:
			c.execute(sql)
		finally:
			c.close()

	def createSSRTable(self):
		self.createTable('''
			CREATE TABLE IF NOT EXISTS ssr(
				ID INTEGER PRIMARY KEY,
				chrom TEXT,
				start INTEGER,
				end INTEGER,
				repeat INTEGER,
				length INTEGER,
				motif INTEGER,
				smotif INTEGER
			)
		''')

	def isTableExists(self, name):
		sql = (
			"SELECT * FROM sqlite_master"
			" WHERE name=? and type='table'"
		)
		return self.get(sql, name)
			
	def iter(self, sql):
		c = self.conn.cursor()
		try:
			for row in c.execute(sql):
				yield row
		finally:
			c.close()

	def get(self, sql, *para):
		c = self.conn.cursor()
		try:
			c.execute(sql, para)
			row = c.fetchone()
			return row[0] if row else None
		finally:
			c.close()

	def insert(self, sql, data):
		self.cursor.execute(sql, data)

	def submit(self):
		self.cursor.close()
		self.conn.commit()