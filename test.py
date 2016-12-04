#!/usr/bin/env python
import sqlite3

'''
conn = sqlite3.connect(':memory:')
cursor = conn.cursor()
cursor.execute("CREATE TABLE ssr (ID INTEGER PRIMARY KEY, code TEXT)")
cursor.execute("INSERT INTO ssr VALUES (?,?)", (None, 123))
cursor.execute("INSERT INTO ssr VALUES (?,?)", (None, 5678))
cursor.execute("ATTACH DATABASE 'test.db' AS 'fdb'")
cursor.execute("CREATE TABLE fdb.ssr AS SELECT * FROM ssr")
cursor.close()
conn.commit()
conn.close()
'''
class A:
	a = 1
	b =2

class B:
	def __init__(self):
		self.c = A()

	def __enter__(self):
		return self.c
