#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

#cache path for plots or temp data
CACHE_PATH = os.path.join(ROOT_PATH, 'cache')

#Default sqlite3 database file
DATABASE = os.path.join(CACHE_PATH, 'ssr.db')

#primer3 exetutable
PRIMER3_EXE = os.path.join(ROOT_PATH, 'primer3_core.exe')

#primer3 temp configure file
PRIMER3_SETTINGS = os.path.join(CACHE_PATH, 'primer3.conf')

#statistical json data store file
STAT_JSON = os.path.join(CACHE_PATH, 'stat.json')

#max table row display in statistical table
MAX_ROWS = 20

#create tables in database
CREATE_TABLES_SQL = """
CREATE TABLE IF NOT EXISTS `ssr`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	standard TEXT,
	motif TEXT,
	type INTEGER,
	repeat INTEGER,
	start INTEGER,
	end INTEGER,
	length INTEGER
);
CREATE TABLE IF NOT EXISTS `cssr`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	start INTEGER,
	end INTEGER,
	motif TEXT,
	standard TEXT,
	complexity INTEGER,
	length INTEGER,
	component TEXT,
	structure TEXT
);

CREATE TABLE IF NOT EXISTS `vntr`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	start INTEGER,
	end INTEGER,
	motif TEXT,
	type INTEGER,
	repeat INTEGER,
	length INTEGER
);

CREATE TABLE IF NOT EXISTS `issr`(
	id INTEGER PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS `fasta`(
	id INTEGER PRIMARY KEY,
	path TEXT
);

CREATE TABLE IF NOT EXISTS `seq`(
	id INTEGER PRIMARY KEY,
	name TEXT,
	fid INTEGER
);

CREATE TABLE IF NOT EXISTS `meta`(
	name TEXT,
	value TEXT
);

"""

