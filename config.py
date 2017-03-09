#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

#cache path for plots or temp data
CACHE_PATH = os.path.join(ROOT_PATH, 'cache')

#Default sqlite3 database file
DATABASE = 'file:ssr?mode=memory&cache=shared'

#statistical json data store file
STAT_JSON = os.path.join(CACHE_PATH, 'stat.json')

#max table row display in statistical table
MAX_ROWS = 20

#create tables in database
CREATE_TABLES_SQL = """
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

