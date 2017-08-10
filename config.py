#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from jinja2 import Environment, FileSystemLoader

VERSION = '0.8.4'

BUILD = '20170810'

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

#cache path for plots or temp data
CACHE_PATH = os.path.join(ROOT_PATH, 'cache/')

#download fasta file directory
DOWNLOAD_PATH = os.path.join(ROOT_PATH, 'download')

#Default sqlite3 database file
DATABASE = os.path.join(CACHE_PATH, 'ssr.db')

#primer3 exetutable
PRIMER3_EXE = os.path.join(ROOT_PATH, 'primer3_core.exe')

#primer3 temp configure file
PRIMER3_SETTINGS = os.path.join(CACHE_PATH, 'primer3.conf')

#primer3 config folder
PRIMER3_CONFIG = os.path.join(ROOT_PATH, 'primer3_config/')

#statistical json data store file
STAT_JSON = os.path.join(CACHE_PATH, 'stat.json')

#max table row display in statistical table
MAX_ROWS = 20

#create jinja template reander
TEMPLATE_DIR = os.path.join(ROOT_PATH, 'template')

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
	complexity INTEGER,
	length INTEGER,
	gap INTEGER,
	component TEXT,
	structure TEXT
);

CREATE TABLE IF NOT EXISTS `vntr`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	motif TEXT,
	type INTEGER,
	repeat INTEGER,
	start INTEGER,
	end INTEGER,
	length INTEGER
);

CREATE TABLE IF NOT EXISTS `issr`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	standard TEXT,
	motif TEXT,
	type INTEGER,
	start INTEGER,
	end INTEGER,
	length INTEGER,
	match INTEGER,
	subsitution INTEGER,
	insertion INTEGER,
	deletion INTEGER,
	score INTEGER
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

CREATE TABLE IF NOT EXISTS `primer`(
	id INTEGER PRIMARY KEY,
	target TEXT,
	entry INTEGER,
	product INTEGER,
	forward TEXT,
	tm1 REAL,
	gc1 REAL,
	stability1 REAL,
	reverse TEXT,
	tm2 REAL,
	gc2 REAL,
	stability2 REAL
);

CREATE TABLE IF NOT EXISTS `primer_meta`(
	pid INTEGER,
	start1 INTEGER,
	length1 INTEGER,
	start2 INTEGER,
	length2 INTEGER
);

CREATE TABLE IF NOT EXISTS `location`(
	id INTEGER PRIMARY KEY,
	category TEXT,
	target INTEGER,
	gene_id TEXT,
	gene_name TEXT,
	feature TEXT
);

CREATE TABLE IF NOT EXISTS `option`(
	id INTEGER PRIMARY KEY,
	name TEXT,
	value TEXT
);

"""

