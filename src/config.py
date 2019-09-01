#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import appdirs
from jinja2 import Environment, FileSystemLoader

VERSION = "0.11.6"

BUILD = '20190901'

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

APPDATA_PATH = appdirs.AppDirs('Krait')

#config files
CONFIG_FILE = os.path.join(APPDATA_PATH.user_data_dir, 'config.ini')

#cache path for plots or temp data
CACHE_PATH = "%s/" % APPDATA_PATH.user_data_dir

#primer3 config folder
PRIMER3_CONFIG = os.path.join(ROOT_PATH, 'primer3_config/')

#max table row display in statistical table
MAX_ROWS = 20

#create jinja template reander
TEMPLATE_DIR = os.path.join(ROOT_PATH, 'template/')

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
	fid INTEGER,
	size INTEGER,
	gc INTEGER,
	ns INTEGER
);

CREATE TABLE IF NOT EXISTS `meta`(
	name TEXT,
	value TEXT
);

CREATE TABLE IF NOT EXISTS `primer`(
	id INTEGER PRIMARY KEY,
	category TEXT,
	target INTEGER,
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

CREATE TABLE IF NOT EXISTS `gene`(
	id INTEGER PRIMARY KEY,
	sequence TEXT,
	start INTEGER,
	end INTEGER,
	geneid TEXT,
	genename TEXT,
	biotype TEXT
);

CREATE TABLE IF NOT EXISTS `location`(
	id INTEGER PRIMARY KEY,
	reptype INTEGER,
	target INTEGER,
	feature INTEGER,
	gid INTEGER
);

CREATE TABLE IF NOT EXISTS `option`(
	id INTEGER PRIMARY KEY,
	name TEXT,
	value TEXT
);

"""

STYLE_QSS = """
*{
	font-family: Roboto;
	font-size: 15px;
}

/* main windows */
SSRMainWindow{
	background:#fff;
}

QTableView{
	border: 0;
	selection-background-color: #F6F6F6;
	selection-color: #000000;
}

QHeaderView::section {
    background-color: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                      stop: 0 #fefefe, stop: 0.5 #fbfbfb,
                                      stop: 0.6 #f8f8f8, stop: 1 #f2f3f3);
    padding: 5px 0;
    border:none;
    border-right: 1px solid #d8d8d8;
    border-bottom: 1px solid #afb0b1;
}

QTextBrowser{
	border: 0;
}

/* tool bar */
QToolBar{
	background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #ffffff, stop: 1 #f2f2f2);
	border-bottom:1px solid #a6a6a6;
	spacing: 8px;
}
SSRFilterInput{
	border:1px solid #a9a9a9;
	padding:4px 18px 4px 4px;
	border-radius: 2px;
	margin-left:50px;
	margin-right:5px;
	background: #fff url(:/icons/filter.png);
	background-position: right center;
	background-repeat: none;
	width:100%;
}

/* status bar */
QStatusBar{
	background:qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #ffffff, stop: 1 #ededed);
	border-top:1px solid #ccc;
}
QProgressBar{
	text-align: right;
	max-width: 120px;
	min-width: 120px;
	max-height: 15px;
	min-height: 15px;
}
"""

