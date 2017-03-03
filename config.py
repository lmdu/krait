#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

#cache path for plots or temp data
CACHE_PATH = os.path.join(ROOT_PATH, 'cache')

#Default sqlite3 database file
SSR_DB = os.path.join(CACHE_PATH, 'ssr.db')

#statistical json data store file
STAT_JSON = os.path.join(CACHE_PATH, 'stat.json')

#max table row display in statistical table
MAX_ROWS = 20

