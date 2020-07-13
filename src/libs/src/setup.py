#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
from distutils.core import setup, Extension
from Cython.Build import cythonize

PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
PRIMER_PATH = os.path.join(PACKAGE_PATH, 'primer3')
PRIMER_SRC_PATH = os.path.join(PRIMER_PATH, 'src')
LIBPRIMER3_PATH = os.path.join(PRIMER_SRC_PATH, 'libprimer3')
THERMO_PARAMS_PATH = os.path.join(LIBPRIMER3_PATH, 'primer3_config')
KLIB_PATH = os.path.join(LIBPRIMER3_PATH, 'klib')
NCLS_PATH = os.path.join(PACKAGE_PATH, 'ncls', 'src')

libprimer3_paths = [
	os.path.join(LIBPRIMER3_PATH, 'thal.c'),
	os.path.join(LIBPRIMER3_PATH, 'oligotm.c'),
	os.path.join(LIBPRIMER3_PATH, 'p3_seq_lib.c'),
	os.path.join(LIBPRIMER3_PATH, 'libprimer3.c'),
	os.path.join(LIBPRIMER3_PATH, 'dpal.c'),
	os.path.join(PRIMER_SRC_PATH, 'primerdesign_helpers.c')
]

extensions = [
	Extension('tandem', ['tandem.c'],
		extra_compile_args=['-DMS_WIN64'],
		extra_link_args=['-DMS_WIN64']
	),
	#Extension('intersection', ['intersection.pyx']),
	#Extension('kseq', ['kseq.c'],
	#	extra_compile_args=['-DMS_WIN64'],
	#	extra_link_args=['-lz', '-DMS_WIN64']
	#),
	#Extension('primerdesign',
	#	sources=[os.path.join('primer3','src','primerdesign_py.c')] + libprimer3_paths,
	#	include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
	#	extra_compile_args=['-DMS_WIN64'],
	#	extra_link_args=['-DMS_WIN64']
	#),
	Extension('ncls', ['ncls/src/ncls.pyx', 'ncls/src/intervaldb.c'],
		include_dirs = [os.path.join(PACKAGE_PATH, "ncls/src"), PACKAGE_PATH],
		extra_compile_args=['-DMS_WIN64'],
		extra_link_args=['-DMS_WIN64']
	)
	#Extension(
	#	'thermoanalysis',
	#	sources=[os.path.join('primer3','thermoanalysis.pyx')] + libprimer3_paths,
	#	include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
	#	extra_compile_args=extra_compile_args
	#)
]

setup(
	name = 'libs',
	version = '1.0.0',
	ext_modules = cythonize(extensions, language_level=2)
)
