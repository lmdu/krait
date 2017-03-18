#!/usr/bin/env python
from distutils.core import setup, Extension
setup(
	name='primerdesign',
	version='1.0',
	ext_modules=[
		Extension('primerdesign',
			sources = [
				'primerdesign_py.c',
				'libprimer3/thal.c',
				'libprimer3/oligotm.c',
				'libprimer3/p3_seq_lib.c',
				'libprimer3/libprimer3.c',
				'libprimer3/dpal.c',
				'primerdesign_helpers.c'
			],
			library_dirs = [
			],
			include_dirs = [
				'libprimer3',
			]
		)
	]
)
