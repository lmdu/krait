#!/usr/bin/env python
from distutils.core import setup, Extension
setup(
	name='primerdesign',
	version='1.0',
	ext_modules=[
		Extension('primerdesign',
			sources = [
				'primerdesign_py.c',
				'thal.c',
				'oligotm.c',
				'p3_seq_lib.c',
				'libprimer3.c',
				'dpal.c',
				'primerdesign_helpers.c'
			],
			include_dirs = ['.', r'C:\winbuild\include', r'C:\winbuild\include\c++\4.8.3\ext', r'C:\winbuild\include\c++\4.8.3\backward']
		)
	]
)
