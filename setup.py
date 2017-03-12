#!/usr/bin/env python
from distutils.core import setup, Extension
setup(name='tandem', version='1.0', ext_modules=[Extension('tandem', ['tandem.c'])])
