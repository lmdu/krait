import os
import sys

filename = os.path.join(sys.exec_prefix, 'Lib', 'distutils', 'cygwinccompiler.py')
target = 'raise ValueError("Unknown MS Compiler version %s " % msc_ver)'

with open(filename) as fp:
	lines = fp.read().replace(target, 'return ["vcruntime140"]')

with open(filename, "w") as fw:
	fw.write(lines)
