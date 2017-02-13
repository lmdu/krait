#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import psutil
import subprocess

def format_boulder_io(boulder_list):
	records = []
	for boulder in boulder_list:
		boulder = ''.join(['{}={}\n'.format(k,v) for k,v in boulder.items()])
		boulder += '=\n'
		records.append(boulder)
	return ''.join(records)

def parse_boulder_io(boulder_lines):
	res = {}
	for line in boulder_lines:
		k,v = line.strip().split('=')
		res[k] = v
	return res

class PrimerDesign:
	'''
	Invoke primer3_core to design primers for ssr with flanking sequences
	@para setting_file, the file path of p3_setting_file used by primer3
	'''
	def __init__(self, setting_file, boulder_list):
		self.setting_file = setting_file
		self.boulder_list = boulder_list

	def __iter__(self):
		return self.designer()

	@property 
	def cmd(self):
		command = ['primer3_core']
		if self.setting_file is not None:
			command.append('-p3_settings_file=%s' % self.setting_file)

	def designer(self):
		proc = psutil.Popen(self.cmd, 
			stdout=subprocess.PIPE,
			stdin=subprocess.PIPE, 
			stderr=subprocess.STDOUT
		)
		boulder_str = format_boulder_io(self.boulder_list)
		proc.stdin.write(boulder_str)
		proc.stdin.flush()
		record = []
		while 1:
			line = proc.stdout.readline()
			if not line:
				if proc.poll() is not None:
					break
				else:
					time.sleep(0.01)

			if line.strip() == '=':
				yield parse_boulder_io(record)
				record = []
			else:
				record.append(line.strip())


