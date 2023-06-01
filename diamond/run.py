#!/usr/bin/env python

import sh
import sys
import fileinput


class Handler(object):
	def __init__(self):
		self.cmd = sh.Command('/usr/local/bin/root')
		self.cmd = self.cmd.bake('-l')
		self.cfile = 'pypad2D.C'


	def modifyCode(self, l):
		for item in l:
			for i, line in enumerate(fileinput.input(self.cfile, inplace=1)):
				sys.stdout.write(item[1] if item[0] in line else line)


	def runCode(self):
		self.cmd("pycmd")


def main():
	h = Handler()

	for y in range(50, 499, 50):
		ypos = y
		repl = [
		["Int_t entryPointY", "Int_t entryPointY = " + str(ypos) + ";\n"],
		["const char *currentProfile=", 'const char *currentProfile="../Results/prof_y=' + str(ypos) + '.dat";\n']
		]
		h.modifyCode(repl)
		h.runCode()
		print('y=', y)















if __name__ == '__main__':
	main()