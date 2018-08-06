#!/usr/bin/env python3
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/

import os, time
import math
import sys
import glob
import argparse
import warnings
import multiprocessing
import subprocess
import astropy.io.fits as pyfits
import numpy as np
import re


def main():
	command='qstat -u fstars'
	proc=subprocess.Popen(['qstat','-u','fstars'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = proc.communicate()

	count=0
	for line in stdout.split(os.linesep):
		job=line.split()
		print(job)
		if(len(job) > 0):
			if(re.match('\d',job[0])):
				if(re.match('hsc*',job[3])):
					if(job[9] == 'R'):
						count=count+1
	print(count)


if __name__ == '__main__':
    main()	