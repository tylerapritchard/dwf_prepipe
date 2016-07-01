#!/usr/bin/env python3

import os, time
import math
import sys
import glob
import argparse
import warnings
import multiprocessing
import subprocess
import astropy.io.fits as pyfits

Test_Data_Repository='/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/test_data/'
Push_Directory='/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'

file_list = os.listdir (Test_Data_Repository)
n=0
wtime=30
for f in file_list:
	print('Transfering File: '+f)
	subprocess.run(['cp',Test_Data_Repository+f,Push_Directory])
	print('Waiting '+str(wtime)+' seconds')
	time.sleep(wtime)