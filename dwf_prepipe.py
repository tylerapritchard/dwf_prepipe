#!/usr/bin/env python 

import os, time
import sys
import glob
import argparse
import warnings
import multiprocessing

def dwf_prepipe_unpack(file_name,push_path,untar_path,processes):
	#placeholder
	#untar file: tar -xf DECam_00504634.tar -C untar/
	#parrallel qsub
	pool = multiprocessing.Pool(processes=processes)
	pool.map(dwf_prepipe_qsubccd,ccd_files)
	pool.close()
	pool.join()

def dwf_prepipe_qsubccd():
	#placeholder

def dwf_prepipe_processccd(file_name):
	#Uncompress CCD
	#Copy CCD to Armins Directory with Proper name
	#Clean up files
	#Call Pipeloop
	#Diagnose & Fix Pipeloop

DWF_Push = '/lustre/projects/p025_swin/tap/DWF_Unpack_Test/push/' #"/lustre/projects/p025_swin/pipes/DWF_PIPE/CTIO_PUSH/"

parser = argparse.ArgumentParser(description='Handle File Ingests for the DWF pipeline', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
		help='Directory where tarballs of compressed files are placed')
	parser.add_argument('-p', dest='processes', type=int, default=multiprocessing.cpu_count(),
		help='Number of processes to plot with (default: #CPU cores)')
	args = parser.parse_args()

path_to_watch = args.push_dir
path_to_untar = args.push_dir,join('untar/')
before = dict ([(f, None) for f in os.listdir (path_to_watch.join('*.tar'))])

while 1:
  time.sleep (10)
  after = dict ([(f, None) for f in os.listdir (path_to_watch.join('*.tar'))])
  added = [f for f in after if not f in before]
  removed = [f for f in before if not f in after]
  if added: print "Added: ", ", ".join (added)
  if removed: print "Removed: ", ", ".join (removed)

  for f in added:  dwf_prepipe_unpack(f,path_to_watch,path_to_untar,args.processes)

  before = after



