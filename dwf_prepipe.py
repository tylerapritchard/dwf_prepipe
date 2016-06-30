#!/usr/bin/env python 

import os, time
import sys
import glob
import argparse
import warnings
import multiprocessing
import subprocess
import pyfits

def dwf_prepipe_unpack(file_name,push_path,untar_path,qsub_path,processes):
	subprocess.run(['tar','-xvf',push_path+file_name,'-C',untar_path])

	DECam_Root=file_name.split('.')[0]
	Exposure=DECam_Root.split('_')[1]

	infoname=untar_path+DECam_Root+'_fileinfo.fits'
	j2fpath='/lustre/projects/p025_swin/dvohl/kerlumph_7_4_decam/bin/Linux-x86-64-gcc/'
	subprocess.run(['j2f','-i',untar_path+DECam_Root+'_1.jp2','-o',infoname])
	#subprocess.run(j2f+'j2f -i '+untar_path+DECam_Root+'_1.jp2 -o '+infoname, shell=True)
	exp=pyfits.getval( filename+'[1]',"EXPNUM")
	Field=pyfits.getval(infoname+'[1]',"OBJECT")
	Filter=pyfits.getval( filename+'[1]',"FILTER")[0]
	ut='ut'+pyfits.getval( filename+'[1]',"OBSID")[6:12]
	mjd=pyfits.getval(filename+'[1]',"MJD-OBS")
	date=pyfits.getval(filename+'[1]',"DATE-OBS")
	print(Field+'.'+Filter+'.'+UT+'.'+EXP+'.'+UT+'_')
	#Filename_Root=
	#placeholder
	#untar file: tar -xf DECam_00504634.tar -C untar/
	#parrallel ccd script write
	#pool = multiprocessing.Pool(processes=processes)
	#pool.map(dwf_prepipe_qsubccd,ccd_files)
	#pool.close()
	#pool.join()

def dwf_prepipe_qsubccd():
	a='placeholder'

	#placeholder

def dwf_prepipe_processccd(file_name):
	a='placeholder'
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
path_to_untar = args.push_dir+'untar/'
path_to_qsub = args.push_dir+'qsub/'
before = dict ([(f, None) for f in os.listdir (path_to_watch)])

while 1:
  time.sleep (10)
  after = dict ([(f, None) for f in os.listdir (path_to_watch)])
  added = [f for f in after if not f in before]
  removed = [f for f in before if not f in after]
  if added: print("Added: ", ", ".join (added))
  if removed: print("Removed: ", ", ".join (removed))
  for f in added:  
  	print('Unpacking: '+f)
  	dwf_prepipe_unpack(f,path_to_watch,path_to_untar,path_to_qsub,args.processes)

  before = after



