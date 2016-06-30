#!/usr/bin/env python3
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/

import os, time
import sys
import glob
import argparse
import warnings
import multiprocessing
import subprocess
import astropy.io.fits as pyfits

def dwf_prepipe_unpack(file_name,push_path,untar_path,qsub_path,processes):
	ccdlist=['38','24','14','56','7','33','11','5','34','58','10','41','16','54','17','18','27','2','3','30','22','23','21','43','6','37','25','19','36','44','32','26','4','15','29','50','12','13','1','42','48','52','55','31','20','59','39','57','53','9','28','35','8','46','45','49','40','51','47']
	DECam_Root=file_name.split('.')[0]

	print('Unpacking:'+file_name)
	subprocess.run(['tar','-xf',push_path+file_name,'-C',untar_path])

	Exposure=DECam_Root.split('_')[1]

	infoname=untar_path+DECam_Root+'_fileinfo.fits'
	subprocess.run(['j2f','-i',untar_path+DECam_Root+'_1.jp2','-o',infoname])

	exp=pyfits.getval(infoname,"EXPNUM")
	Field=pyfits.getval(infoname,"OBJECT")
	Filter=pyfits.getval(infoname,"FILTER")[0]
	ut='ut'+pyfits.getval(infoname,"OBSID")[6:12]
	mjd=pyfits.getval(infoname,"MJD-OBS")
	date=pyfits.getval(infoname,"DATE-OBS")

	process_root=Field+'.'+Filter+'.'+ut+'.'+str(exp)
	print('Processing:'+process_root)

	#parrallel ccd script write
	ccd_files=process_root.join(ccdlist)
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

DWF_Push = '/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/' #"/lustre/projects/p025_swin/pipes/DWF_PIPE/CTIO_PUSH/"

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


dwf_prepipe_unpack('DECam_00504110.tar',path_to_watch,path_to_untar,path_to_qsub,args.processes)

#while 1:
#  after = dict ([(f, None) for f in os.listdir (path_to_watch)])
#  added = [f for f in after if not f in before]
#  removed = [f for f in before if not f in after]
#  if added: print("Added: ", ", ".join (added))
#  if removed: print("Removed: ", ", ".join (removed))
#  for f in added:  
#  	dwf_prepipe_unpack(f,path_to_watch,path_to_untar,path_to_qsub,args.processes)

#  before = after
#  time.sleep (5)




