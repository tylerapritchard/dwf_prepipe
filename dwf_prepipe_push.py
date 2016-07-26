#!/usr/bin/env python3
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/
import os, time
import math
import sys
import glob
import warnings
import argparse
import subprocess

#Process New Raw fits file & Ship to g2
def dwf_prepipe_pushfile(file_name,data_dir,Qs):
	file_name=file_name.split('/')[-1].split('.')[0]

	jp2_dir=data_dir+"jp2/"
	user='fstars'
	host='g2.hpc.swin.edu.au'
	reciever=user+'@'+host
	push_dir='/lustre/projects/p025_swin/fstars/push/'
	target_dir='/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'

	print('Unpacking:'+file_name)
	print(file_name)
	print(data_dir+file_name+'.fits.fz')
	subprocess.run(['funpack',data_dir+file_name+'.fits.fz'])
	if not os.path.isdir(jp2_dir+file_name): 
		print('Creating Directory: '+jp2_dir+file_name)	
		os.makedirs(jp2_dir+file_name)
	
	print('Compressing:'+file_name)
	subprocess.run(['f2j','-i',data_dir+file_name+'.fits','-o',jp2_dir+file_name+'/'+file_name+'.jp2','Qstep='+str(Qs),'-num_threads',str(1)])
	
	print('Packaging:'+jp2_dir+file_name+'.tar')
	subprocess.run(['tar','-cf',jp2_dir+file_name+'.tar','-C',jp2_dir+file_name+'/','.'])
	
	print('Shipping:'+jp2_dir+file_name+'.tar')
	command="scp "+jp2_dir+file_name+".tar "+reciever+":"+push_dir+"; ssh "+reciever+" 'mv "+push_dir+file_name+".tar "+target_dir+"'"
	subprocess.Popen(command,shell=True)
	print('Returning to watch directory')

#Input Keyword Default Values
DWF_PID = "/home4/images/fits/2016A-0095/"
Qs_Def=0.000038

#Parse Inputs
parser = argparse.ArgumentParser(description='DWF_Prepipe push script for raw data from CTIO', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-d','--data_dir',metavar='DIRECTORY',type=str,default=DWF_PID,
	help='Directory where tarballs of compressed files are placed')
parser.add_argument('-q','--Qs',metavar='DIRECTORY',type=float,default=Qs_Def,
	help='Qstep for fits2jpeg compression')
args = parser.parse_args()

path_to_watch=args.data_dir
Qs=args.Qs

#Begin Monitoring Directory
print('Monitoring:'+path_to_watch)
before = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])
while 1:
  after = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])
  added = [f for f in after if not f in before]
  removed = [f for f in before if not f in after]
  if added: print("Added: ", ", ".join (added))
  if removed: print("Removed: ", ", ".join (removed))
  for f in added: 
  	time.sleep(10) #Make sure file finishes writing
  	dwf_prepipe_pushfile(f,path_to_watch,Qs)
  before = after
  time.sleep (2)
