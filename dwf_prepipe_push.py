#!/usr/bin/env python3
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/

import os, time
import math
import sys
import glob
import warnings
import multiprocessing
import subprocess
import astropy.io.fits as pyfits

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
	#subprocess.Popen(['scp',jp2_dir+file_name+'.tar',reciever+':'+push_dir,';','ssh',reciever,"'mv "+push_dir+file_name+".tar "+target_dir+"'"])
	command="scp "+jp2_dir+file_name+".tar "+reciever+":"+push_dir+"; ssh "+reciever+" 'mv "+push_dir+file_name+".tar "+target_dir+"'"
	#print(command)
	subprocess.Popen(command,shell=True)
	print('Returning to watch directory')

path_to_watch = "/home4/images/fits/2016A-0095/"
Qs=0.0000001#0.000038
print('Monitoring:'+path_to_watch)
before = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])
while 1:
  after = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])
  added = [f for f in after if not f in before]
  removed = [f for f in before if not f in after]
  if added: print("Added: ", ", ".join (added))
  if removed: print("Removed: ", ", ".join (removed))
  for f in added: 
  	time.sleep(10 )
  	dwf_prepipe_pushfile(f,path_to_watch,Qs)
  before = after
  time.sleep (2)
