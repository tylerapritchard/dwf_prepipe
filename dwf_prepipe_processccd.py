#!/usr/bin/env python3

import os
import datetime
import time
import shutil
import argparse
import subprocess

#~/.astropy/config/astropy.cfg was getting messed up - seperate default (used by pipeloop?) and this
os.environ['XDG_CONFIG_HOME']='/home/fstars/.python3_config/'

import astropy.io.fits as pyfits



DWF_Push = '/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'
use_local=1

parser = argparse.ArgumentParser(description='DWF_Prepipe process for single compressed CCD', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-p','--push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
	help='Directory where taballs of compressed files are placed')
parser.add_argument('-l','--local',metavar='DIRECTORY',type=int,default=use_local,
	help='Use node local storage for jpg to fits conversion? (1 or 0)')
parser.add_argument('-i','--input_file',type=str,help='input .jp2 file')

args = parser.parse_args()

#Set local Directory and check to see if it exists
local_dir='/lfs/data0/dwf/'
local_convert=args.local

if(local_convert):
	if not os.path.isdir(local_dir): 
		print('Creating Directory: '+local_dir)	
		os.makedirs(local_dir)



photepipe_rawdir= '/projects/p025_swin/pipes/arest/DECAM/DEFAULT/rawdata/'
push_dir=args.push_dir
untar_path=push_dir+'untar/'

file_name=args.input_file
DECam_Root=file_name.split('.')[0]
ccd_num=DECam_Root.split('_')[2]

if(local_convert):
	#Move .jp2 to local directory
	print('Moving '+untar_path+file_name+' to '+local_dir+file_name)
	shutil.move(untar_path+file_name,local_dir+file_name)
	untar_path=local_dir

#Uncompress Fits on local Directory
uncompressed_fits=untar_path+DECam_Root+'.fits'
print('Uncompressing: '+file_name+' in path: '+untar_path)
subprocess.run(['j2f','-i',untar_path+file_name,'-o',uncompressed_fits,'-num_threads',str(1)])

#Extract nescessary information from file for naming scheme
exp=pyfits.getval(uncompressed_fits,"EXPNUM")
Field=pyfits.getval(uncompressed_fits,"OBJECT")
Filter=pyfits.getval(uncompressed_fits,"FILTER")[0]

#FOR Chile!
timestamp=datetime.datetime.utcnow().time()
if(timestamp > datetime.time(22,30)):
	ut='ut'+str(int(pyfits.getval(uncompressed_fits,"OBSID")[6:12])+1)
else:
	ut='ut'+pyfits.getval(uncompressed_fits,"OBSID")[6:12]
ut='ut160730'

obstype=pyfits.getval(uncompressed_fits,"OBSTYPE")

newname=Field+'.'+Filter+'.'+ut+'.'+str(exp)+'_'+ccd_num+'.fits'
if((obstype == 'dome flat') or (obstype == 'domeflat')):
	newname='domeflat.'+Filter+'.'+ut+'.'+str(exp)+'_'+ccd_num+'.fits'
if((obstype == 'zero') or (obstype == 'bias')):
	newname='bias.'+ut+'.'+str(exp)+'_'+ccd_num+'.fits'

ut_dir=photepipe_rawdir+ut+'/'
dest_dir=ut_dir+ccd_num+'/'


#Check to see if UT date & CCD Directories have been created
if not os.path.isdir(ut_dir): 
	print('Creating Directory: '+ut_dir)	
	os.makedirs(ut_dir)
else:
	print('Directory Exists: '+ut_dir)

if not os.path.isdir(dest_dir): 
	print('Creating Directory: '+dest_dir)	
	os.makedirs(dest_dir)
else:
	print('Directory Exists: '+dest_dir)	

#Move Uncompressed Fits File
print('Renaming '+file_name+' to '+newname)
print('And moving to: '+dest_dir)
shutil.move(uncompressed_fits,dest_dir+newname)

#Call Pipeloop for default CCD reduction
subprocess.run(['pipeloop.pl','-red',ut,ccd_num,'-id',str(exp)])

#Remove unescessary .jp2
print('Deleting: '+untar_path+file_name)
subprocess.run(['rm',untar_path+file_name])

