#!/bin/env python3

import os
import shutil
import argparse
import subprocess
import astropy.io.fits as pyfits


DWF_Push = '/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'

parser = argparse.ArgumentParser(description='DWF_Prepipe process for single compressed CCD', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-p','--push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
	help='Directory where tarballs of compressed files are placed')
parser.add_argument('-i','--input_file',type=str,help='input .jp2 file')
args = parser.parse_args()


photepipe_rawdir= '/projects/p025_swin/pipes/arest/DECAM/DEFAULT/rawdata/'

push_dir=args.push_dir
untar_path=push_dir+'untar/'

file_name=args.input_file
DECam_Root=file_name.split('.')[0]
ccd_num=DECam_Root.split('_')[2]
uncompressed_fits=untar_path+DECam_Root+'.fits'
print('Uncompressing: '+file_name+' in path: '+untar_path)
subprocess.run(['j2f','-i',untar_path+file_name,'-o',uncompressed_fits])

#Extract nescessary information from file for naming scheme
exp=pyfits.getval(uncompressed_fits,"EXPNUM")
Field=pyfits.getval(uncompressed_fits,"OBJECT")
Filter=pyfits.getval(uncompressed_fits,"FILTER")[0]
ut='ut'+pyfits.getval(uncompressed_fits,"OBSID")[6:12]

newname=Field+'.'+Filter+'.'+ut+'.'+str(exp)+'_'+ccd_num+'.fits'

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
subprocess.run(['pipeloop.pl','-red',ut,ccd_num,'-redobad'])#,'-im',newname])

#Remove unescessary .jp2
print('Deleting: '+untar_path+file_name)
subprocess.run(['rm',untar_path+file_name])

