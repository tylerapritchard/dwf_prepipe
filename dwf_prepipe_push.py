#!/usr/bin/env python3
#-W error:"WARNING: File may have been truncated:*""
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/
import os, time
import math
import sys
import glob
import warnings
import argparse
import subprocess
import astropy.io.fits as pyfits

def dwf_prepipe_validatefits(file_name,data_dir):
	warnings.filterwarnings('error','.*File may have been truncated:.*',UserWarning)
	valid=0
	while(not valid):
		try:
			test=pyfits.open(data_dir+file_name+'.fits.fz')
		#except OSError:
		#	print('OS')
		#	print(file_name+' still writing ...')
		#	time.sleep(3)
		except UserWarning:
			print('user')
			print(file_name+' still writing ...')
			time.sleep(0.5)
		except IOError:
			print('io')
			print(file_name+' still writing ...')
			time.sleep(0.5)
		else:
			print('pass')
			valid=1

#Package new raw .fits.fz file 
def dwf_prepipe_packagefile(file_name,data_dir,Qs):
	file_name=file_name.split('/')[-1].split('.')[0]
	jp2_dir=data_dir+"jp2/"
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

#Parallel Ship file to G2
def dwf_prepipe_parallel_pushfile(file_name,data_dir):
	#g2 configuration
	user='fstars'
	host='g2.hpc.swin.edu.au'
	reciever=user+'@'+host
	push_dir='/lustre/projects/p025_swin/fstars/push/'
	target_dir='/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'

	jp2_dir=data_dir+"jp2/"

	print('Shipping:'+jp2_dir+file_name+'.tar')
	command="scp "+jp2_dir+file_name+".tar "+reciever+":"+push_dir+"; ssh "+reciever+" 'mv "+push_dir+file_name+".tar "+target_dir+"'"
	subprocess.Popen(command,shell=True)
	print('Returning to watch directory')

#Serial Ship to g2
def dwf_prepipe_serial_pushfile(file_name,data_dir):
	#g2 configuration
	user='fstars'
	host='g2.hpc.swin.edu.au'
	reciever=user+'@'+host
	push_dir='/lustre/projects/p025_swin/fstars/push/'
	target_dir='/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/'

	jp2_dir=data_dir+"jp2/"

	print('Shipping:'+jp2_dir+file_name+'.tar')
	command="scp "+jp2_dir+file_name+".tar "+reciever+":"+push_dir+"; ssh "+reciever+" 'mv "+push_dir+file_name+".tar "+target_dir+"'"
	subprocess.run(command,shell=True)
	print('Returning to watch directory')

#Input Keyword Default Values
DWF_PID = "/home4/images/fits/2016A-0095/"
Qs_Def=0.000038
method_def='p'
nbundle_def=4

#Parse Inputs
parser = argparse.ArgumentParser(description='DWF_Prepipe push script for raw data from CTIO', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-d','--data_dir',metavar='DIRECTORY',type=str,default=DWF_PID,
	help='Directory where tarballs of compressed files are placed')
parser.add_argument('-q','--Qs',metavar='DIRECTORY',type=float,default=Qs_Def,
	help='Qstep for fits2jpeg compression')
parser.add_argument('--method',metavar='PROTOCOL',type=str,default=method_def,
	help='File Transfer method:(s)erial, (p)arrallel, (b)undle, (l)ist')
parser.add_argument('--nbundle',metavar='NUMBER',type=str,default=nbundle_def,
	help='Number of Files to bundle together')

args = parser.parse_args()

path_to_watch=args.data_dir
Qs=args.Qs
method=args.method
nbundle=args.nbundle
#Begin Monitoring Directory
print('Monitoring:'+path_to_watch)
before = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])

while 1:
	after = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits.fz')])
	added = [f for f in after if not f in before]
	removed = [f for f in before if not f in after]
	if added: print("Added: ", ", ".join (added))
	if removed: print("Removed: ", ", ".join (removed))
	
	if ((method == 'p') and added):
		for f in added: 
			print('Processing: '+f)
			dwf_prepipe_validatefits(f,path_to_watch)
			dwf_prepipe_packagefile(f,path_to_watch,Qs)
			dwf_prepipe_parallel_pushfile(f,path_to_watch)

	if ((method == 's') and added):
		dwf_prepipe_validatefits(f,path_to_watch)
		print('Processing: '+added[-1])
		dwf_prepipe_packagefile(added[-1],path_to_watch,Qs)
		dwf_prepipe_serial_pushfile(added[-1],path_to_watch)
	
	if ((method == 'b') and added):
		if(len(added) > nbundle):
			added.sort()
			bundle=added[-1*nbundle:]
		else:
			bundle=added
		print(['Bundling:'+str(f) for f in bundle])
		for f in bundle: 
			print('Processing: '+f)
			dwf_prepipe_validatefits(f,path_to_watch)
			dwf_prepipe_packagefile(f,path_to_watch,Qs)
			#do all but the last scp in parallel; then force python to wait until the final transfer is complete
			if(f == bundle[len(bundle)-1]):
				dwf_prepipe_serial_pushfile(f,path_to_watch)
			else:
				dwf_prepipe_parallel_pushfile(f,path_to_watch)
	
	before = after
	time.sleep (1)
