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
#~/.astropy/config/astropy.cfg was getting messed up - seperate default (used by pipeloop?) and this
os.environ['XDG_CONFIG_HOME']='/home/fstars/.python3_config/'

import astropy.io.fits as pyfits

#Uncompress new file + create & submit assosciated qsub scripts
def dwf_prepipe_unpack(file_name,push_path,untar_path,qsub_path):
	ccdlist=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59']
	#Truncated CCD list for testing
	#ccdlist=['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'] #truncated version for test
	
	DECam_Root=file_name.split('.')[0]

	#Untar new file
	print('Unpacking:\t'+file_name)
	subprocess.run(['tar','-xf',push_path+file_name,'-C',untar_path])

	Exposure=DECam_Root.split('_')[1]

	#Create Qsub scripts for new file with n_per_ccd jobs per script
	n_per_ccd=15
	n_scripts=math.ceil(len(ccdlist)/n_per_ccd)
	print('Writing '+str(n_scripts)+' qsub scripts for '+file_name)
	for n in range(n_scripts):
		dwf_prepipe_qsubccds(DECam_Root,DECam_Root+'_q'+str(n+1),ccdlist[n_per_ccd*n:(n+1)*n_per_ccd],qsub_path,push_path)

#Write Qsub Script & submit to queue
def dwf_prepipe_qsubccds(filename_root,qroot,ccds,qsub_path,push_path):
	image_list=[filename_root+'_'+f+'.jp2' for f in ccds]

	qsub_out_dir=qsub_path+'out/'

	qsub_name=qsub_path+qroot+'.qsub'

	print('Creating Script: '+qsub_name+' for CCDs '+min(ccds)+' to '+max(ccds))
	qsub_file=open(qsub_name,'w')

	walltime='00:05:00'
	queue='sstar'
	nodes='1'
	ppn='16'

	qsub_file.write('#!/bin/csh \n')
	qsub_file.write('#PBS -N {0}\n'.format(qroot))
	qsub_file.write('#PBS -o {0}{1}.stdout\n'.format(qsub_path+'out/',qroot))
	qsub_file.write('#PBS -e {0}{1}.stderr\n'.format(qsub_path+'out/',qroot))
	qsub_file.write('#PBS -l walltime={0}\n'.format(walltime))
	qsub_file.write('#PBS -q {0}\n'.format(queue))
	qsub_file.write('#PBS -A p025_swin\n')
	qsub_file.write('#PBS -l nodes={0}:ppn={1}\n'.format(nodes,ppn))
	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo Automated script by dwf_prepipe\n')
	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write("echo -n 'Job is running on node '; cat $PBS_NODEFILE\n")
	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo PBS: qsub is running on $PBS_O_HOST\n')
	qsub_file.write('echo PBS: originating queue is $PBS_O_QUEUE\n')
	qsub_file.write('echo PBS: executing queue is $PBS_QUEUE\n')
	qsub_file.write('echo PBS: working directory is $PBS_O_WORKDIR\n')
	qsub_file.write('echo PBS: execution mode is $PBS_ENVIRONMENT\n')
	qsub_file.write('echo PBS: job identifier is $PBS_JOBID\n')
	qsub_file.write('echo PBS: job name is $PBS_JOBNAME\n')
	qsub_file.write('echo PBS: current home directory is $PBS_O_HOME\n')
	qsub_file.write('echo PBS: PATH = $PBS_O_PATH\n')
	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo PBS: PATH = $PBS_O_PATH\n')

	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo Just in case, do the below:\n')

	#This Shit keeps getting corrupted, write from know good ones 
	qsub_file.write('cp /home/fstars/.astropy/config/astropy.cfg.good /home/fstars/.astropy/config/astropy.cfg\n')
	qsub_file.write('cp /home/fstars/.python3_config/astropy/astropy.cfg.good /home/fstars/.python3_config/astropy/astropy.cfg\n')
	
	#Create the local directory if its not allready there
	#and delete everything inside since we're taking a full node
	qsub_file.write('mkdir /lfs/data0/dwf/\n')
	qsub_file.write('rm /lfs/data0/dwf/*.jp2\n')
	qsub_file.write('rm /lfs/data0/dwf/*.fits\n')
	qsub_file.write('echo ------------------------------------------------------\n')

	for f in image_list:	
		qsub_file.write('~/dwf_prepipe/dwf_prepipe_processccd.py -i {0} &\n'.format(f))
	qsub_file.write('wait\n')

	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo Safety Cleanup for the local disk:\n')

	qsub_file.write('rm /lfs/data0/dwf/*.jp2\n')
	qsub_file.write('rm /lfs/data0/dwf/*.fits\n')
	qsub_file.write('echo ------------------------------------------------------\n')

	qsub_file.close()
	subprocess.run(['qsub',qsub_name])

def main():
	DWF_Push = '/lustre/projects/p025_swin/fstars/DWF_Unpack_Test/push/' #"/lustre/projects/p025_swin/pipes/DWF_PIPE/CTIO_PUSH/"
	parser = argparse.ArgumentParser(description='Handle File Ingests for the DWF pipeline', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
	help='Directory where tarballs of compressed files are placed')
	#parser.add_argument('-p', dest='processes', type=int, default=multiprocessing.cpu_count(),
	#	help='Number of processes to plot with (default: #CPU cores)')
	args = parser.parse_args()

	path_to_watch = args.push_dir
	path_to_untar = args.push_dir+'untar/'
	path_to_qsub = args.push_dir+'qsub/'
	before = dict ([(f, None) for f in os.listdir (path_to_watch)])

	#dwf_prepipe_unpack('DECam_00504110.tar',path_to_watch,path_to_untar,path_to_qsub)
	print('Monitoring Directory:'+path_to_watch)
	while 1:
		after = dict ([(f, None) for f in os.listdir (path_to_watch)])
		added = [f for f in after if not f in before]
		removed = [f for f in before if not f in after]
		
		if added: print("Added: ", ", ".join (added))
		if removed: print("Removed: ", ", ".join (removed))
		
		for f in added: 
			dwf_prepipe_unpack(f,path_to_watch,path_to_untar,path_to_qsub)
		
		before = after
		time.sleep (5)


if __name__ == '__main__':
    main()	


