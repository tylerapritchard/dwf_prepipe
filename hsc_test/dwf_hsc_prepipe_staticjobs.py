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
import astropy.io.fits as pyfits
import numpy as np
import re
from itertools import islice

def nth_index(file_list, exp, n):
	ind=0
	while(n > 0):
		if(file_list[ind].endswith(exp)):
			n=n-1
		ind=ind+1
	return(ind-1)

def dwf_hsc_getjobs():
	proc=subprocess.Popen(['qstat','-u','fstars'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,universal_newlines=True)
	stdout, stderr = proc.communicate()
	count=0
	for line in stdout.split(os.linesep):
		job=line.split()
		if(len(job) > 0):
			if(re.match('\d',job[0])):
				if(re.match('processccd*',job[3])):
					#if(job[9] == 'R'):
					count=count+1
	return(count)

def dwf_prepipe_getseeing(ccd,mjd_id):
	fwhm_path='/lustre/projects/p025_swin/pipes/DWF_PIPE/FWHM/'
	files=glob.glob(fwhm_path+'COSMOS_2236*')
	files=sorted(files)
	if(len(files) > 0):
		file=files[-1]
		mjd, ccds, fwhm = np.loadtxt(file,unpack=True)
		inds=np.where(ccds == int(ccd))#int(ccd))
		ind=inds[0]
	else:
		ind=[]
	try:
		default_path='/lustre/projects/p025_swin/pipes/HSC_Data/hsc_mary/dwf_prepipe_psf_default.dat'
		default_psf=np.loadtxt(default_path)
	except IOError:
		default_psf=1
		print('WARNING: No File for default PSF found, using default_psf={0}'.format(default_psf))
		print('WARNING: default PSF path: {0}'.format(default_path))		
	finally:
		if(len(ind) > 0):
			return(fwhm[ind[0]]*0.168)
		else:
			return(default_psf)

def dwf_prepipe_getvisitccd(file):
	ccd=pyfits.getval(file,"DET-ID")#,ext=1)
	visit=pyfits.getval(file,"EXP-ID")#,ext=1)
	visit=visit[5:]
	mjd=pyfits.getval(file,"MJD")#,ext=1)
	filt=pyfits.getval(file,"FILTER01")#,ext=1)
	filt=filt.upper()

	return(visit,ccd,mjd,filt)

#Write Qsub Script & submit to queue

def dwf_prepipe_qsubccds(visit,hsc_root,ccds,qsub_path,rerun,n):
	qroot='processccd-'+visit+'-q'+str(n)+'-'+rerun
	config_file='/projects/p025_swin/pipes/HSC_Data/hsc/minimaltest/hsc_fullmin.py'
	qsub_name=qsub_path+qroot+'.qsub'

	print('Writing Script: '+qsub_name)
	qsub_file=open(qsub_name,'w')

	walltime='00:30:00'
	queue='sstar'
	nodes='1'
	ppn='16'

	qsub_file.write('#!/bin/bash \n')
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

	qsub_file.write('echo $HOST\n')

	#Create the local directory if its not allready there
	#and delete everything inside since we're taking a full node
	qsub_file.write('source /projects/p025_swin/pipes/hscpipe/4.0.5/bashrc\n')
	qsub_file.write('setup-hscpipe\n')
	#qsub_file.write('setup -j astrometry_net_data sdss-dr9-fink-v5b\n')
	qsub_file.write('setup -j astrometry_net_data ps1_pv3_3pi_20170110-and\n')

	clobber=''
	
	default_path='/lustre/projects/p025_swin/pipes/HSC_Data/hsc_mary/dwf_prepipe_psf_default.dat'
	default_psf=np.loadtxt(default_path)

	j=0
	for ccd in ccds:
		psf_val=dwf_prepipe_getseeing(ccd,'2235')#to fix PSF to good measured value
		qsub_file.write('sleep '+str(0.5*j)+ ' ; '+'time hscProcessCcd.py {0} --rerun {1} -c isr.fwhm={4} -j 1 --id visit={2} ccd={3}  -C {5} --clobber-config --no-backup-config &\n'.format(hsc_root,rerun,visit,ccd,psf_val,config_file)) # Do the Full processing
		j=j+1

	qsub_file.write('wait\n')
	qsub_file.close()

	subprocess.run(['qsub',qsub_name])

def main():
	DWF_Push='/lustre/projects/p025_swin/DWF_HSC_PushTest/fz/'
	unzip_dir = '/projects/p025_swin/pipes/HSC_Data/push/'
	HSC_DIR='/lustre/projects/p025_swin/pipes/HSC_Data/hsc/'
	MARY_DIR='/lustre/projects/p025_swin/pipes/HSC_Mary/'
	rerun='DWF_201802'
	hsc_root='/projects/p025_swin/pipes/HSC_Data/hsc/'
	qsub_path='/lustre/projects/p025_swin/pipes/HSC_Data/qsub/'

	parser = argparse.ArgumentParser(description='Handle File Ingests for the DWF pipeline', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
	help='Directory where tarballs of compressed files are placed')

	args = parser.parse_args()

	path_to_watch = args.push_dir
	path_to_unzip = unzip_dir
	path_to_mary=MARY_DIR
	path_to_qsub = args.push_dir+'qsub/'

	ccdlist=[str(i) for i in range(104)]
	n_per_ccd=15
	n_scripts=math.ceil(len(ccdlist)/n_per_ccd)

	old_visit=''

	while 1:
		job_count=dwf_hsc_getjobs()
		print('Job Count: '+str(job_count))
		if(job_count < 15):
			file_list=glob.glob(unzip_dir+'*.fits')
			file_list=sorted(file_list,reverse=True)
			ind=nth_index(file_list,'00.fits',3)
			visit, ccd, mjd, filt = dwf_prepipe_getvisitccd(file_list[ind])
			if(visit != old_visit):
				print('Writing '+str(n_scripts)+' qsub scripts for '+visit)
				for n in range(n_scripts):
					print('Creating Script: '+visit+' for CCDs '+str(n_per_ccd*n)+' to '+str((n+1)*n_per_ccd))
					dwf_prepipe_qsubccds(visit,hsc_root,ccdlist[n_per_ccd*n:(n+1)*n_per_ccd],qsub_path,rerun,n)
				old_visit=visit
			else:
				print('Previous Visit=Current Visit, no new data?')

		time.sleep (20)


if __name__ == '__main__':
    main()	
