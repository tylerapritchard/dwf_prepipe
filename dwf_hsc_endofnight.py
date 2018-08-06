#!/usr/bin/env python3

import numpy as np
import time
import os
import astropy.io.fits as pyfits
import math
import glob
import re
import subprocess
def dwf_hsc_getjobs():
	proc=subprocess.Popen(['qstat','-u','fstars'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT,universal_newlines=True)
	stdout, stderr = proc.communicate()
	count=0
	for line in stdout.split(os.linesep):
		job=line.split()
		if(len(job) > 0):
			if(re.match('\d',job[0])):
				if(re.match('process*',job[3])):
					#if(job[9] == 'R'):
					count=count+1
	return(count)

def dwf_prepipe_getvisitccd(file):
	ccd=pyfits.getval(file,"DET-ID",ext=1)
	visit=pyfits.getval(file,"EXP-ID",ext=1)
	visit=visit[5:]
	mjd=pyfits.getval(file,"MJD",ext=1)
	filt=pyfits.getval(file,"FILTER01",ext=1)
	filt=filt.upper()
	ccd_num=str(ccd)
	if(len(ccd_num)==1): ccd_num='00'+str(ccd)
	if(len(ccd_num)==2): ccd_num='0'+str(ccd)

	return(visit,ccd_num,mjd,filt)

def dwf_prepipe_getseeing(ccd,mjd_id):
	fwhm_path='/lustre/projects/p025_swin/pipes/DWF_PIPE/FWHM/'
	files=glob.glob(fwhm_path+'COSMOS_2234*')
	files=sorted(files)
	if(len(files) > 0):
		file=files[-1]
		mjd, ccds, fwhm = np.loadtxt(file,unpack=True)
		inds=np.where(ccds == int(ccd))
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

def dwf_prepipe_qsubccds(files,path_to_watch,path_to_unzip,path_to_mary,path_to_qsub,HSC_DIR,rerun,eon):
	config_file='/projects/p025_swin/pipes/HSC_Data/hsc/minimaltest/hsc_fullmin.py'
	qsub_path='/lustre/projects/p025_swin/pipes/HSC_Data/qsub/'

	qsub_out_dir=qsub_path+'out/'

	qroot='process-eon-'+str(eon)+'-'+rerun

	qsub_name=qsub_path+qroot+'.qsub'
	print(qsub_name)
	print('Creating Script: '+qsub_name+' for '+str(len(files))+' files starting with:'+ files[0])
	print(files)

	qsub_file=open(qsub_name,'w')

	walltime='00:20:00'
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
	qsub_file.write('mkdir $PBS_JOBFS/dwf/\n')
	qsub_file.write('source /projects/p025_swin/pipes/hscpipe/4.0.5/bashrc\n')
	qsub_file.write('setup-hscpipe\n')
	#qsub_file.write('source /home/fstars/dwf_prepipe/hsc_test/envvar/hscpipe_env_var.sh\n')
	qsub_file.write('setup -j astrometry_net_data ps1_pv3_3pi_20170110-and\n')
	#qsub_file.write('setup -j sdss-dr9-fink-v5b\n')

	clobber=''
	psf_val='1.0'
	n=0
	for f in files:
		image_in=f
		visit, ccd, mjd, filt = dwf_prepipe_getvisitccd(image_in)
		ccd_num=str(ccd)
		if(len(ccd_num)==1): ccd_num='00'+str(ccd)
		if(len(ccd_num)==2): ccd_num='0'+str(ccd)
		mjd_id=str(mjd-55927)	

		psf_val=dwf_prepipe_getseeing(ccd,mjd_id)#to fix PSF to good measured value

		image=f.split('/')[-1][:-3]
		print(image)
		image_out=path_to_unzip+image

		com_script='('
		com_script=com_script+'sleep '+str(1*n)+ ' ; '

		image=path_to_unzip+f[:-3]
		com_script=com_script+'time funpack -O {1} {0}; sleep 1 ;'.format(image_in,image_out)

		com_script=com_script+'time hscIngestImages.py {0} {1}  -j 1 --mode=link --no-backup-config'.format(HSC_DIR,image_out)

		com_script=com_script+' ; time hscProcessCcd.py {0} --rerun {1} --id visit={2} ccd={3} -C {4} --no-backup-config --clobber-config'.format(HSC_DIR,rerun,visit,ccd,config_file)#--no-backup-config --clobber-config
		
		com_script=com_script+') &\n'
		qsub_file.write(com_script)
		n=n+1

	qsub_file.write('wait\n')

	qsub_file.write('echo ------------------------------------------------------\n')
	qsub_file.write('echo Safety Cleanup for the local disk:\n')

	qsub_file.write('rm $PBS_JOBFS/dwf/*.jp2\n')
	qsub_file.write('rm $PBS_JOBFS/dwf/*.fits\n')
	qsub_file.write('echo ------------------------------------------------------\n')

	qsub_file.close()
	subprocess.run(['qsub',qsub_name])

def main():

	fz_dir='/lustre/projects/p025_swin/DWF_HSC_PushTest/fz/'
	temp_dir='/lustre/projects/p025_swin/pipes/HSC_Data/test_temp'
	file_list=glob.glob(fz_dir+'*.fz')
	files=sorted(file_list,reverse=True)

	path_to_watch = fz_dir
	path_to_unzip = '/projects/p025_swin/pipes/HSC_Data/push/'
	path_to_mary='/lustre/projects/p025_swin/pipes/HSC_Mary/'
	path_to_qsub = fz_dir+'qsub/'
	HSC_DIR='/lustre/projects/p025_swin/pipes/HSC_Data/hsc/'
	rerun='DWF_201802'


	minimum=0#139922
	maximum=14183800 #14900000#14032525#14056004#14060600#14059310#14071436#14091600#14075600
	visit,ccd_num,mjd,filt = dwf_prepipe_getvisitccd(files[0])
	print(visit,minimum)
	n=0
	while((len(files) > 1) and (int(visit) > minimum)):
		index=0
		goodfiles=0

		batch=list()
		while(goodfiles < 15):
			fid=files[index][-16:-8]
			if(int(fid) < maximum):
				visit,ccd_num,mjd,filt = dwf_prepipe_getvisitccd(files[index])
				corr_name='CORR-'+visit+'-'+ccd_num+'.fits'
				if((not os.path.isfile('/projects/p025_swin/pipes/HSC_Data/hsc/rerun/DWF_201802/02235/HSC-G/corr/'+corr_name)) and (int(visit) < maximum)):
					print('Adding '+corr_name+' / '+files[index]+' to processing list')
					batch.append(files[index])
					goodfiles=goodfiles+1
			index=index+1
		while(dwf_hsc_getjobs() > 20):
			print("Current Jobs: "+ str(dwf_hsc_getjobs()))
			time.sleep(10)
		dwf_prepipe_qsubccds(batch,path_to_watch,path_to_unzip,path_to_mary,path_to_qsub,HSC_DIR,rerun,n)
		time.sleep(1)
		n=n+1
		files=files[index:]

if __name__ == '__main__':
    main()	
