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

os.environ['XDG_CONFIG_HOME']='/home/fstars/.python3_config/'

import astropy.io.fits as pyfits
def dwf_prepipe_getseeing(ccd,mjd_id):
	fwhm_path='/lustre/projects/p025_swin/pipes/DWF_PIPE/FWHM/'
	files=glob.glob(fwhm_path+'COSMOS_2234*')
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
	ccd=pyfits.getval(file,"DET-ID",ext=1)
	visit=pyfits.getval(file,"EXP-ID",ext=1)
	visit=visit[5:]
	mjd=pyfits.getval(file,"MJD",ext=1)
	filt=pyfits.getval(file,"FILTER01",ext=1)
	filt=filt.upper()

	return(visit,ccd,mjd,filt)

#Write Qsub Script & submit to queue
def dwf_prepipe_qsubccds(files,path_to_watch,path_to_unzip,path_to_mary,path_to_qsub,HSC_DIR, rerun):
	config_file='/projects/p025_swin/pipes/HSC_Data/hsc/minimaltest/hsc_minimal_processCcd_t5.py'
	qsub_path='/lustre/projects/p025_swin/pipes/HSC_Data/qsub/'
	qsub_out_dir=qsub_path+'out/'

	qroot='hsc-processccd-'+files[0][:-3]+'-'+rerun

	qsub_name=qsub_path+qroot+'.qsub'

	print('Creating Script: '+qsub_name+' for '+str(len(files))+' files starting with:'+ files[0])
	print(files)

	qsub_file=open(qsub_name,'w')

	walltime='00:40:00'
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
	#qsub_file.write('setup -j astrometry_net_data ps1_pv3_3pi_20170110-and\n')
	qsub_file.write('setup -j sdss-dr9-fink-v5b\n')

	clobber=''
	psf_val='1.0'
	n=0
	for f in files:
		#print(f)
		im=path_to_watch+f
		visit, ccd, mjd, filt = dwf_prepipe_getvisitccd(im)
		ccd_num=str(ccd)
		if(len(ccd_num)==1): ccd_num='00'+str(ccd)
		if(len(ccd_num)==2): ccd_num='0'+str(ccd)

		mjd_id=str(mjd-55927)
		

		psf_val=dwf_prepipe_getseeing(ccd,mjd_id)#to fix PSF to good measured value


		#print([visit,ccd])
		com_script='('
		com_script=com_script+'sleep '+str(1*n)+ ' ; '
		#if(image[-3:] == '.fz'): 
		#path_to_unzip='$PBS_JOBFS/dwf/'
		image=path_to_unzip+f[:-3]
		com_script=com_script+'time funpack -O {1} {0}; sleep 5 ;'.format(im,image)
		#if(image[-4:] == '.jp2'): 
		#	com_script=com_script+'j2f {0} ; '.format(image)
		#	image=image[:-4]

		com_script=com_script+'time hscIngestImages.py {0} {1}  -j 1 --mode=link --no-backup-config'.format(HSC_DIR,image)
		
		#com_script=com_script+'; time hscProcessCcd.py {0} --rerun {1} -C {4} -c isr.fwhm={5} -c isr.doBrighterFatter=False --id visit={2} ccd={3} --no-backup-config --clobber-config'.format(HSC_DIR,rerun,visit,ccd,config_file,psf_val)#  --clobber-config
		#com_script=com_script+'; time hscProcessCcd.py {0} --rerun {1} --id visit={2} ccd={3}  --clobber-config'.format(HSC_DIR,rerun,visit,ccd)
		com_script=com_script+' ; time hscProcessCcd.py {0} --rerun {1} --id visit={2} ccd={3} --no-backup-config --clobber-config'.format(HSC_DIR,rerun,visit,ccd)#--no-backup-config --clobber-config
		
		#corr_dir=HSC_DIR+'rerun/'+rerun+'/'+mjd_id+'/'+filt+'/corr/'
		#corr_name='CORR-'+visit+'-'+ccd_num+'.fits'
		#com_script=com_script+' ./home/fstars/dwf_prepipe/hsc_test/dwf_hsc_prepipe_convertCorrCorr2.py {0} {1} ) &'.format(corr_dir+corr_name, path_to_mary+mjd_id+'-'+corr_name)
		
		#print(com_script)
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
	DWF_Push='/lustre/projects/p025_swin/DWF_HSC_PushTest/fz/'
	unzip_dir = '/projects/p025_swin/pipes/HSC_Data/push/'
	HSC_DIR='/lustre/projects/p025_swin/pipes/HSC_Data/hsc/'
	MARY_DIR='/lustre/projects/p025_swin/pipes/HSC_Mary/'
	rerun='DWF_201802'
	
	parser = argparse.ArgumentParser(description='Handle File Ingests for the DWF pipeline', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--push_dir',metavar='DIRECTORY',type=str,default=DWF_Push,
	help='Directory where tarballs of compressed files are placed')
	#parser.add_argument('-p', dest='processes', type=int, default=multiprocessing.cpu_count(),
	#	help='Number of processes to plot with (default: #CPU cores)')
	args = parser.parse_args()

	path_to_watch = args.push_dir
	path_to_unzip = unzip_dir#args.push_dir+'untar/'
	path_to_mary=MARY_DIR
	path_to_qsub = args.push_dir+'qsub/'



	before = dict ([(f, None) for f in os.listdir (path_to_watch)])
	#before = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fz')])
	#dwf_prepipe_unpack('DECam_00504110.tar',path_to_watch,path_to_untar,path_to_qsub)
	print('Monitoring Directory:'+path_to_watch)
	print('Uncompressing To: '+path_to_unzip)
	print('Writing Scripts to: '+path_to_qsub)
	print('Rerun Name: '+rerun)

	remaining=[]
	while 1:
		included_extenstions = ['fz', 'jp2', 'fits']
		file_names = [fn for fn in os.listdir(path_to_watch) if any(fn.endswith(ext) for ext in included_extenstions)]
		#file_names=glob.glob(path_to_watch+'*.fz')

		after = dict ([(f, None) for f in file_names])

		#after = dict ([(f, None) for f in os.listdir (path_to_watch)])
		added = [f for f in after if not f in before]
		removed = [f for f in before if not f in after]
		
		if added: print("Added: ", ", ".join (added))
		if removed: print("Removed: ", ", ".join (removed))
		#if remaining: print("Remaining: ",",".join(remaining))

		if(remaining != []): 
			files=remaining+added
		else:
			files=added

		while(len(files) > 13):
			#print(files)
			time.sleep(2)
			dwf_prepipe_qsubccds(files[:15],path_to_watch,path_to_unzip,path_to_mary,path_to_qsub,HSC_DIR,rerun)
			files=files[16:]
			time.sleep(2)

		remaining=files

		before = after
		time.sleep (5)

if __name__ == '__main__':
    main()	

