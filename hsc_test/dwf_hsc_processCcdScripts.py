#!/usr/bin/env python3
import argparse
import math
import subprocess
import numpy as np

def dwf_prepipe_getseeing(visit,ccd):
	ccd_num=ccd
	if(len(ccd_num)==1): ccd_num='00'+ccd
	if(len(ccd_num)==2): ccd_num='0'+ccd

	if (visit=='0001888'):
		qa_root='/projects/p025_swin/pipes/HSC_Data/hsc/rerun/cosmos_full/00823/HSC-G/qa/'
	else:
		qa_root='/projects/p025_swin/pipes/HSC_Data/hsc/rerun/cosmos_full/01052/HSC-G/qa/'

	try:
		ccd,x,y,pix,ell,pa,a,b,e1,e2,ell_e1e2 = np.loadtxt(qa_root+'seeingGrid-'+visit+'-'+ccd_num+'.txt', dtype='f', unpack=True)
	except FileNotFoundError:
		return(1.0)
	#print(np.nanmean(pix),np.nanmean(pix)*0.17)
	return(np.nanmean(pix)*0.17)

def dwf_prepipe_qsubccds(visit,hsc_root,ccds,qsub_path,rerun,n):
	qroot='hsc-processccd-'+visit+'-q'+str(n)+'-'+rerun
	config_file='/projects/p025_swin/pipes/HSC_Data/hsc/minimaltest/hsc_minimal_processCcd_t5.py'
	qsub_name=qsub_path+qroot+'.qsub'

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
	qsub_file.write('setup -j astrometry_net_data ps1_pv3_3pi_20170110-and\n')

	clobber=''
	psf_val='1.0'
	if(visit=='0001888'):
		psf_val='1.445'
		psf_val='1.25'
	if(visit=='0011690'):
		psf_val='0.748'
		psf_val='0.6'
	if(visit=='0011706'):
		psf_val='0.9'
		psf_val='1.1'
	for ccd in ccds:
		psf_val=dwf_prepipe_getseeing(visit,ccd)#to fix PSF to good measured value
		#print(psf_val) 
		#qsub_file.write('time hscProcessCcd.py {0} --rerun {1} --id visit={2} ccd={3}  --clobber-config &\n'.format(hsc_root,rerun,visit,ccd)) # Do the Full processing
		#qsub_file.write('time hscProcessCcd.py {0} --rerun {1} -C {4} --id visit={2} ccd={3}  &\n'.format(hsc_root,rerun,visit,ccd,config_file)) # --clobber-config --no-backup-config
		#qsub_file.write('time hscProcessCcd.py {0} --rerun {1} -C {4} -c isr.fwhm={5} --id visit={2} ccd={3}  --clobber-config &\n'.format(hsc_root,rerun,visit,ccd,config_file,psf_val)) # FixPSF
		#qsub_file.write('time hscProcessCcd.py {0} --rerun {1} -C {4} -c isr.fwhm={5} -c isr.doBrighterFatter=False --id visit={2} ccd={3}  --clobber-config &\n'.format(hsc_root,rerun,visit,ccd,config_file,psf_val)) # FixPSF, BrighterFatter Off
		qsub_file.write('time hscProcessCcd.py {0} --rerun {1} -C {4} -c isr.fwhm={5} --id visit={2} ccd={3}  --clobber-config &\n'.format(hsc_root,rerun,visit,ccd,config_file,psf_val)) # FixPSF, BrighterFatter ON

	qsub_file.write('wait\n')
	qsub_file.close()
	subprocess.run(['qsub',qsub_name])

def main():
	hsc_root='/projects/p025_swin/pipes/HSC_Data/hsc/'
	rerun='cosmos_min5_fixpsfccd'
	qsub_path='/projects/p025_swin/pipes/HSC_Data/hsc/queuetest/queuescript/'

	parser = argparse.ArgumentParser(description='Handle File Ingests for the DWF pipeline', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--visit',metavar='Visit Number',type=str, help='Visit Number for Files to be processed')
	#parser.add_argument('-p', dest='processes', type=int, default=multiprocessing.cpu_count(),
	#	help='Number of processes to plot with (default: #CPU cores)')
	args = parser.parse_args()
	visit=args.visit

	ccdlist=[str(i) for i in range(104)]
	n_per_ccd=15
	n_scripts=math.ceil(len(ccdlist)/n_per_ccd)
	
	print('Writing '+str(n_scripts)+' qsub scripts for '+visit)
	for n in range(n_scripts):
		print('Creating Script: '+visit+' for CCDs '+str(n_per_ccd*n)+' to '+str((n+1)*n_per_ccd))
		dwf_prepipe_qsubccds(visit,hsc_root,ccdlist[n_per_ccd*n:(n+1)*n_per_ccd],qsub_path,rerun,n)

if __name__ == '__main__':
    main()

