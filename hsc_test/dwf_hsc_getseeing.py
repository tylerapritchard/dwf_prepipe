#!/usr/bin/env python3

import os, time
import os, time
import math
import sys
import glob
import argparse
import warnings
import multiprocessing
import subprocess
import numpy as np
import astropy.io.fits as pyfits
#~/.astropy/config/astropy.cfg was getting messed up - seperate default (used by pipeloop?) and this
os.environ['XDG_CONFIG_HOME']='/home/fstars/.python3_config/'

def getSeeing(rerun,mjd_id,visit,ccd,default_fwhm):
	ccd_num=ccd
	if(len(ccd_num)==1): ccd_num='00'+ccd
	if(len(ccd_num)==2): ccd_num='0'+ccd
	qa_root='/projects/p025_swin/pipes/HSC_Data/hsc/rerun/{0}/{1}/HSC-G/qa/'.format(rerun,mjd_id)
	try:
		ccd,x,y,pix,ell,pa,a,b,e1,e2,ell_e1e2 = np.loadtxt(qa_root+'seeingGrid-'+visit+'-'+ccd_num+'.txt', dtype='f', unpack=True)
		return(np.nanmean(pix)*0.17)
	except FileNotFoundError:
		return(default_fwhm)

def main():
	rerun='DWF_201802'
	mjd_id='02234'
	visit='00000'
	default_fwhm='1.0'

	parser = argparse.ArgumentParser(description='Get Seeing Estimates from HSC Data', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--rerun', metavar='RERUN_ID',type=str,default=rerun)
	parser.add_argument('--mjd_id',metavar='MJD_ID',type=str,default=mjd_id)
	parser.add_argument('--visit',metavar='VISIT',type=str,default=visit)
	parser.add_argument('--default',metavar='DEFAULT_FWHM',type=str,default=default_fwhm)

	#parser.add_argument('-p', dest='processes', type=int, default=multiprocessing.cpu_count(),
	#	help='Number of processes to plot with (default: #CPU cores)')

	args = parser.parse_args()

	ccdlist=[str(i) for i in range(105)]	
	for ccd in ccdlist:
		print(visit,ccd,getSeeing(args.rerun,args.mjd_id,args.visit,ccd,args.default))

if __name__ == '__main__':
    main()	

