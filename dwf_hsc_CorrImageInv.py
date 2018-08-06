#!/usr/bin/env python3
import os, time
import math
import sys
import glob
import argparse
import warnings
import astropy.io.fits as pyfits
from collections import Counter

def get_fieldname(visit,im_dir):
	obs=glob.glob(im_dir+'CORR*'+visit+'*')
	if(len(obs) > 0):
		field_name=pyfits.getval(obs[0],"OBJECT")
		return(field_name)
	else:
		return("NoField")

def get_count(im_dir):
	files=glob.glob(im_dir+'CORR-*')
	visits=[f.split('-')[-2] for f in files]
	visit_counter=Counter(sorted(visits))
	return(visit_counter)

def print_count(visit_counter,im_path):
	print("Visit CCD Count in {0} :".format(im_path))
	print("----------------")
	for visit,count in sorted(visit_counter.items()):
		print(get_fieldname(visit,im_path),visit,count)	
	print('')

def main():
	HSC_DIR='/lustre/projects/p025_swin/pipes/HSC_Data/hsc/'
	mjd_id='02236'
	rerun='DWF_201802'
	filt='HSC-G'
	parser = argparse.ArgumentParser(description='Check the number of CORR images per visit for a given night', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--hsc_dir',metavar='DIRECTORY',type=str,default=HSC_DIR,help='HSC Root Directory')
	parser.add_argument('--rerun',metavar='RERUN',type=str,default=rerun,help='hscpipe rerun name')
	parser.add_argument('--mjd_id',metavar='MDID',type=str,default=mjd_id,help='hscpipe mjd id')
	parser.add_argument('--filter',metavar='FILTER',type=str,default=filt,help='filter name')
	parser.add_argument('--path',metavar='PATH',type=str,default='',help='filter name')
	parser.add_argument('--repeat',metavar='REPEAT',type=bool,default=False,help='Infinitely Repeat? True/False')

	args = parser.parse_args()
	if(args.path == ''):
		im_dir=args.hsc_dir+'rerun/'+args.rerun+'/'+args.mjd_id+'/'+args.filter+'/corr/'
	else:
		im_dir=args.path

	print_count(get_count(im_dir),im_dir)

	while(args.repeat):
		print_count(get_count(im_dir),im_dir)
		time.sleep(60)

if __name__ == '__main__':
	main()	
