#!/usr/bin/env python3
#example usage ./dwf_prepipe.py /projects/p025_swin/fstars/DWF_Unpack_Test/push/
import astropy.io.fits as pyfits
import glob

def main():
	file_list=glob.glob('/lustre/projects/p025_swin/DWF_HSC_PushTest/fz/*.fz')
	file_list=sorted(file_list,reverse=True)
	files=file_list[:2000]
	for f in files:
		ccd=pyfits.getval(f,"DET-ID",ext=1)
		visit=pyfits.getval(f,"EXP-ID",ext=1)
		print(f,visit,ccd)

if __name__ == '__main__':
    main()	

