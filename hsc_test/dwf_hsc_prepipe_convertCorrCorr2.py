#!/usr/bin/env python2

import os
import sys
import numpy
from math import *

from astropy.io import fits
import sip_tpv 
import pdb	
	
argvs = sys.argv # for run-time inputs

fname = argvs[0] # input image 
fname_out = argvs[1] # output image


fin = fits.open(fname)
fin_header = fin[1].header
sip_tpv.sip_to_pv(fin_header)
fin.writeto(fname_out,clobber="TRUE")
