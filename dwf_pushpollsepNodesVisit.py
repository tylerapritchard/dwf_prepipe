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
#import astropy.io.fits as pyfits
import datetime
import pyfits

def dwf_prepipe_validatefits(file_name):
	#print file_name
	warnings.filterwarnings('error','.*File may have been truncated:.*',UserWarning)
	valid=0
	propid = ''
	#while(not valid):
	if not valid:
		try:
			test=pyfits.open(file_name)
			propid = test[0].header['PROP-ID'].rstrip()
		except OSError:
			print('OS Error:')
			print(file_name+' still writing ...')
			time.sleep(3)
		except UserWarning:
			print('User Warning: ')
			print(file_name+' still writing ...')
			time.sleep(0.5)
		except IOError:
			print('IO Error:')
			print(file_name+' still writing ...')
			time.sleep(0.5)
		else:
			#print(file_name+' pass!')
			if propid == PROP_ID:
				valid=1
	return valid

def fzfiles(visit,inode = None):
	fitsSent = []
	for i in range(2):
		for ccd in validccds:
			if fzcomp == 'fz':
				fits0 = 'HSCA%06d%02d.fits.fz' %(visit+i,ccd)
			elif fzcomp == 'f2j':
				fits0 = 'HSCA%06d%02d.jp2' %(visit+i,ccd)
			else:
				fits0 = 'HSCA%06d%02d.fits' %(visit+i,ccd)
			fitsSent.append(fits0)
	nlen = len(fitsSent)

	if inode == None:
		inode = int(os.uname()[1].split('.')[0].replace('hsck-kn0',''))-iskipnode		

	if inode >= 0:
		nb = []
		for i in range(nnodes):
			nb.append(nlen/nnodes*i)
		nb.append(nlen)

		fitsSentUse = fitsSent[nb[inode]:nb[inode+1]]
	else:
		fitsSentUse = fitsSent

	return fitsSentUse

def fzfilesVisit(visit,inode = None):
	fitsSent = []
	for i in range(2):
		for ccd in validccds:
			if fzcomp == 'fz':
				fits0 = 'HSCA%06d%02d.fits.fz' %(visit+i,ccd)
			elif fzcomp == 'f2j':
				fits0 = 'HSCA%06d%02d.jp2' %(visit+i,ccd)
			else:
				fits0 = 'HSCA%06d%02d.fits' %(visit+i,ccd)
			fitsSent.append(fits0)
	nlen = len(fitsSent)

	if inode == None:
		inode = int(os.uname()[1].split('.')[0].replace('hsck-kn0',''))-iskipnode		

	if int(visit)%(nnodes*2)/2 == inode:
		fitsSentUse = fitsSent
	else:
		fitsSentUse = []

	#print visit,inode,fitsSentUse

	return fitsSentUse

#Package new raw .fits.fz file 
def dwf_prepipe_packagefile(file,data_dir,Qs):
	visit = file
	file_name = str(visit)

	fitsSentUse = fzfilesVisit(visit)

	inode = int(os.uname()[1].split('.')[0].replace('hsck-kn0',''))-iskipnode		

	if fzcomp == 'fz' or fzcomp == 'f2j':
		print('Compressing:'+file_name)
		procs = []
		for fits in fitsSentUse:
			if fzcomp == 'fz':
				command = 'fpack -S %s%s > fz/%s' %(data_dir,fits.replace('.fz',''),fits)
			elif fzcomp == 'f2j':
				command = ''
			#print command
			procs.append(subprocess.Popen(command,shell=True))
		#stop

		while True:
			finish = True
			for proc in procs:
				#print proc.poll()
				if proc.poll() != 0:
					finish = False
			if finish:
				break
			time.sleep (1)

	print('Packaging:'+jp2_dir+file_name+'.tar')
	if fzcomp == 'fz':
		command = 'tar cf %s%s-%d.tar -C %s %s; rm -f fz/{%s}' %(jp2_dir,file_name,inode,'fz/',' '.join(fitsSentUse),','.join(fitsSentUse))
	else:
		command = 'tar cf %s%s-%d.tar -C %s %s' %(jp2_dir,file_name,inode,data_dir,' '.join(fitsSentUse))
	procs = subprocess.Popen(command,shell=True)
	procs.wait()
	return 

#Parallel Ship file to G2
def dwf_prepipe_parallel_pushfileNodes(file):
	#file_name=file.split('/')[-1].split('.')[0]
	file_name=str(file)

	#g2 configuration

	jp2_dir="./"#data_dir+"/"

	inode = int(os.uname()[1].split('.')[0].replace('hsck-kn0',''))-iskipnode		

	print('Shipping(parallel):'+jp2_dir+file_name)
	fitsSentUse = fzfilesVisit(file)
	outdir = fz_dir
	funpack = ''
	if fzcomp == 'fz':
		for fits in fitsSentUse:
			funpack += 'funpack -O %s/%s %s/%s; ' %(target_dir,fits.replace('.fz',''),fz_dir,fits)
	else:
		outdir = target_dir
	funpack = ''
	command="scp %s%s-%d.tar %s:%s; ssh %s 'tar xf %s%s-%d.tar -C %s; %s rm %s%s-%d.tar'; rm %s%s-%d.tar" %(jp2_dir,file_name,inode,reciever,push_dir,reciever,push_dir,file_name,inode,outdir,funpack,push_dir,file_name,inode,jp2_dir,file_name,inode)
	print command
	procs = subprocess.Popen(command,shell=True)
	procs.wait()
	print('Returning to watch directory')
	return True

def dwf_prepipe_parallel_pushfile(file,data_dir):
	#file_name=file.split('/')[-1].split('.')[0]
	file_name=str(file)

	#g2 configuration

	jp2_dir="./"#data_dir+"/"

	print('Shipping(parallel):'+jp2_dir+file_name)
	procs = []
	for i in range(nnodes):
		fitsSentUse = fzfilesVisit(file,inode=i)
		funpack = ''
		for fits in fitsSentUse:
			funpack += 'funpack -O %s/%s %s/%s; ' %(target_dir,fits.replace('.fz',''),fz_dir,fits,) 
		command="scp %s%s-%d.tar %s:%s; ssh %s 'tar xf %s%s-%d.tar -C %s; %s rm %s%s-%d.tar'; rm %s%s-%d.tar" %(jp2_dir,file_name,i,reciever,push_dir,reciever,push_dir,file_name,i,fz_dir,funpack,push_dir,file_name,i,jp2_dir,file_name,i)
		print command
		procs.append(subprocess.Popen(command,shell=True))
	print('Returning to watch directory')
	return procs

#Serial Ship to g2
def dwf_prepipe_serial_pushfile(file,data_dir):
	#file_name=file.split('/')[-1].split('.')[0]
	file_name=str(file)+'.tar'
	#g2 configuration
	user='azenteno'
	host='g2.hpc.swin.edu.au'

	jp2_dir="./"#data_dir+"/"

	print('Shipping(serial):'+jp2_dir+file_name)
	command="scp "+jp2_dir+file_name+" "+reciever+":"+push_dir+"; ssh "+reciever+" 'tar xf "+push_dir+file_name+" -C "+target_dir+"; rm "+push_dir+file_name+"'; rm "+jp2_dir+file_name
	print command
	procs = subprocess.Popen(command,shell=True)
	procs.wait()
	print('Returning to watch directory')

def dwf_prepipe_cleantemp(file,data_dir):
	##Clean Temperary files - unpacked .fits, bundler .tar and individual .jp2
	file_name=file.split('/')[-1].split('.')[0]
	jp2_dir=data_dir+"jp2/"+file_name
	fits_name=file_name+'.fits'
	#remove funpacked .fits file 
	print('Removing: '+data_dir+fits_name)
	os.remove(data_dir+fits_name)
	#remove excess .tar
	#print('Removing: '+jp2_dir+file_name+'.tar')
	#os.remove(jp2_dir+'.tar')
	#Remove .jp2 files
	print('Cleaning: '+jp2_dir+'/')
	[os.remove(jp2_dir+'/'+jp2) for jp2 in os.listdir(jp2_dir) if jp2.endswith(".jp2")]

def dwf_prepipe_endofnight(data_dir,exp_min,Qs):

	#Get list of files in remote target directory & list of files in local directory
	remote_list=subprocess.getoutput("ssh dwf 'ls "+target_dir+"*.tar'")
	sent_files=[file.split('/')[-1].split('.')[0] for file in remote_list.splitlines() if file.endswith(".tar")]
	obs_list=[f.split('/')[-1].split('.')[0] for f in glob.glob(data_dir+'*.fits')]

	obs_list.sort(reverse=True)
	sent_files.sort(reverse=True)

	missing=[f for f in obs_list if not f in sent_files]

	print('Starting end of night transfers for general completion')
	print('Missing Files: '+str(len(missing))+'/'+str(len(obs_list))+' ('+str(len(sent_files))+' sent)')

	for f in missing:
		exp=int(f.split('_')[1])
		if(exp > exp_min):
			print('Processing: '+f)
			dwf_prepipe_packagefile(f,data_dir,Qs)
			dwf_prepipe_serial_pushfile(f,data_dir)
			dwf_prepipe_cleantemp(f,data_dir)		

def visitUpdate(path_to_watch,visits,sentvisits,taringvisits,taredvisits,pollingvisits):

	for visit in sentvisits:
		if visit in visits:
			print 'remove(s) ',visit,'from visits'
			visits.remove(visit)
	for visit in taringvisits:
		if visit in visits:
			print 'remove(t) ',visit,'from visits'
			visits.remove(visit)
	for visit in taredvisits:
		if visit in visits:
			print 'remove(t) ',visit,'from visits'
			visits.remove(visit)
	for visit in pollingvisits:
		if visit in visits:
			print 'remove(p) ',visit,'from visits'
			visits.remove(visit)

	fitsfiles = glob.glob(path_to_watch+'HSCA*00.fits')
	for fits1 in fitsfiles:
		visit =int(os.path.basename(fits1).replace('HSCA','').replace('.fits',''))/100
		if visit%2 == 0 and not visit in visits and not visit in sentvisits and not visit in pollingvisits\
			    and not visit in taringvisits and not visit in taredvisits \
			    and visit >= startvisit and visit <= endvisit:
			validAll = True
			fitsSentUse = fzfilesVisit(visit)
			if len(fitsSentUse) != 0:
				for fits in fitsSentUse:
					fits0 = fits.replace('.fz','')
					valid = dwf_prepipe_validatefits(path_to_watch+fits0)
					if not valid:
						validAll = False
			else:
				validAll = False
			if validAll:
				visits.append(visit)
				processCheck[visit] = dict()		
				processCheck[visit]['atHilo'] = datetime.datetime.now()


	print 'visits',visits #,processCheck
	print 'sentvisits',sentvisits
	print 'taringvisits',taringvisits
	print 'taredvisits',taredvisits
	print 'pollingvisits',pollingvisits

	return visits

def checkTarball(visit):
	while True:
		finish = True
		for inode in range(nnodes):
			dt = datetime.datetime.now() - datetime.datetime.fromtimestamp(os.stat('%s%d-%d.tar' %(jp2_dir,visit,inode)).st_mtime)
			if dt.seconds < 5.:
				finish = False
		if finish:
			break
		time.sleep(1)
	return finish

def main():
	#Input Keyword Default Values
	Qs_Def=0.000055
	method_def='b'
	exp_min=-1
	#Parse Inputs
	parser = argparse.ArgumentParser(description='DWF_Prepipe push script for raw data from CTIO', formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-d','--data_dir',metavar='DIRECTORY',type=str,default=DWF_PID,
		help='Directory where tarballs of compressed files are placed')
	parser.add_argument('-q','--Qs',metavar='DIRECTORY',type=float,default=Qs_Def,
		help='Qstep for fits2jpeg compression')
	parser.add_argument('--method',metavar='PROTOCOL',type=str,default=method_def,
		help='File Transfer method:(s)erial, (p)arrallel, (b)undle, (l)ist, (e)nd of night')
	parser.add_argument('--nbundle',metavar='NUMBER',type=int,default=nbundle_def,
		help='Number of Files to bundle together')
	parser.add_argument('--exp_min',metavar='NUMBER',type=int,default=exp_min,
		help='Exposure Number Start for end of night file transfer catchup')

	args = parser.parse_args()

	path_to_watch=args.data_dir
	Qs=args.Qs
	method=args.method
	nbundle=args.nbundle
	fout = open('send-%s.log' %(os.uname()[1].split('.')[0]),'a')
	timestart = datetime.datetime.now()
	print >> fout, '======== NodesVisit: start file transfer at %s ========' \
	    %(timestart.strftime('%Y/%m/%d %H:%M:%S'))
	print >> fout, 'Polling directory: %s' %(path_to_watch)
	print >> fout, 'PROP_ID: %s' %(PROP_ID)
	#Begin Monitoring Directory
	print('Monitoring:'+path_to_watch)
	#before = dict ([(f, None) for f in glob.glob(path_to_watch+'*.fits')])
	before = []
	visits = []
	sentvisits = []
	taringvisits = []
	taredvisits = []
	pollingvisits = []

	if(method == 'e'):
		dwf_prepipe_endofnight(path_to_watch, exp_min, Qs)
		return

	while 1:
		#sys.exit()

		visits = visitUpdate(path_to_watch,visits,sentvisits,taringvisits,taredvisits,pollingvisits)

		sortadd=visits
		sortadd.sort()		
			#print sortadd
		if(len(sortadd) > nbundle):
			bundle=sortadd[:nbundle]
		else:
			bundle=sortadd
		print(['Bundling:'+str(f) for f in bundle])
		if os.uname()[1].split('.')[0] != 'hsck-kn01':
			for visit in bundle[:nbundle]: 
				processCheck[visit]['tarstart'] = datetime.datetime.now()
				processCheck[visit]['status'] = 'taring'
				taringvisits.append(visit)
				dwf_prepipe_packagefile(visit,path_to_watch,Qs)
				processCheck[visit]['tarend'] = datetime.datetime.now()
				processCheck[visit]['status'] = 'tared'
				taringvisits.remove(visit)
				taredvisits.append(visit)
				print >> fout, '%d: (tar) %f sec' \
				    %(visit, (processCheck[visit]['tarend'] - processCheck[visit]['tarstart']).seconds)
				fout.flush()
				processCheck[visit]['scpstart'] = datetime.datetime.now()
				taredvisits.remove(visit)
				pollingvisits.append(visit)
				processCheck[visit]['status'] = 'polling'
				processCheck[visit]['scp'] = dwf_prepipe_parallel_pushfileNodes(visit)
				sentvisits.append(visit)
				pollingvisits.remove(visit)
				processCheck[visit]['scpend'] = datetime.datetime.now()
				processCheck[visit]['status'] = 'polled'
				print >> fout, '%d: (scp) %f sec' \
				    %(visit, (processCheck[visit]['scpend'] - processCheck[visit]['scpstart']).seconds)
				print >> fout, '%d: (at Hilo) %s' \
				    %(visit, processCheck[visit]['atHilo'].strftime('%Y/%m/%d %H:%M:%S'))
				print >> fout, '%d: (tarstart) %s' \
				    %(visit, processCheck[visit]['tarstart'].strftime('%Y/%m/%d %H:%M:%S'))
				print >> fout, '%d: (tarend) %s' \
				    %(visit, processCheck[visit]['tarend'].strftime('%Y/%m/%d %H:%M:%S'))
				print >> fout, '%d: (scpstart) %s' \
				    %(visit, processCheck[visit]['scpstart'].strftime('%Y/%m/%d %H:%M:%S'))
				print >> fout, '%d: (scpend) %s' \
				    %(visit, processCheck[visit]['scpend'].strftime('%Y/%m/%d %H:%M:%S'))
				fout.flush()
		time.sleep (1)
		if len(sentvisits) >= (endvisit-startvisit)/2+1 or len(taredvisits) >= (endvisit-startvisit)/2+1:
			print >> fout, 'Total: %f sec (%f min/54exp)' %((datetime.datetime.now()-timestart).seconds,(datetime.datetime.now()-timestart).seconds/60.) 
			break

if __name__ == '__main__':
	startvisit = 140988
	#startvisit = 137142
	endvisit = 137066
	endvisit = 999999
	## Don't miss the last "/"
	#DWF_PID = "/data/data2/rawData/HSC/poll_qbalance/refileDone.20180113HST/"
	DWF_PID = "/data/data2/rawData/HSC/poll_qbalance/" # for DWF run
	#DWF_PID = "/data/data2/rawData/HSC/poll_qbalance/refileDone/"
	user='azenteno'
	host='g2.hpc.swin.edu.au'
	reciever='dwf' #user+'@'+host
	push_dir='/lustre/projects/p025_swin/DWF_HSC_PushTest/'
	fz_dir='/lustre/projects/p025_swin/DWF_HSC_PushTest/fz/'
	target_dir='/lustre/projects/p025_swin/DWF_HSC_PushTest/push/'
	jp2_dir="./" #data_dir+"jp2/"

	PROP_ID = "o18139" # S18A-161 (Jeff's run)
	#PROP_ID = "o17142" # S17B-055I (Suzuki-san's run)

	fzcomp = 'fz' # None, fz, f2j

	if fzcomp == 'f2j':
		print 'f2j is not available yet.'
		exit

	nsep = 4

	nnodes = 6
	iskipnode = 9-nnodes

	nbundle_def=1
	nbundle_scp= 1

	validccds = range(58)
	validccds.remove(48)
	validccds.remove(49)
	validccds.remove(50)
	validccds.remove(51)
	validccds.remove(52)
	validccds.remove(57)
	
	processCheck = dict()

	main()	
