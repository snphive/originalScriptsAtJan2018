import os,sys,subprocess,time,glob

if os.path.exists('Runs'):
	subprocess.call('rm -rf Runs',shell=True)
