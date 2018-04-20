import os,sys,glob,subprocess,time

start = os.getcwd()
print(start)

pdbs = sorted(glob.glob('RepairPDBs/RepairPDB_*.pdb'))
os.chdir('SequenceDetailFiles')
print(len(pdbs))

pdb_chunks = [pdbs[i:i + 300] for i in range(0, len(pdbs), 300)]

for x,pdb_chunk in enumerate(pdb_chunks):
	jobname = 'job_'+str(x)+'.q'
	g = open(jobname,'w')
	g.write('#!/bin/bash\n')
	g.write('#$ -N sd'+str(x)+'\n')
	g.write('#$ -cwd\n')
	g.write('#$ -V\n')
	g.write('#$ -q all.q\n')
        g.write('#$ -l h_vmem=10G\n')
        if (x % 2 == 0):
                g.write('#$ -l hostname=hodor1.vib\n')
        else:
                g.write('#$ -l hostname=hodor2.vib\n')
	g.write('source ~/.bash_profile\n')
	for pdb in pdb_chunk:
		pdbname = pdb.split('/')[-1]
		namenaked = pdbname.split('.')[0]
		g.write('/switchlab/group/tools/FoldX_2015/FoldX -manual '+start+'/'+pdb+' options_seqdet.txt commands_seqdet.txt\n')
		g.write('mv '+start+'/RepairPDBs/SequenceDetail_'+namenaked+'.fxout '+start+'/SequenceDetails/.\n')
	g.close()
	subprocess.call('qsub '+jobname,shell=True)
