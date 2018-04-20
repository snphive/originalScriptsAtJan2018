import os,sys,glob,subprocess,time

start = os.getcwd()
print(start)

pdbs = sorted(glob.glob('PDBs/*.pdb'))
os.chdir('RepairFiles')
print(len(pdbs))

pdb_chunks = [pdbs[i:i + 20] for i in range(0, len(pdbs), 20)]

for x,pdb_chunk in enumerate(pdb_chunks):
	jobname = 'job_'+str(x)+'.q'
	g = open(jobname,'w')
	g.write('#!/bin/bash\n')
	g.write('#$ -N j'+str(x)+'\n')
	g.write('#$ -cwd\n')
	g.write('#$ -V\n')
	g.write('#$ -q all.q\n')
	g.write('source ~/.bash_profile\n')
	for pdb in pdb_chunk:
		pdbname = pdb.split('/')[-1]
		namenaked = pdbname.split('.')[0]
		g.write('/switchlab/group/tools/FoldX_2015/FoldX -manual '+start+'/'+pdb+' options_repair.txt commands_repair.txt\n')
		g.write('mv '+start+'/PDBs/RepairPDB_'+pdbname+' '+start+'/RepairPDBs/.\n')
		g.write('rm '+start+'/PDBs/RepairPDB_'+namenaked+'.fxout\n')
	g.close()
	subprocess.call('qsub '+jobname,shell=True)
