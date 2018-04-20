import os,sys,glob,subprocess,time

start = os.getcwd()
print(start)

pdbs = sorted(glob.glob('PDBs/*.pdb'))
repairs = sorted(glob.glob('RepairPDBs/RepairPDB_*.pdb'))
os.chdir('RepairFiles')
print(len(pdbs))
pdbstobedone = []
repairsnkd = []
for pdb in repairs:
	name = pdb.split('/')[-1].split('_')[-1]
	repairsnkd.append(name)

for pdb in pdbs:
	name = pdb.split('/')[-1]
	if name not in repairsnkd:
		pdbstobedone.append(pdb)

pdb_chunks = [pdbstobedone[i:i + 2] for i in range(0, len(pdbstobedone), 2)]

for x,pdb_chunk in enumerate(pdb_chunks):
	jobname = 'job2_'+str(x)+'.q'
	g = open(jobname,'w')
	g.write('#!/bin/bash\n')
	g.write('#$ -N rp'+str(x)+'\n')
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
		g.write('/switchlab/group/tools/FoldX_2015/FoldX -manual '+start+'/'+pdb+' options_repair.txt commands_repair.txt\n')
		g.write('mv '+start+'/PDBs/RepairPDB_'+pdbname+' '+start+'/RepairPDBs/.\n')
		g.write('rm '+start+'/PDBs/RepairPDB_'+namenaked+'.fxout\n')
	g.close()
	subprocess.call('qsub '+jobname,shell=True)
