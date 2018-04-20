import os,sys,glob,subprocess,time

if not os.path.exists('PDBs'):
	os.makedirs('PDBs')

pdbs = sorted(glob.glob('./Extracted/*/target.pdb'))

print(len(pdbs))

g = open("mapping.txt",'w')
g.write('Number\thash\n')

for x,pdb in enumerate(pdbs):
	hash = pdb.split('/')[2]
	num = x+1
	num_str = str(num)
	print(num_str)
	g.write(num_str+'\t'+hash+'\n')
	subprocess.call('cp '+pdb+' ./PDBs/'+num_str+'.pdb',shell=True)
g.close()

subprocess.call('python /switchlab/group/OptimizeProtein/SourceFiles/Scripts/pdb2fasta.py',shell=True)