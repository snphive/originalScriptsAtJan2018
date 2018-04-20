# Create cysteine bridges
import os,sys,glob,subprocess,time,yasara

Mols = sys.argv[1].split('_')
R_Path = sys.argv[2]
FoldX_Path = sys.argv[3]
Agadir_Path = sys.argv[4]

yasara.info.mode = 'txt'

aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H','H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

startingDir = os.getcwd()
pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
name = pdb.split('/')[-1].split('.')[0]
rawpdb = pdb.split('/')[-1]


# SequenceDetail to get all residues
#seqdet = open('SequenceDetail_' + pdbfile[:-4] + '.fxout','r').read().splitlines()
seqdet = open('./Repair/SequenceDetail_'+name+'.fxout').readlines()
# put all residues in a list
allres = []
for line in seqdet:
	line = line.split('\t')
	if line[0] in aa_dict_3to1.keys():
		resname = line[0]
		resname1 = aa_dict_3to1[resname]
		chain = line[1]
		resnum = line[2]
		allres.append([resname,resname1,chain,resnum])

# Do PrintNetworks to find original cys bridges
cys_orig = []
cyslines = open('./Repair/AllAtoms_Disulfide_PrintNetworks_'+name+'.fxout','r').readlines()
for line in cyslines:
	if line[0:3] == 'CYS':
		line = line.split('\t')
		cys_orig.append(['C'+line[0][3:],'C'+line[2][3:]]) # [[CA23,CA88],[CA88,CA23],...]

yasara.run('DelObj all')
yasara.run('LoadPDB '+startingDir+'/Repair/'+name+'.pdb')
yasara.run('DelWater')
newCysBridges = []
for Mol in Mols:
	resnums = yasara.run('ListRes all and Mol '+Mol+',Format=RESNUM')
	residues = yasara.run('ListRes all and Mol '+Mol+',Format=RESNAME')
	for y,resnum in enumerate(resnums):
		contacts = yasara.run('ListConAtom CA and Res '+str(residues[y])+' '+str(resnum)+' and Mol '+Mol+',CA,Cutoff=7,Results=2')
		print contacts
		contacts = contacts[1::2]
		print contacts
		for x,contact in enumerate(contacts):
			print x
			print contact
			res1 = yasara.run('ListAtom '+str(contact)+',Format=RESNAME')[0]
			res1_1 = aa_dict_3to1[res1]
			num1 = yasara.run('ListAtom '+str(contact)+',Format=RESNUM')[0]
			mol1 = yasara.run('ListAtom '+str(contact)+',Format=MOLNAME')[0]
			print res1 + str(num1)
			newCysBridgeString = aa_dict_3to1[str(residues[y])]+Mol+str(resnum)+'C,'+res1_1+mol1+str(num1)+'C'
			newCysBridges.append(newCysBridgeString)

os.chdir('Runs/CysScan')
indiv = open('individual_list.txt','w')

combos = newCysBridges
name = name.split('_')[-1]
for mutation in combos:
	mutations = mutation.split(",")
	mutation = "_".join(mutations)
	if not os.path.exists(mutation):
		os.makedirs(mutation)
	os.chdir(mutation)
	if not os.path.exists('Fasta'):
		os.makedirs('Fasta')
	if not os.path.exists('Agadir'):
		os.makedirs('Agadir')
	if not os.path.exists('Agadir/Options.txt'):
		subprocess.call('cp '+startingDir+'/../../SourceFiles/AgadirFiles/* ./Agadir/.',shell=True)
	if os.path.exists(startingDir+'/Repair/RepairPDB_'+name+'.pdb'):
		subprocess.call('cp '+startingDir+'/Repair/RepairPDB_'+name+'.pdb .',shell=True)
		subprocess.call('cp '+startingDir+'/../../SourceFiles/FoldXFiles/* .',shell=True)
	else:
		print 'Something is wrong'
	targets = sorted(glob.glob('./PDBs/*.pdb'))
	for pdb in targets:
		name = pdb.split('/')[-1].split('.')[0]
		f = open(pdb).readlines()
		atomlines = []
		mols = []
		for line in f:
			if 'ATOM' == line[0:4]:
				mol = line[21]
				atomlines.append(line)
				if mol not in mols:
					mols.append(mol)
		for mol in mols:
			fastalist = []
			resnum = '0'
			for line in atomlines:
				if line[21] == mol and resnum!=line[22:26].strip(' '):
                        		resnum = line[22:26].strip(' ')
                        		aa = line[17:20]
					fastalist.append(aa_dict_3to1[aa])
			fasta = "".join(fastalist)
			#print name+'_'+mol
			#print fasta
			f = open(fastafolder+'/'+name+'_'+mol+'.fasta','w')
			f.write('>'+name+'_'+mol+'\n')
			f.write(fasta)
			f.close()
	f = open('runscript.txt','w')
	f.write('<TITLE>FOLDX_runscript;\n')
	f.write('<JOBSTART>#;\n')
	f.write('<PDBS>RepairPDB_'+name+'.pdb;\n')
	f.write('<BATCH>#;\n')
	f.write('<COMMANDS>FOLDX_commandfile;\n')
	f.write('<BuildModel>#,individual_list.txt;\n')
	f.write('<END>#;\n')
	f.write('<OPTIONS>FOLDX_optionfile;\n')
	f.write('<Temperature>298;\n')
	f.write('<IonStrength>0.05;\n')
	f.write('<ph>7;\n')
	f.write('<moveNeighbours>true;\n')
	f.write('<VdWDesign>2;\n')
	f.write('<numberOfRuns>3;\n')
	f.write('<OutPDB>#;\n')
	f.write('<END>#;\n')
	f.write('<JOBEND>#;\n')
	f.write('<ENDFILE>#;\n')
	f.close()
	h = open('individual_list.txt','w')
	if '_' in mutation:
		tempmuts = mutation.split('_')
		tempmuts = ",".join(tempmuts)
		h.write(tempmuts+';\n')
	else:
		h.write(mutation+';\n')
	h.close()
	g = open('./job.q','w')
	g.write('#!/bin/bash\n')
	g.write('#$ -N CS_'+mutation+'\n')
        g.write('#$ -V\n')
        g.write('#$ -cwd\n')
        g.write('source ~/.bash_profile\n')
	g.write(FoldX_Path+' -runfile runscript.txt\n')
	#g.write('python '+startingDir+'/../../SourceFiles/Scripts/agadir.py\n')
	g.close()
	if '_' in mutation:
		tempmuts = mutation.split('_')
		tempmuts = ",".join(tempmuts)
		indiv.write(mutation+';\n')
	else:
		indiv.write(mutation+';\n')
	subprocess.call('qsub job.q',shell=True)
	os.chdir('./..')

indiv.close()
check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'CS_' in output_qstat:
        print 'Waiting for all CysScan jobs to finish'
        time.sleep(10)
        check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()

dirs = sorted(glob.glob('./*'))
name_agad = name
print 'Name agad:\t'+name_agad
name = 'RepairPDB_'+name
print 'Name:\t'+name
for path in dirs:
	if os.path.isdir(path):
		os.chdir(path.split('/')[-1])
		mutation = path.split('/')[-1]
		subprocess.call('cp runscript.txt runscript_build.txt',shell=True)
		subprocess.call('rm runscript.txt',shell=True)
		f = open('runscript.txt','w')
		f.write('<TITLE>FOLDX_runscript;\n')
		f.write('<JOBSTART>#;\n')
		f.write('<PDBS>'+name+'_1_0.pdb,'+name+'_1_1.pdb,'+name+'_1_2.pdb,WT_'+name+'_1_0.pdb,WT_'+name+'_1_1.pdb,WT_'+name+'_1_2.pdb,;\n')
		f.write('<BATCH>#;\n')
		f.write('<COMMANDS>FOLDX_commandfile;\n')
		f.write('<AnalyseComplex>#;\n')
		f.write('<END>#;\n')
		f.write('<OPTIONS>FOLDX_optionfile;\n')
		f.write('<Temperature>298;\n')
		f.write('<IonStrength>0.05;\n')
		f.write('<ph>7;\n')
		f.write('<moveNeighbours>true;\n')
		f.write('<VdWDesign>2;\n')
		f.write('<numberOfRuns>3;\n')
		f.write('<OutPDB>#;\n')
		f.write('<END>#;\n')
		f.write('<JOBEND>#;\n')
		f.write('<ENDFILE>#;\n')
		f.close()
		pdb = name+'_1_0.pdb'
		print os.getcwd()
		f = open(pdb).readlines()
		atomlines = []
		mols = []
		for line in f:
			if 'ATOM' == line[0:4]:
				mol = line[21]
				atomlines.append(line)
				if mol not in mols:
					mols.append(mol)
		for mol in mols:
			fastalist = []
			resnum = '0'
			for line in atomlines:
				if line[21] == mol and resnum!=line[22:26].strip(' '):
					resnum = line[22:26].strip(' ')
					aa = line[17:20]
					fastalist.append(aa_dict_3to1[aa])
			fasta = "".join(fastalist)
			print name+'_'+mol
			print fasta
			f = open('Fasta/'+name+'_'+mol+'.fasta','w')
			f.write('>'+name+'_'+mol+'\n')
			f.write(fasta)
			f.close()
		g = open('./job.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N AC_'+mutation+'\n')
		g.write('#$ -V\n')
        	g.write('#$ -cwd\n')
        	g.write('source ~/.bash_profile\n')
		g.write(FoldX_Path+' -runfile runscript.txt\n')
		g.write('python '+startingDir+'/../../SourceFiles/Scripts/agadir.py\n')
		g.close()
		subprocess.call('qsub job.q',shell=True)
		os.chdir('./..')
	else:
		print path

check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'AC_' in output_qstat:
	print 'Waiting for all AnalyseComplex jobs to finish'
	time.sleep(10)
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()
os.chdir(startingDir)
g = open('SummaryCysScan.txt','w')
g.write('Mutation\tMol\tStretch\tddG\tdTANGO\tComplexSum\t')
analyseComplex = False
if os.path.isfile('./Repair/Interaction_AnalyseComplex_'+name+'.fxout'): 
	analyseComplex = True
	h = open('./Repair/Interaction_AnalyseComplex_'+name+'.fxout').readlines()
	Interactions = []
	for line in h[9:]:
		pieces = line.split('\t')
		mol1 = pieces[1]
		mol2 = pieces[2]
		molcomp = mol1+'_'+mol2
		Interactions.append(molcomp)
	interactionlist = []
	for i in Interactions:
		interactionlist.append('Complex_'+i)
	interactionstring = "\t".join(interactionlist)
	g.write(interactionstring+'\n')
else:
	g.write('\n')
dirs = sorted(glob.glob('./Runs/CysScan/*'))

path_agad_list = glob.glob('./Agadir/*')
TangoWT = 0
for i in path_agad_list:
	if os.path.isdir(i):
		print i
		f = open(i+'/PSX_globaltotal.out','r').readlines()
		pieces = f[1].split()
		TangoMol = float(pieces[2])
		TangoWT = TangoWT + TangoMol
print 'TangoTotalWT = '+str(TangoWT)


for path in dirs:
	if os.path.isdir(path):
		mut = path.split('/')[-1]
		mol = mut[1]
		f = open(path+'/Average_BuildModel_'+name+'.fxout','r').readlines()
		ddG = f[9].split()[2]
		path_agad_list = glob.glob(path+'/Agadir/*')
		TangoMut = 0
		for i in path_agad_list:
			if os.path.isdir(i):
				f = open(i+'/PSX_globaltotal.out','r').readlines()
				pieces = f[1].split()
				TangoMol = float(pieces[2])
				TangoMut = TangoMut + TangoMol
				path_agad = i
		print 'TangoTotalMut = '+str(TangoMut)
		print path_agad
		stretch = path_agad.split('/')[-1].split('_')[-2]
		f = open(path_agad+'/PSX_globaltotal.out','r').readlines()
		print f[1]
		mol = mut[1]
		if analyseComplex:
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_0.fxout','r').readlines()
			complex_Mut1 = []
			for line in f[9:]:
					complex_Mut1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_1.fxout','r').readlines()
			complex_Mut2 = []
			for line in f[9:]:
					complex_Mut2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_2.fxout','r').readlines()
			complex_Mut3 = []
			for line in f[9:]:
					complex_Mut3.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_0.fxout','r').readlines()
			complex_WT1 = []
			for line in f[9:]:
					complex_WT1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_1.fxout','r').readlines()
			complex_WT2 = []
			for line in f[9:]:
					complex_WT2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_2.fxout','r').readlines()
			complex_WT3 = []
			for line in f[9:]:
					complex_WT3.append(float(line.split()[5]))
			ComplexList=[]
			ComplexSum = 0
			for x,line in enumerate(f[9:]):
				complex_WT = float((complex_WT1[x]+complex_WT2[x]+complex_WT3[x])/3)
				complex_Mut = float((complex_Mut1[x]+complex_Mut2[x]+complex_Mut3[x])/3)
				Complex = complex_Mut-complex_WT
				ComplexSum = ComplexSum + Complex
				ComplexList.append(str(Complex))
			Complex = "\t".join(ComplexList)
		else:
			Complex = ""
			ComplexSum = 0
		dTango = float(TangoMut)-float(TangoWT)
		sumstring = mut+'\t'+mol+'\t'+stretch+'\t'+ddG+'\t'+str(dTango)+'\t'+str(ComplexSum)+'\t'+Complex+'\n'
		g.write(sumstring)
g.close()



