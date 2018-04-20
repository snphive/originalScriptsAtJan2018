# Supercharging
# Search in Molecule given for mutations that increase the netcharge
# Go negative and positive
# After scan of repaired pdb, get best netcharger that don't destabilize the complex, either adding a charge or removing a residue of the opposite charge.
# Continue until there are no residues left to mutate that can increase the netcharge.
import os,sys,glob,subprocess,time,yasara
startingDir = os.getcwd()

Mols = sys.argv[1].split('_')
R_Path = sys.argv[2]
FoldX_Path = sys.argv[3]
Agadir_Path = sys.argv[4]
Charge = sys.argv[5]

# amino acid dictionaries
yasara.info.mode = 'txt'

aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H','H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

pos = ['R','K']
neg = ['E','D']
pdb_start = glob.glob(startingDir+'/Repair/RepairPDB*.pdb')[0]
name = pdb_start.split('/')[-1].split('.')[0].split('_')[-1]
print name
print pdb_start

mol = Mols[0]

oldfasta_lines = open(startingDir+'/Fasta/'+name+'_'+mol+'.fasta','r').readlines()
oldfasta = []
if len(oldfasta_lines) == 2:
	for a in oldfasta_lines[1]:
		oldfasta.append(a)
else:
	for oldfasta_line in oldfasta_lines[1:]:
		for a in oldfasta_line:
			oldfasta.append(a)
oldfasta = "".join(oldfasta)

netcharge = 0
for a in oldfasta:
	if a in pos:
		netcharge = netcharge +1
		print "positive residue: "+a
	if a in neg:
		netcharge = netcharge -1
		print "negative residue: "+a
print "netcharge = "+str(netcharge)
print Charge
if Charge == 'S':
	if netcharge > 0:
		charges = pos
		anticharges = neg
	if netcharge < 0:
		charges = neg
		anticharges = pos
	if netcharge == 0:
		print 'No net charge in molecule'
if Charge == 'P':
	charges = pos
	anticharges = neg
if Charge == 'N':
	charges = neg
	anticharges = pos
round = 1
h = open('./Repair/Interaction_AnalyseComplex_RepairPDB_'+name+'.fxout').readlines()
Interactions = []
for line in h[9:]:
		pieces = line.split('\t')
		mol1 = pieces[1]
		mol2 = pieces[2]
		molcomp = mol1+'_'+mol2
		Interactions.append(molcomp)

summary = open('SummarySupercharge.txt','w')
header = 'Round\tMutation\tNetcharge\tChargeDif\tddG\tComplexSum\t'
for i in Interactions:
	header = header +'Complex_'+i+'\t'
header = header + '\n'
summary.write(header)
superchargefolder = startingDir+'/Runs/Supercharge'
os.chdir(superchargefolder)

while round < 2:
	RoundString = 'Round_'+str(round)
	print RoundString
	PrevRoundString = 'Round_'+str(round-1)
	if round == 1:
		seqdet_lines = open(startingDir+'/Repair/SequenceDetail_RepairPDB_'+name+'.fxout').readlines()
		pdb_path = pdb_start
		print pdb_path
	else:
		seqdet_lines = open(startingDir+'/Repair/SequenceDetail_'+PrevRoundString+'_supercharged.fxout').readlines()
		pdb_path = startingDir+'/Repair/'+PrevRoundString+'_supercharged.pdb'
	if not os.path.exists(RoundString):
		os.makedirs(RoundString)
	os.chdir(RoundString)
	mutations = []
	for seqdet_line in seqdet_lines[9:]:
		pieces = seqdet_line.split('\t')
		if len(pieces) == 1:
			continue
		if len(pieces) == 0:
			continue
		if pieces[1] == 'GAP':
			continue
		if pieces[1] not in aa_dict_3to1.keys():
			continue
		res3 = pieces[1]
		res1 = aa_dict_3to1[res3]
		resnum = pieces[3]
		mol = pieces[2]
		if res1 in charges:
			continue
		if not mol == Mols[0]:
			continue
		if res1 in anticharges:
			for mut in AAlist:
				if res1 == mut:
					continue
				else:
					mutation = [res1,mol,resnum,mut]
					mutation = "".join(mutation)
					mutations.append(mutation)
		else:
			for mut in charges:
				mutation = [res1,mol,resnum,mut]
				mutation = "".join(mutation)
				mutations.append(mutation)
	indiv_round = open(superchargefolder+'/'+RoundString+'_individual_list.txt','w')
	for mutation in mutations:
		indiv_round.write(mutation+';\n')
		if not os.path.exists(mutation):
			os.makedirs(mutation)
		os.chdir(mutation)
		if os.path.exists(pdb_path):
			subprocess.call('cp '+pdb_path+' .',shell=True)
			subprocess.call('cp '+startingDir+'/../../SourceFiles/FoldXFiles/* .',shell=True)
			name_pdb = pdb_path.split('/')[-1]
		else:
			print 'Something is wrong'
		f = open('runscript.txt','w')
		f.write('<TITLE>FOLDX_runscript;\n')
		f.write('<JOBSTART>#;\n')
		f.write('<PDBS>'+name_pdb+';\n')
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
		h.write(mutation+';\n')
		h.close()
		g = open('./job.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N SC_'+mutation+'\n')
		g.write('#$ -V\n')
        	g.write('#$ -cwd\n')
        	g.write('source ~/.bash_profile\n')
		g.write(FoldX_Path+' -runfile runscript.txt\n')
		g.close()
		subprocess.call('qsub job.q',shell=True)
		os.chdir('./..')
	indiv_round.close()
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()
	while 'SC_' in output_qstat:
		print 'Waiting for all '+RoundString+' Supercharge jobs to finish'
		time.sleep(10)
		check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
		output_qstat = check_qstat.stdout.read()

	dirs = sorted(glob.glob('./*'))
	for path in dirs:
		if os.path.isdir(path):
			os.chdir(path.split('/')[-1])
			subprocess.call('cp runscript.txt runscript_build.txt',shell=True)
			subprocess.call('rm runscript.txt',shell=True)
			f = open('runscript.txt','w')
			f.write('<TITLE>FOLDX_runscript;\n')
			f.write('<JOBSTART>#;\n')
			f.write('<PDBS>'+name_pdb.split('.')[0]+'_1_0.pdb,'+name_pdb.split('.')[0]+'_1_1.pdb,'+name_pdb.split('.')[0]+'_1_2.pdb,WT_'+name_pdb.split('.')[0]+'_1_0.pdb,WT_'+name_pdb.split('.')[0]+'_1_1.pdb,WT_'+name_pdb.split('.')[0]+'_1_2.pdb,;\n')
			f.write('<BATCH>#;\n')
			f.write('<COMMANDS>FOLDX_commandfile;\n')
			f.write('<AnalyseComplex>#;\n')
			f.write('<SequenceDetail>#;\n')
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
			subprocess.call('qsub job.q',shell=True)
			os.chdir('./..')
		else:
			print path
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()
	while 'SC_' in output_qstat:
			print 'Waiting for all '+RoundString+' AnalyseComplex jobs to finish'
			time.sleep(10)
			check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
			output_qstat = check_qstat.stdout.read()
	g = open(superchargefolder+'/'+RoundString+'_Summary.txt','w')
	g.write(header)
	dirs = sorted(glob.glob('./*'))
	ddG_min = float(10)
        best_mut = ""
        best_path = ""
        best_summaryString = ""
        best_netcharge = 0
	for path in dirs:
		if os.path.isdir(path):
			mut = path.split('/')[-1]
			print RoundString+'\n'+path
			f = open(path+'/Average_BuildModel_'+name_pdb.split('.')[0]+'.fxout','r').readlines()
			ddG = f[9].split()[2]
			f = open(path+'/Interaction_AnalyseComplex_'+name_pdb.split('.')[0]+'_1_0.fxout','r').readlines()
			complex_Mut1 = []
			for line in f[9:]:
				complex_Mut1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name_pdb.split('.')[0]+'_1_1.fxout','r').readlines()
			complex_Mut2 = []
			for line in f[9:]:
				complex_Mut2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name_pdb.split('.')[0]+'_1_2.fxout','r').readlines()
			complex_Mut3 = []
			for line in f[9:]:
				complex_Mut3.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name_pdb.split('.')[0]+'_1_0.fxout','r').readlines()
			complex_WT1 = []
			for line in f[9:]:
				complex_WT1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name_pdb.split('.')[0]+'_1_1.fxout','r').readlines()
			complex_WT2 = []
			for line in f[9:]:
				complex_WT2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name_pdb.split('.')[0]+'_1_2.fxout','r').readlines()
			complex_WT3 = []
			for line in f[9:]:
				complex_WT3.append(float(line.split()[5]))
			ComplexList=[]
			for x,line in enumerate(f[9:]):
				complex_WT = float((complex_WT1[x]+complex_WT2[x]+complex_WT3[x])/3)
				complex_Mut = float((complex_Mut1[x]+complex_Mut2[x]+complex_Mut3[x])/3)
				Complex = complex_Mut-complex_WT
				ComplexList.append(str(Complex))
			ComplexSum = 0
			Complex = "\t".join(ComplexList)
			for comp in ComplexList:
				ComplexSum = ComplexSum+float(comp)
			if mut[0] in neg and mut[-1] in pos:
				ChargeDif = 2
			if mut[0] in neg and mut[-1] not in pos:
				ChargeDif = 1
			if mut[0] not in neg and mut[0] not in pos and mut[-1] in pos:
				ChargeDif = 1
			if mut[0] in pos and mut[-1] in neg:
				ChargeDif = -2
			if mut[0] in pos and mut[-1] not in neg:
				ChargeDif = -1
			if mut[0] not in pos and mut[0] not in neg and mut[-1] in neg:
				ChargeDif = -1
			netcharge_mut = netcharge+ChargeDif
			summaryString = RoundString+'\t'+mut+'\t'+str(netcharge_mut)+'\t'+str(ChargeDif)+'\t'+ddG+'\t'+str(ComplexSum)+'\t'+Complex+'\n'
			g.write(summaryString)
			if path == dirs[0]:
				ddG_min = float(ddG)
			if float(ddG) < ddG_min:
				if ComplexSum < 0.1:
					print 'We have a winner:\n'+summaryString
					ddG_min = float(ddG)
					best_mut = mut
					best_path = path
					best_summaryString = summaryString
					best_netcharge = netcharge_mut
					best_complex = ComplexSum
	g.close()
	summary.write(best_summaryString)
	netcharge = best_netcharge
	subprocess.call('cp '+superchargefolder+'/'+RoundString+'/'+best_mut+'/SequenceDetail_'+name_pdb.split('.')[0]+'_1_0.fxout '+startingDir+'/Repair/SequenceDetail_'+RoundString+'_supercharged.fxout',shell=True)
	subprocess.call('cp '+superchargefolder+'/'+RoundString+'/'+best_mut+'/'+name_pdb.split('.')[0]+'_1_0.pdb '+startingDir+'/Repair/'+RoundString+'_supercharged.pdb',shell=True)
	print 'End of round :'+str(round)
	round = round + 1
	if best_complex > 0.1:
                print 'This one disturbs the Complex interactions: '
                print best_summaryString
                print 'Execution halted....'
                break
	os.chdir('./..')
summary.close()
os.chdir(startingDir)
