import os,sys,glob,subprocess,time

startingDir = os.getcwd()

if not os.path.exists("Summary"):
	os.makedirs("Summary")

f = open('./Summary/SequenceDetailAll.txt','w')
seqdets = glob.glob('./Results/*/Repair/SequenceDetail_RepairPDB_*.fxout')

for x,seqdet in enumerate(seqdets):
	lines = open(seqdet,'r').readlines()
	if x == 0:
		f.write('number	Pdb	RepairPDB	amino acid	chain	oldnumber	omega angle	phi angle	psi angle	Secondary Structure	total energy	Backbone Hbond	Sidechain Hbond	Van der Waals	Electrostatics	Solvation Polar	Solvation Hydrophobic	Van der Waals clashes	entropy sidechain	entropy mainchain	sloop_entropy	mloop_entropy	cis_bond	torsional clash	backbone clash	helix dipole	water bridge	disulfide	electrostatic kon	partial covalent bonds	energy Ionisation	Entropy Complex	heteroBackHbond	heteroSideHbond	sidechain burial	mainchain burial	sidechain Occ	mainchain Occ	index\n')
	chains = []
	chain = ""
	for line in lines[9:-1]:
		pieces = line.split('\t')
		chain = pieces[2]
		if chain not in chains:
			count = 0
			chains.append(chain)
		count = count + 1
		finalline = line
		finalline = finalline.replace('RepairPDB_','')
		pdb = seqdet.split('_')[-1].split('.')[0]
		f.write(str(count)+'\t'+pdb+'\t'+finalline)
f.close()

f= open('./Summary/StabilityAll.txt','w')
stabs = glob.glob('./Results/*/Repair/Stability_RepairPDB_*.fxout')
print len(stabs)

for x,stab in enumerate(stabs):
	lines = open(stab,'r').readlines()
	if x== 0:
		f.write('Pdb\tRepairPDB'+lines[8])
	finalline = lines[9]
	finalline = finalline.replace('RepairPDB_','')
	#finalline = finalline.replace('_chothia.pdb','')
	pdb = finalline.split('.')[0]
	f.write(pdb+'\t'+finalline)
f.close()

f= open('./Summary/AnalyseComplex.txt','w')
ACs = glob.glob('./Results/*/Repair/Interaction_AnalyseComplex_RepairPDB*.fxout')
print len(ACs)

for x,AC in enumerate(ACs):
	lines = open(AC,'r').readlines()
	if x== 0:
		f.write('Pdb\t'+lines[8])
	finalline = lines[9]
	finalline = finalline.replace('RepairPDB_','')
	#finalline = finalline.replace('_chothia.pdb','')
	pdb = finalline.split('.')[0]
	f.write(pdb+'\t'+finalline)
f.close()

f=open('./Summary/AgadAll.txt','w')
g=open('./Summary/Windows.txt','w')
g.write('Pdb\tPos\tNGK\tAPR\tCGK\tTANGO\tLength\n')
f.write('Pdb\tName\tMol\tResnum\tRes1\tTANGO\tAPRcount\n')
agadfolders = []
agads = glob.glob('./Results/*/Agadir/*')
for agad in agads:
	temp = agad.split('/')[-1]
	if temp == 'Options.txt':
		continue
	agadfolders.append(agad)

APRcount = 0
APR = False

for agad in agadfolders:
	print agad
	agad_lines = open(agad+'/PSX_globalresidue.out','r').readlines()[1:]
	agad_window_lines = open(agad+'/PSX_tangowindow.out','r').readlines()[1:]
	for agad_window_line in agad_window_lines:
		g.write(agad_window_line)
	for agad_res_line in agad_lines:
		pieces = agad_res_line.split('\t')
		nametemp = pieces[0]
		moltemp = nametemp.split('_')[-1]
		resnum = pieces[1]
		res1 = pieces[2]
		TANGOstring = pieces[3]
		TANGOscore = float(TANGOstring)
		if TANGOscore > 5:
			if APR == False:
				APR = True
				APRnew = True
				APRcount += 1
			else:
				APRnew = False
		if TANGOscore < 5 and APR == True:
			EndOfAPR = 1
			APR = False
		else:
			EndOfAPR = 0
		if APR == False:
			APRcountstring = '0'
		else:
			APRcountstring = str(APRcount)
		pdb = nametemp.split('_')[0]
		f.write(pdb+'\t'+nametemp+'\t'+moltemp+'\t'+resnum+'\t'+res1+'\t'+TANGOstring+'\t'+APRcountstring+'\n')
f.close()
g.close()

os.chdir('Summary')
subprocess.call('R < ./../SourceFiles/Scripts/StretchPlotsAndScoringAll.R --no-save',shell=True)
os.chdir(startingDir)
