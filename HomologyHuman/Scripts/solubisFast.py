# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
import os,sys,glob,subprocess,time,yasara,datetime

aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H','H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

if not os.path.exists('Solubis'):
	os.makedirs('Solubis')

#subprocess.call('cp ./RepairFiles/rotabase.txt ./Solubis/.',shell=True)
slct = open('selection.txt','r').readlines()

startingDir = os.getcwd()

os.chdir('Solubis')
if not os.path.exists('Runscripts'):
	os.makedirs('Runscripts')
if not os.path.exists('Indivs'):
	os.makedirs('Indivs')
if not os.path.exists('AverageBuild'):
	os.makedirs('AverageBuild')
if not os.path.exists('Fasta'):
	os.makedirs('Fasta')
os.chdir(startingDir)
#prog = open('ProgressSolubis.txt','w')
#prog.write('Solubis job started\n')
#pdblist = sorted(glob.glob('./RepairPDBs/*.pdb'))
#totalpdbs = len(slct)
#timestart = datetime.datetime.now()
#prog.write('Total amount of RepairPDBs to be done lies around:\t'+str(totalpdbs)+' (will be a bit less, percentage displayed will be accurate)\n')
#prog.write('Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())+'\n')
#prog.flush()
mapping = open('mapping.txt','r').readlines()[1:]

tmIDs = []
tmlines = open('TmIDs.txt','r').readlines()
for tm in tmlines:
	ID = tm.split('\t')[0]
	tmIDs.append(ID)

gatekeepers = ['R','P','K','E','D']
"""
g = open('./Solubis/job_s0.q','w')
g.write('#!/bin/bash\n')
g.write('#$ -N SB_0\n')
g.write('#$ -q all.q\n')
g.write('#$ -cwd\n')
#g.write('#$ -l h_vmem=2G\n')
g.write('source ~/.bash_profile\n')

indivAll = open('IndividualListAllSolubis.txt','w')
"""
summary = open('Summary/FinalSummary.txt','r').readlines()
splitfactor = 20
tempw = 0
totalpdbs = 0
fastas = []
#doubles = open('DuplicatedSeqs.txt','w')

"""

for mapline in slct:
	name = mapline.split('\t')[0]
	pdb = './RepairPDBs/RepairPDB_'+name+'.pdb'
	if not os.path.exists(pdb):
		continue
	if os.path.exists('./Solubis/AverageBuild/Average_BuildModel_RepairPDB_'+name+'.fxout'):
		continue
	fasta = './Fasta/'+name+'_A.fasta'
	if not os.path.exists(fasta):
		continue
	fastalines = open(fasta,'r').readlines()
	sequenceFasta = fastalines[-1].strip('\n')
	print(sequenceFasta)
	if sequenceFasta in fastas:
		doubles.write(mapline)
		continue
	else:
		fastas.append(sequenceFasta)
	print(name)
	sumlines = []
	for sumline in summary:
		pieces = sumline.split('\t')
		if pieces[0] == name:
			sumlines.append(sumline)
	if len(sumlines) == 0:
		continue
	if name in tmIDs:
                print(name)
                continue
	tempw = tempw+1
	os.chdir('Solubis')
	if tempw % splitfactor == 0:
		g.close()
		check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
		output_qstat = check_qstat.stdout.read()
		while output_qstat.count('robkan') > 120:
			print('Waiting for all Solubis jobs to finish')
			time.sleep(600)
			averages = len(glob.glob('./AverageBuild/Average_BuildModel*.fxout'))
			percentage = str(float((averages/len(slct)*100)))
			timenow = datetime.datetime.now()
			timedif = timenow-timestart
			prog.write('We are at '+percentage+'%, time is Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())+'\t after '+str(timedif)+' (estimate, waiting for room to submit jobs)\n')
			prog.flush()
			check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
			output_qstat = check_qstat.stdout.read()
		subprocess.call('qsub job_s'+str(tempw/splitfactor-1)+'.q',shell=True)
		newjobnumber = str(tempw/splitfactor)
		g = open('job_s'+newjobnumber+'.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N SB_'+newjobnumber+'\n')
		g.write('#$ -q all.q\n')
#		g.write('#$ -l h_vmem=2G\n')
		g.write('#$ -cwd\n')
		g.write('source ~/.bash_profile\n')
	subprocess.call('cp ./../'+pdb+' .',shell=True)
	indiv = open('individual_list_'+name+'.txt','w')
	run = open('runscript_'+name+'.txt','w')
	run.write('<TITLE>FOLDX_runscript;\n')
	run.write('<JOBSTART>#;\n')
	run.write('<PDBS>RepairPDB_'+name+'.pdb;\n')
	run.write('<BATCH>#;\n')
	run.write('<COMMANDS>FOLDX_commandfile;\n')
	run.write('<BuildModel>#,individual_list_'+name+'.txt;\n')
	run.write('<END>#;\n')
	run.write('<OPTIONS>FOLDX_optionfile;\n')
	run.write('<Temperature>298;\n')
	run.write('<IonStrength>0.05;\n')
	run.write('<ph>7;\n')
	run.write('<moveNeighbours>true;\n')
	run.write('<VdWDesign>2;\n')
	run.write('<numberOfRuns>1;\n')
	run.write('<OutPDB>#;\n')
	run.write('<END>#;\n')
	run.write('<JOBEND>#;\n')
	run.write('<ENDFILE>#;\n')
	run.close()
	Mols = ['A']
	#agads = sorted(glob.glob(startingDir+'/Agadir/'+name+'_*'))
	#for agad in agads:
	#	print(agad)
	#	moltempje = agad.split('_')[-1]
	#	if moltempje == "A":
	#		Mols.append(moltempje)
	#	else:
	#		continue
	nfasta = open('./Fasta/'+name+'.fasta','w')
	for mol in Mols:
		oldfasta_lines = open(startingDir+'/Fasta/'+name+'_'+mol+'.fasta','r').readlines()
		nfasta.write('>'+name+'_'+mol+'_WT\n')
		nfasta.write(oldfasta_lines[1]+'\n\n')
		oldfasta = oldfasta_lines[1].strip('\n')
		# get stretches
		mutation = ""
		for line in sumlines:
			pieces = line.split('\t')
			if pieces[14] != mol:
				continue
			name = pieces[0]
			index_agadir = pieces[2]
			aa = pieces[3]
			mol = pieces[14]
			index_mutation_foldX = pieces[15]
			for gate in gatekeepers:
				newfastalist = list(oldfasta)
				try:
					newfastalist[int(index_agadir)-1] = gate
				except:
					print(name+'\t'+mol+'\t'+index_agadir+'\t'+index_mutation_foldX+'\t'+aa+'\t'+gate+'\n')
					print(sumline)
				newfasta = "".join(newfastalist)
				mutation = aa+mol+index_mutation_foldX+gate
				newfastaname = name+'_'+mol+'_'+mutation
				nfasta.write('>'+newfastaname+'\n')
				nfasta.write(newfasta+'\n\n')
				indiv.write(mutation+';\n')
				indivAll.write(name+'\t'+mol+'\t'+mutation+'\t'+index_agadir+'\t'+index_mutation_foldX+'\t'+aa+'\t'+gate+'\n')
	nfasta.close()
	indiv.close()
	totalpdbs += 1
	if mutation != "":
		g.write('/switchlab/group/tools/FoldX_2015/FoldX -runfile runscript_'+name+'.txt\n')
		g.write('mv individual_list_'+name+'.txt ./Indivs/.\n')
	else:
		g.write('rm individual_list_'+name+'.txt\n')
	g.write('rm runscript_'+name+'.txt\n')
	if mutation != "":
		g.write('mv Average_BuildModel_RepairPDB_'+name+'.fxout ./AverageBuild/.\n')
		g.write('rm RepairPDB_'+name+'_*.pdb\n')
		g.write('rm WT_RepairPDB_'+name+'_*.pdb\n')
		g.write('rm Dif_BuildModel_RepairPDB_'+name+'.fxout\n')
        	g.write('rm PdbList_BuildModel_RepairPDB_'+name+'.fxout\n')
        	g.write('rm Raw_BuildModel_RepairPDB_'+name+'.fxout\n')	
	g.write('rm RepairPDB_'+name+'.pdb\n')
	g.write('rm BuildModel_RepairPDB_'+name+'.fxout\n')
	os.chdir(startingDir)
indivAll.close()
g.close()
os.chdir('Solubis')
subprocess.call('qsub job_s'+newjobnumber+'.q',shell=True)

agadjob = open('agadjob.q','w')
agadjob.write('#!/bin/bash\n')
agadjob.write('#$ -N agadir\n')
agadjob.write('#$ -cwd\n')
agadjob.write('#$ -V\n')
agadjob.write('#$ -l h_vmem=3G\n')
agadjob.write('#$ -q all.q\n')
agadjob.write('source ~/.bash_profile\n')
agadjob.write('python '+startingDir+'/Scripts/agadir.py\n')
agadjob.close()
subprocess.call('qsub agadjob.q',shell=True)

os.chdir(startingDir)

check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'SB_' in output_qstat:
	print('Waiting for all Solubis jobs to finish')
	time.sleep(600)
	averages = len(glob.glob('Solubis/AverageBuild/Average_BuildModel*.fxout'))
	percentage = str(float((averages/totalpdbs)*100))
	timenow = datetime.datetime.now()
	timedif = timenow-timestart
	prog.write('We are at '+percentage+'%, time is Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())+'\t after '+str(timedif)+'\n')
	prog.flush()
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()


sumshort = open('FinalSummarySolubisShort.txt','w')
sumshort.write('Pdb\tMol\tMutation\tddG\tTANGOold\tTANGOnew\tdTANGO\n')

for mapline in mapping:
	sumshort.flush()
	name = mapline.split('\t')[0]
	pdb = './RepairPDBs/RepairPDB_'+name+'.pdb'
	avbuild = './Solubis/AverageBuild/Average_BuildModel_RepairPDB_'+name+'.fxout'
	indiv = './Solubis/Indivs/individual_list_'+name+'.txt'
	agadir = './Solubis/Agadir/'+name+'/PSX_globaltotal.out'
	if not os.path.exists(pdb):
		print('No PDB: '+name)
		continue
	if not os.path.exists(avbuild):
		print('No Build: '+name)
		continue
	if not os.path.exists(indiv):
		print('No Indiv: '+name)
		continue
	if not os.path.exists(agadir):
		print('No Agadir: '+name)
		continue
	avbuild_lines = open(avbuild,'r').readlines()[9:]
	indiv_lines = open(indiv,'r').readlines()
	agadir_lines = open(agadir,'r').readlines()[1:]
	print(str(len(avbuild_lines))+'\t'+str(len(indiv_lines))+'\t'+str(len(agadir_lines)))
	WTtango = agadir_lines[0].split('\t')[2]
	agadir_lines = agadir_lines[1:]
	ddGs = []
	muts = []
	oldTANGOs = []
	newTANGOs = []
	dTANGOs = []
	for indiv_num,indiv_line in enumerate(indiv_lines):
		mut = indiv_line.strip(';\n')
		gatemut = mut[-1]
		ddG = avbuild_lines[indiv_num].split()[2]
		ddGs.append(ddG)
		muts.append(mut)
		for agadir_line in agadir_lines:
			if mut == agadir_line.split('\t')[0].split('_')[-1]:
				tango = agadir_line.split('\t')[2]
				dtango = str(float(tango)-float(WTtango))
				dTANGOs.append(dtango)
				oldTANGOs.append(WTtango)
				newTANGOs.append(tango)
	for lala,mut in enumerate(muts):
		sumshort.write(name+'\tA\t'+muts[lala]+'\t'+ddGs[lala]+'\t'+oldTANGOs[lala]+'\t'+newTANGOs[lala]+'\t'+dTANGOs[lala]+'\n')
sumshort.close()

"""
summaryfinal = open('FinalSummarySolubisExt.txt','w')
summaryfinal.write('Pdb\tMol\tRes\tTANGO\tAPRnumber\tResnumTANGO\tResnumFoldX\tMutation\tddG_R\tdTANGO_R\toldTANGO_R\tnewTANGO_R\tddG_P\tdTANGO_P\toldTANGO_P\tnewTANGO_P\tddG_K\tdTANGO_K\toldTANGO_K\tnewTANGO_K\tddG_E\tdTANGO_E\toldTANGO_E\tnewTANGO_E\tddG_D\tdTANGO_D\toldTANGO_D\tnewTANGO_D\n')

for mapline in mapping:
	summaryfinal.flush()
	name = mapline.split('\t')[0]
	pdb = './RepairPDBs/RepairPDB_'+name+'.pdb'
	avbuild = './Solubis/AverageBuild/Average_BuildModel_RepairPDB_'+name+'.fxout'
	indiv = './Solubis/Indivs/individual_list_'+name+'.txt'
	agadir = './Solubis/Agadir/'+name+'/PSX_globaltotal.out'
	if not os.path.exists(pdb):
		print('No PDB: '+name)
		continue
	if not os.path.exists(avbuild):
		print('No Build: '+name)
		continue
	if not os.path.exists(indiv):
		print('No Indiv: '+name)
		continue
	if not os.path.exists(agadir):
		print('No Agadir: '+name)
		continue
	print name
	avbuild_lines = open(avbuild,'r').readlines()[9:]
	indiv_lines = open(indiv,'r').readlines()
	agadir_lines = open(agadir,'r').readlines()[1:]
	sumlines = []
	for sumline in summary:
		pieces = sumline.split('\t')
		if pieces[0] == name and pieces[14] == "A":
			sumlines.append(sumline)
	if len(sumlines) == 0:
		continue
	WTtango = {}
	for agadline in agadir_lines:
		if '_WT' in agadline:
			mol = agadline.split('\t')[0].split('_')[1]
			tango = agadline.split('\t')[2]
			print(mol+'\t'+tango+'\t'+agadline)
			WTtango[mol] = tango
	for line in sumlines:
		pieces = line.split('\t')
		name = pieces[0]
		index_agadir = pieces[2]
		aa = pieces[3]
		mol = pieces[14]
		tangores = pieces[4]
		index_mutation_foldX = pieces[15]
		APRnumber = pieces[50]
		mutation = aa+mol+index_mutation_foldX
		ddG_R = ''
		ddG_K = ''
		ddG_P = ''
		ddG_E = ''
		ddG_D = ''
		ddGs = []
		muts = []
		for indiv_num,indiv_line in enumerate(indiv_lines):
			if mutation in indiv_line:
				mut = indiv_line.strip(';\n')
				gatemut = mut[-1]
				ddG = avbuild_lines[indiv_num].split()[2]
				#ddGs.append(ddG)
				#muts.append(mut)
				if gatemut == 'R':
					ddG_R = ddG
					continue
				if gatemut == 'K':
					ddG_K = ddG
					continue
				if gatemut == 'E':
					ddG_E = ddG
					continue
				if gatemut == 'P':
					ddG_P = ddG
					continue
				if gatemut == 'D':
					ddG_D = ddG
					continue
		dtango_R = ''
		dtango_E = ''
		dtango_K = ''
		dtango_D = ''
		dtango_P = ''
		old_tango_R = ''
		old_tango_E = ''
		old_tango_K = ''
		old_tango_D = ''
		old_tango_P = ''
		new_tango_R = ''
		new_tango_E = ''
		new_tango_K = ''
		new_tango_D = ''
		new_tango_P = ''
		#dTANGOs = []
		for agadir_line in agadir_lines:
			if mutation in agadir_line:
				mut = agadir_line.split('\t')[0].split('_')[-1]
				gatemut = mut[-1]
				#print(mut+'\t'+gatemut)
				tango = agadir_line.split('\t')[2]
				dtango = str(float(tango)-float(WTtango[mol]))
				#dTANGOs.append(dtango)
				#print('Dtango = '+dtango)
				if gatemut == 'R':
					dtango_R = dtango
					old_tango_R = WTtango[mol]
					new_tango_R = tango
					continue
				if gatemut == 'K':
					dtango_K = dtango
					old_tango_K = WTtango[mol]
					new_tango_K = tango
					continue
				if gatemut == 'E':
					dtango_E = dtango
					old_tango_E = WTtango[mol]
					new_tango_E = tango
					continue
				if gatemut == 'P':
					dtango_P = dtango
					old_tango_P = WTtango[mol]
					new_tango_P = tango
					continue
				if gatemut == 'D':
					dtango_D = dtango
					old_tango_D = WTtango[mol]
					new_tango_D = tango
					continue
		summaryfinal.write(name+'\t'+mol+'\t'+aa+'\t'+tangores+'\t'+APRnumber+'\t'+index_agadir+'\t'+index_mutation_foldX+'\t'+mutation+'\t'+
							ddG_R+'\t'+dtango_R+'\t'+old_tango_R+'\t'+new_tango_R+'\t'+
							ddG_P+'\t'+dtango_P+'\t'+old_tango_P+'\t'+new_tango_P+'\t'+
							ddG_K+'\t'+dtango_K+'\t'+old_tango_K+'\t'+new_tango_K+'\t'+
							ddG_E+'\t'+dtango_E+'\t'+old_tango_E+'\t'+new_tango_E+'\t'+
							ddG_D+'\t'+dtango_D+'\t'+old_tango_D+'\t'+new_tango_D+
							'\n')
summaryfinal.close()
