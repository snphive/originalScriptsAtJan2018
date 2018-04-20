import os,sys,glob,subprocess,time,re

startingDir = os.getcwd()

if not os.path.exists('Summary'):
	os.makedirs('Summary')

aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H','H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

mapping = open('mapping.txt','r').readlines()
agad_header = open('./Agadir/1_A/PSX_globalresidue.out','r').readlines()[0].strip('\n')
agad_header = agad_header.replace('\tTANGO\t','\taa\tTANGO\t')
seqdet_header = open('./SequenceDetails/SequenceDetail_RepairPDB_1.fxout','r').readlines()[8].strip('\n')

final_header = 'PdbName\t'+agad_header+'\t'+seqdet_header+'\tGAP\tAPRnumber\n'
final = open('./Summary/FinalSummary.txt','w')
final.write(final_header)

lengths = open('./Summary/Lengths.txt','w')
lengths.write('Pdb\tMol\tAAs\n')

totalAgad = 0
totalSeqdet = 0
totalMatch = 0
totalGaps = 0
totalIncomplete = 0
APRcount = 0
lenprint = len(mapping[1:])
print float(lenprint)
EndOfAPR = 0
for w,mapline in enumerate(mapping):
	print "%.2f" % float(float(float(w)/float(lenprint))*100)
	pieces_mapping = mapline.split('\t')
	pdbraw = pieces_mapping[0]
	molsAll = []
	agads = sorted(glob.glob('./Agadir/'+pdbraw+'_*'))
	for agad in agads:
		moltempje = agad.split('_')[1]
		molsAll.append(moltempje)
	mols_seqdet = molsAll
	mols_agad = molsAll
	#print pdbmapped
	for molnum,mol in enumerate(mols_seqdet):
		agad_bool = False
		seqdet_bool = False
		agad_res = './Agadir/'+pdbraw+'_'+mols_agad[molnum]+'/PSX_globalresidue.out'
		seqdet = './SequenceDetails/SequenceDetail_RepairPDB_'+pdbraw+'.fxout'
		if os.path.exists(agad_res):
			agad_res_lines = open(agad_res,'r').readlines()[1:]
			agad_bool = True
		else:
			#print agad_res
		if os.path.exists(seqdet):
			seqdet_lines = open(seqdet,'r').readlines()[9:]
			seqdet_lines_mol = []
			max_gapsize = 0
			for seqdet_line in seqdet_lines:
				pieces_seqdet = seqdet_line.split('\t')
				if 'GAP' in seqdet_line or len(pieces_seqdet)<3:
					max_gapsize += 1
					if 'GAP' in seqdet_line:
						totalGaps += 1
					if len(pieces_seqdet)<3:
						totalIncomplete += 1
					continue
				if pieces_seqdet[2] == mol:
					seqdet_lines_mol.append(seqdet_line)
				if pieces_seqdet[2] not in mols_seqdet:
					#print seqdet_line
					max_gapsize += 1
			seqdet_bool = True
		else:
			#print seqdet
		if seqdet_bool == True and agad_bool == True:
			lengths.write(pdbraw+'\t'+mol+'\t'+str(len(seqdet_lines_mol))+'\n')
			agad_x = 0
			outofrange = False
			matches = 0
			totalAgad = totalAgad + len(agad_res_lines)
			totalSeqdet = totalSeqdet + len(seqdet_lines_mol)
			APR = False
			GAPinAPR = 0
			APRnew = False
			for agad_res_line in agad_res_lines:
				res1_agad = agad_res_line.split('\t')[2]
				resnum_agad = agad_res_line.split('\t')[1]
				TANGOscore = float(agad_res_line.split('\t')[3])
				if TANGOscore > 5:
					if APR == False:
						APR = True
						GAPinAPR = 0
						APRnew = True
					else:
						APRnew = False
				if TANGOscore < 5 and APR == True:
					EndOfAPR = 1
					APR = False
					#APRcount += 1
				else:
					EndOfAPR = 0
				try:
					res1_seqdet = aa_dict_3to1[seqdet_lines_mol[agad_x].split('\t')[1]]
				except:
					outofrange = True
					continue
				if not res1_agad == res1_seqdet:
					if APR == True:
						GAPinAPR += 1
					continue
				if outofrange == False:
					if APRnew == True:
						APRcount += 1
					if APR == True:
						APRcountTemp = APRcount
					else:
						APRcountTemp = 0
					if APRcountTemp > 0:
						final.write(mapline.split('\t')[0]+'\t'+agad_res_line.strip('\n')+'\t'+seqdet_lines_mol[agad_x].strip('\n')+'\t'+str(GAPinAPR)+'\t'+str(APRcountTemp)+'\t'+str(EndOfAPR)+'\n')
					matches += 1
					totalMatch += 1
					agad_x += 1
					GAPinAPR = 0
				else:
					break
#print str(len(agad_res_lines))+'\t'+str(len(seqdet_lines_mol))+'\t'+str(matches)
print 'Total Agadlines:\t'+str(totalAgad)
print 'Total SequenceDetailLines:\t'+str(totalSeqdet)
print 'Total Matches:\t'+str(totalMatch)
print 'Total Gaps:\t'+str(totalGaps)
print 'Total Incomplete SeqdetLines:\t'+str(totalIncomplete)
final.close()
lengths.close()
