import os,glob,subprocess
fastafolder = 'Fasta'

agadir = 'Agadir'
if not os.path.exists(agadir):
	os.makedirs(agadir)
g = open('./Agadir/Options.txt','w')
g.write('<TITLE>AGADIR_optionfile;\n')
g.write('<Temperature>298.;\n')
g.write('<pH>7.5;')
g.write('<IonStrength>0.150;\n')
g.write('<TfeConc>0.;\n')
g.write('<Stability>0.;')
g.write('<Concentration>1.;\n')
g.write('<Nterm>#;\n')
g.write('<Cterm>#;\n')
g.write('<global_total>true;\n')
g.write('<tango_window>false;\n')
g.write('<waltz_window>false;\n')
g.write('<limbo_window>false;\n')
g.write('<agadir_window>false;\n')
g.write('<casablanca_window>false;\n')
g.write('<complex_window>false;\n')
g.write('<repeat_window>false;\n')
g.write('<patentTango_window>false;\n')
g.write('<tango_residue>false;\n')
g.write('<waltz_residue>false;\n')
g.write('<limbo_residue>false;\n')
g.write('<complex_residue>false;\n')
g.write('<agadir_residue>false;\n')
g.write('<casablanca_residue>false;\n')
g.write('<repeat_residue>false;\n')
g.write('<windows_file_per_sequence>false;\n')
g.write('<residue_file_per_sequence>false;\n')
g.write('<END>\n')
g.close()
startingDir = os.getcwd()
fastafiles = sorted(glob.glob('./'+fastafolder+'/*'))
fasta_chunks = [fastafiles[i:i + 100] for i in range(0, len(fastafiles), 100)]

for x,fasta_chunk in enumerate(fasta_chunks):
        jobname = 'job_'+str(x)+'.q'
        g = open(jobname,'w')
        g.write('#!/bin/bash\n')
        g.write('#$ -N j'+str(x)+'\n')
        g.write('#$ -cwd\n')
        g.write('#$ -V\n')
        g.write('#$ -q all.q\n')
	g.write('#$ -l h_vmem=3G\n')
	g.write('source ~/.bash_profile\n')
	g.write('cd Agadir\n')
	for j in fasta_chunk:
		name = j.split('/')[-1].split('.')[0]
		print(name)
		if os.path.exists(name+'/PSX_globaltotal.out'):
			continue
		if not os.path.exists(name):
			g.write('mkdir '+name+'\n')
		g.write('cd '+name+'\n')
		g.write('agadirwrapper '+startingDir+'/'+j+' '+startingDir+'/Agadir/Options.txt\n')
		g.write('cd ./..\n')
	g.close()
	subprocess.call('qsub job_'+str(x)+'.q',shell=True)
