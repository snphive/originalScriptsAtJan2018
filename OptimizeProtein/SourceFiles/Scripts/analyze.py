import os, sys, glob, subprocess

startingDir = os.getcwd()
targets = sorted(glob.glob('./PDBs/*.pdb'))
figs = 'Figures'
results = 'Results'
scripts = 'Scripts'
agadir = 'Agadir'
agadirfolders_start = sorted(glob.glob(Agadir + '/*'))
for pdb in targets:
    name = pdb.split('/')[-1].split('.')[0]
    agadirfolders_eind = sorted(glob.glob(results + '/' + name + '/' + agadir + '/*'))
    agadfolders = []
    for agad_start in agadirfolders_start:
        if name == agad_start.split('/')[-1].split('_')[0]:
            agadfolders.append(agad_start)
    total_pdb = 0
    totals = []
    for agad in agadfolders:
        f = open(agad + '/PSX_globaltotal.out', 'r').readlines()
        header_pieces = f[0].split()
        name_with_mol = agad.split('/')[-1]
        mol = name_with_mol.split('_')[-1]
        for x, program in enumerate(header_pieces):
            if 'TANGO' == program:
                index_program = x
        total = f[1].split()[index_program]
        totals.append(total)
        total_pdb = total_pdb + float(total)
