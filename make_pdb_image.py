import sys
import argparse
sys.path.append('/home/ml2529/alphafold/pymol/lib/python3.6/site-packages')

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI

import pymol


pymol.finish_launching()


'''
The command pymol.cmd.dss('all') doesn't seem to work so instead get secondary structure prediction
from AlphaFold NN prediction
'''
def add_sec_structure(filename):

	with open(filename, 'r') as f:
		content = f.readlines()
		rows = [ c.strip() for c in content if not c.startswith('#') and c != '\n' ]


	ss_start = 1
	ss_run = ""

	for i in range(len(rows)):
		fields = rows[i].split()

		#print("%s\t%s\t%s" %(fields[0],fields[1],fields[2]))

		if fields[2] == 'H' and ss_run != 'H':
			if ss_run != '':
				residues = "%s-%s/" % (ss_start,i)
				ss = "ss='%s'" % ss_run
				#print("%s,%s" % (residues,ss))		
				pymol.cmd.alter(residues, ss)

			ss_run = 'H'
			ss_start = i+1


		elif fields[2] == 'E' and ss_run != 'S':
			if ss_run != '':
				residues = "%s-%s/" % (ss_start,i)
				ss = "ss='%s'" % ss_run
				#print("%s,%s" % (residues,ss))		
				pymol.cmd.alter(residues, ss)

			ss_run = 'S'
			ss_start = i+1
		
		elif fields[2] != 'H' and fields[2] != 'E':
			if ss_run != '':
				residues = "%s-%s/" % (ss_start,i)
				ss = "ss='%s'" % ss_run
				#print("%s,%s" % (residues,ss))		
				pymol.cmd.alter(residues, ss)

			ss_run = ''
			ss_start = i+1

	if ss_run == 'H' or ss_run == 'S':
		residues = "%s-%s/" % (ss_start,len(rows))
		ss = "ss='%s'" % ss_run
		#print("%s,%s" % (residues,ss))
		pymol.cmd.alter(residues, ss)
			
	pymol.cmd.rebuild()



parser = argparse.ArgumentParser()

parser.add_argument('--pdb', '-p', help='pdb file', required=True)
parser.add_argument('--sec', '-s', help='secondary structure prediction file', required=False)
args = parser.parse_args()

pdb_file = args.pdb
pdb_name = pdb_file.split('.')[0]

pymol.cmd.load(pdb_file, pdb_name)
pymol.cmd.disable("all")
pymol.cmd.enable(pdb_name)

if args.sec:
	add_sec_structure(args.sec)

#print(pymol.cmd.get_names())
#pymol.cmd.dss('all')

#pymol.cmd.hide('all')
pymol.cmd.show('cartoon')

#pymol.cmd.set('ray_opaque_background', 0)
#pymol.cmd.color('red', 'ss h')
#pymol.cmd.color('yellow', 'ss s')
#pymol.cmd.show('cartoon')
#pymol.cmd.color('orange','all')

#pymol.cmd.cartoon('automatic')

#pymol.cmd.rebuild()
#pymol.cmd.show('cartoon')



#pymol.cmd.set('cartoon', 'automatic')


pymol.cmd.png("%s.png"%(pdb_name))
#pymol.cmd.save("check.pdb")
pymol.cmd.quit()

