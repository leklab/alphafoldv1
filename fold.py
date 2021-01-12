'''
Code adapted from the various files from the RaptorX-3DModeling github repository
https://github.com/j3xugit/RaptorX-3DModeling

Namely these files:
https://github.com/j3xugit/RaptorX-3DModeling/blob/master/Folding/Scripts4Rosetta/FoldNRelax.py
https://github.com/j3xugit/RaptorX-3DModeling/blob/master/Folding/Scripts4Rosetta/GeneratePairPotential4Rosetta.py
https://github.com/j3xugit/RaptorX-3DModeling/blob/master/Folding/GenPairwisePotentialFromPrediction.py
https://github.com/j3xugit/RaptorX-3DModeling/blob/master/DL4DistancePrediction4/config.py
'''

import numpy as np
import argparse

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

def ReadFASTAFile(inputFile):
	with open(inputFile, 'r') as fh:  
    		sequence = fh.readlines()
  
	# removing the trailing "\n" and any header lines
	sequence = [line.strip() for line in sequence if not '>' in line]
	sequence = ''.join( sequence )    # combine into a single sequence

	return sequence


def ExtractPhiPsiDistribution(seq, cstFile):

	assert seq is not None
	seqLen = len(seq)	
	PhiDistribution = [None] * seqLen 
	PsiDistribution = [None] * seqLen 
	rows = None

	with open(cstFile, 'r') as f:
		content = f.readlines()
		rows = [ c for c in content if 'AMBER' in c and 'Dihedral' in c ]

	if len(rows) < 5:
		print("ERROR: there are very few Phi/Psi constraints in " % (cstFile))
		exit(1)

	for r in rows:
		fields = r.split()
		assert len(fields)==13
		assert fields[0] == 'Dihedral'
		x0 = np.float32(fields[-3])
		n = np.float32(fields[-2])
		k = np.float32(fields[-1])
		idx = np.int32(fields[4])-1

		if fields[1]=='N' and fields[3]=='CA' and fields[5]=='C' and fields[7]=='N':
			## psi
			PsiDistribution[idx] = (x0, n, k)
		elif fields[1]=='C' and fields[3]=='N' and fields[5]=='CA' and fields[7]=='C':
			##phi
			PhiDistribution[idx] = (x0, n, k)
		else:
			print('WARNING: unknown dihedral type in line: ' % (r))
			continue

	distribution = dict()
	distribution['phi'] = PhiDistribution
	distribution['psi'] = PsiDistribution

	return distribution

def SampleDihedralsByAMBER(distribution):
	## generate a list of discrete angles from 0 to 360
	step = 8.
	anchors_degree = np.arange(step/2.-180, 180.-0.1, step)
	anchors = anchors_degree / 180. * np.pi
	#print anchors_degree
	#print anchors
	
	## random sample pertubations between -step/2 and step/2
	samples = np.random.uniform(-step/2, step/2, size=len(distribution)).astype(np.float32) 

	for idx, d in zip(range(len(distribution)), distribution):
		if d is None:
			samples[idx] = np.random.uniform(-180., 180.)
			continue

		x0, n, k = d
		#print idx, d
		## ll represents log-likelihood
		ll = -k * ( 1 + np.cos( (n*anchors) - x0 ) )
	
		## calculate the probability of the anchors by d
		prob = np.exp(ll)
		prob = prob / np.sum(prob)

		#print idx, prob

		## random sample one anchor by prob
		sample = np.random.choice(anchors_degree, p=prob)
		samples[idx] += sample

	return samples

def SetAngles(pose, phis, psis):
	for i, phi, psi in zip(range(1, pose.total_residue() + 1), phis, psis):
		pose.set_phi(i, phi)
		pose.set_psi(i, psi)


def InitializePose(inputFile, PhiPsiDistribution=None):

	if inputFile.endswith('.fasta') or inputFile.endswith('.seq'):
		sequence = ReadFASTAFile(inputFile)
		pose = pose_from_sequence(sequence)

		if PhiPsiDistribution is None:
			phis = np.random.uniform(-180, 180, pose.total_residue() )
			psis = np.random.uniform(-180, 180, pose.total_residue() )
			SetAngles(pose, phis, psis)
		else:
			PhiDistribution = PhiPsiDistribution['phi']
			phis = SampleDihedralsByAMBER(PhiDistribution)
			PsiDistribution = PhiPsiDistribution['psi']
			psis = SampleDihedralsByAMBER(PsiDistribution)

			SetAngles(pose, phis, psis)
	else:
		pose = pose_from_pdb(inputFile)

	return pose

def RemoveClash(scorefxn, mover, pose):
    for _ in range(0, 5):
        if float(scorefxn(pose)) < 10:
            break
        mover.apply(pose)

## pose shall already contain some constraints
def Fold(pose, ncycles=1000, tolerance=0.0001, UseNBList=True, UsePerturbation=False):
	assert pose is not None

	mmap = MoveMap()
	mmap.set_bb(True)
	mmap.set_chi(False)
	mmap.set_jump(True)
	mmap.show()

	sf = ScoreFunction()
	sf.add_weights_from_file('params/scorefxn.wts')

	sf1 = ScoreFunction()
	sf1.add_weights_from_file('params/scorefxn1.wts')

	sf_vdw = ScoreFunction()
	sf_vdw.add_weights_from_file('params/scorefxn_vdw.wts')

	sf_cart = ScoreFunction()
	sf_cart.add_weights_from_file('params/scorefxn_cart.wts')

	min_mover = MinMover(mmap, sf, 'lbfgs_armijo_nonmonotone', tolerance, True)
	min_mover.max_iter(ncycles)

	min_mover1 = MinMover(mmap, sf1, 'lbfgs_armijo_nonmonotone', tolerance, True)
	min_mover1.max_iter(ncycles)

	min_mover_vdw = MinMover(mmap, sf_vdw, 'lbfgs_armijo_nonmonotone', tolerance, True)
	min_mover_vdw.max_iter(500)

	min_mover_cart = MinMover(mmap, sf_cart, 'lbfgs_armijo_nonmonotone', tolerance, True)
	min_mover_cart.max_iter(ncycles)
	min_mover_cart.cartesian(True)

	## remove clash in the initial pose
	RemoveClash(sf_vdw, min_mover_vdw, pose)

	repeat_mover = RepeatMover(min_mover, 4)
	repeat_mover.apply(pose)

	if UsePerturbation:
        	pose = MinimizeEnergyByPerturbation(pose, min_mover, sf, sigmas=[10, 7.5, 3, 2])

	min_mover_cart.apply(pose)

	RemoveClash(sf_vdw, min_mover1, pose)

	sf.show(pose)

	switch = SwitchResidueTypeSetMover("fa_standard")
	switch.apply(pose)

	return pose



if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('--fasta', '-f', help='fasta file', required=True)
	parser.add_argument('--constraints', '-c', help='Rosetta constraints file', required=True)	
	parser.add_argument('--out', '-o', help='output pdb file for folded structure', required=True)

	args = parser.parse_args()


	seq = ReadFASTAFile(args.fasta)
	print(seq)

	PhiPsiDistribution = ExtractPhiPsiDistribution(seq,args.constraints)
	#print(PhiPsiDistribution)
	pyrosetta.init()

	pose = InitializePose(args.fasta, PhiPsiDistribution)

	if pose is None:
		print('ERROR: the intial pose is None')
		exit(1)

	switch = SwitchResidueTypeSetMover("centroid")
	switch.apply(pose)

	## read in constraints
	constraints = protocols.constraint_movers.ConstraintSetMover()
	constraints.constraint_file(args.constraints)
	constraints.add_constraints(True)
	constraints.apply(pose)


	#pose = Fold(pose, tolerance=tolerance, ncycles=ncycles, UseNBList=UseNBList, UsePerturbation=UsePerturbation)

	pose = Fold(pose, ncycles=1000)

	if pose is None:
		print('ERROR: the folded pose is None')
		exit(1)

	#DSSP = protocols.moves.DsspMover()
	#DSSP.apply(pose)

	pose.dump_pdb(args.out)








