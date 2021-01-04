
import pickle
import tensorflow as tf  
import numpy as np
from scipy.stats import vonmises
import os


#{'min_range': 2.0, 'max_range': 20.0, 'num_bins': 64

distCutoffs = np.array( [0] + np.linspace(2.0, 20.0, num=62).tolist() ).astype(np.float32)
eps = np.finfo(np.float32).eps

## d needs to be positive, cannot be -1
## cutoffs is the distance boundary array
## return the largest index position such that cutoffs[position]<=d and  d<cutoffs[position+1]
def LabelsOfOneDistance(d, cutoffs):
    result = np.digitize(np.array([d]), cutoffs) - 1
    return np.int16(result[0])

def SelectCB(AA, bUseAlternativeCB=True):
	assert len(AA) == 1

	if AA.upper() == 'G':
		if bUseAlternativeCB:
			return 'CA'
		else:
			return ''
	else:
		return 'CB'


def SelectAtomPair(sequence, i, j, atomPairType):

	'''
	if atomPairType == 'CaCa':
		return 'CA', 'CA'

	if atomPairType == 'NO':
		return 'N', 'O'
	'''

	if atomPairType == 'CbCb':
		a1 = SelectCB(sequence[i])
		a2 = SelectCB(sequence[j])
		return a1, a2

	'''
	if atomPairType == 'CaCg':
		a1 = 'CA'
		a2 = SelectCG(sequence[j])
		return a1, a2

	if atomPairType == 'CgCg':
		a1 = SelectCG(sequence[i])
		a2 = SelectCG(sequence[j])
		return a1, a2
	'''

	return None


def CalcPotentialByDFIRE(predDistMatrix, alpha=1.61, largestDistance=18, useWeight=False, minPotential=-30, maxPotential=30):
	
	#potentials = dict()
	## validProbs saves the prob of one atom/residue pair likely have valid coordinates
	#validProbs = dict()

	cutoff = distCutoffs

	rc = min(cutoff[-1], largestDistance) - 0.001
	print("Last distance bin value: %f" % rc)

	#highest index before the rc (last distance bin)
	rc_index = LabelsOfOneDistance(rc, cutoff)
	#print(cutoff[rc_index])
	print("Index of last distance bin: %d with value: %f" % (rc_index,cutoff[rc_index]))


	binwidths = [ d2 - d1 for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]
	bincenters = [ (d2 + d1)/2. for d1, d2 in zip(cutoff[:-1], cutoff[1:]) ]

	## calculate reference potential defined as alpha*log (r/rc) + log(\delta r/ \delta rc)
	## \delta(r) is binwidths and r is the bincenters

	
	refPot = alpha * np.log( bincenters / bincenters[rc_index]) + np.log( binwidths / binwidths[rc_index] )

	## idx is the index for a bin
	def CalcApproxRefPot(idx=0):
		points = np.arange(cutoff[idx] + 0.5/2, cutoff[idx+1], 0.5)
		values = np.power(points / bincenters[rc_index], alpha)
		avg = np.average(values)
		tmpRefPot = np.log(avg) + np.log( binwidths[idx] / binwidths[rc_index] )
		return tmpRefPot

	for i in range(len(binwidths)):
		if binwidths[i] >= 1:
			refPot[i] = CalcApproxRefPot(i)

	## calculate the observed potential defined as log( p(r) /p(rc) ) where p(r) is the predicted distance probability
	predProb = predDistMatrix
	predProbRC = predProb[:, :, rc_index : rc_index+1]
	obsPot = np.log(predProb / predProbRC)

	## calculate the final potential, which is the difference between reference potential and observed potential
	potential = np.zeros_like(predDistMatrix)
	potential[:, :, :rc_index ] = refPot[: rc_index] - obsPot[:, :, :rc_index]

	print("Observed potential")
	print(obsPot.shape)


	print("Predicted probability")
	print(predProb.shape)
	print(predProb[0][1])


	#Sum of probability if the last value is removed - final valid probabilities
	validProb = 1 - predProb[:, :, -1]
	#validProb = np.ones((predProb.shape[0], predProb.shape[1]), dtype=np.float32)

	print("Valid probability")
	print(validProb.shape)
	print(validProb[0][1])


	#normalize potential based on valid probability
	potential *= validProb[:, :, np.newaxis]
	#print(potential[0][1])

	#Remove last value of the potential - final distance potential
	potential = potential[:, :, :-1]
	
	#potentials[response] = potential.astype(np.float32)
	#validProbs[response] = validProb.astype(np.float32)
	print("Potential shape")
	print(potential.shape)
	#return potentials, validProbs
	return potential


## generate distance potential constraints for PyRosetta
## currently topRatio and potThreshold are not used for distance potential
def GenerateSplinePotential4Distance(sequence, pot, minSeqSep=1):

	#target, sequence, potential, distCutoffs = potData[:4]
	allConstraints = []

	x = distCutoffs
	print(x)
	print(len(x))
	print(pot.shape)

	binWidths = [ b-a for a, b in zip(x[1:-1], x[2:]) ]
	binWidth = np.average(binWidths)

	## here we add repulsion to reduce steric clashes
	## for CaCa and CbCb, the minimum distance is 3.6A. The penalty is 3 in [2, 3.6] and 10 in [0, 2]

	firstMinDist = 2
	secondMinDist = 3.6
	yPenalty = [10, 4, 0.5]

	xPrefix = [ 0, firstMinDist, secondMinDist ]

	## find the index of the 2nd min distance in x, i.e., x[secondLabel] <=secondMinDist < x[secondLabel+1]
	secondLabel = LabelsOfOneDistance(secondMinDist + 0.0001, x)
		
	#print 'secondLabel=', secondLabel
	assert secondLabel >= 1
	assert secondLabel < len(distCutoffs)

	print(secondLabel)

	xk = [ (a+b)/2 for a, b in zip(x[secondLabel:-1], x[secondLabel+1:]) ]
	xk.append(x[-1] + binWidth/2.)
	xk = xPrefix + xk
	#print(len(xk))


	size = pot.shape
	residuePairs = []

	for i in range(size[0]):
		jstart = i+minSeqSep
		
		#if not IsSymmetricLabel(labelName):
			#jstart=0

		for j in range(jstart, size[1]):
			offset = abs(i-j)
			if offset < minSeqSep:
				continue
			residuePairs.append( (i, j) )


	#print(size)
	#print(range(size[0]))
	#print(residuePairs)

	
	for i, j in residuePairs:
		y = pot[i, j]

		"""
		## y[0] is the potential for the first interval [0, x[0]). We increase potential for distance < x[0] for every binWidth Angstrom
		yPrefix2 = [ y[0] + ye for ye in yPrefix ]
		yk = yPrefix2 + y[1:].tolist()
		"""
		yPrefix = [ max(y[secondLabel], 0) + ye for ye in yPenalty ]
		y2 = y.tolist()
		yk = yPrefix + y2[secondLabel:]
		#print(yk)
		#print(len(yk))
		#print(len(xk))

		#assert len(xk) == len(yk), 'xk and yk length does not match for ' + labelName + ' and residues ' + str(i) + ' ' + str(j)
		assert len(xk) == len(yk), 'xk and yk length does not match for residues ' + str(i) + ' ' + str(j)


		## when one atom pair is not symmetric (e.g., NO), it appears twice in the constraint set, so we divide its potential by 2
		#if not IsSymmetricLabel(labelName):
			#yk = [ ye/2. for ye in yk]

		atom1, atom2 = SelectAtomPair(sequence, i, j, 'CbCb')

		constraint = dict()
		constraint['x'] = xk
		constraint['y'] = yk
		constraint['response'] = 'CbCb'
		constraint['binWidth'] = binWidth

		constraint['type'] = 'AtomPair'
		constraint['atoms'] = [ atom1, atom2]
		constraint['atomNums'] = [ i+1, j+1]

		#print(constraint)
		allConstraints.append(constraint)

	
	return allConstraints
	

## this function writes the constraints into Rosetta format
## target is the protein name
## constraints is a list of python dict and each dict corresponds to one constraint
def WriteSplineConstraints(constraints, savefile=None, savefolder4histfile=None):
	if savefile is None:
		print('ERROR: please specify the save file for constaints!')
		exit(1)
		
	if savefolder4histfile is None:
		print('ERROR: please specify the save file for constaints!')
		exit(1)

		
	histfileDir = savefolder4histfile
	if not os.path.isdir(histfileDir):
		os.mkdir(histfileDir)

	expVal = 0.
	weight = 1.

	numIgnored = 0

	potStrs = []
	for constraint in constraints:
		## write histogram to histfile
		#response = constraint['response']
		#labelName, _, _ = ParseResponse(response)
		response = constraint['response']
		labelName = response

		x = constraint['x']
		y = constraint['y']
		if not np.isfinite(y).all():
			print('WARNING: ignore one constraint since it may have an NaN or infinite value:' % (constraint))
			numIgnored += 1
			continue
			
		atomNums = [ str(i) for i in constraint['atomNums'] ]
		atomNumStr = '-'.join(atomNums)

		histfile = os.path.join(histfileDir, response + '-' + atomNumStr + '.potential.txt')
		xStr = '\t'.join(['x_axis'] + [ "{:.4f}".format(e) for e in x ] )
		yStr = '\t'.join(['y_axis'] + [ "{:.4f}".format(e) for e in y ] )
		with open(histfile, 'w') as fh:
			fh.write('\n'.join([xStr, yStr]) + '\n')

                #potStr = ' '.join(['Angle', atom1.upper(), str(i+1), atom2.upper(), str(i+2), atom3.upper(), str(j+1), 'SPLINE', description, histfile] + [ "{:.4f}".format(e) for e in [expVal, weight, binWidth] ] )
		potStrList = [ constraint['type'] ]
		for name, number in zip(constraint['atoms'], atomNums):
			potStrList.extend([name.upper(), number])
		potStrList.append('SPLINE')
		potStrList.append(labelName)
		potStrList.append(histfile)

		potStrList.extend( ['0', '1', "{:.6f}".format(constraint['binWidth']) ])
		potStr = ' '.join(potStrList)

		potStrs.append(potStr)

	if numIgnored > 100:
		print('ERROR: too many constraints are ignored:'%(numIgnored))
		exit(1)

	if len(potStrs)>0:
        	with open(savefile, 'w') as fh:
        		fh.write('\n'.join(potStrs) + '\n')

	return potStrs


def GeneratePhiPsiPotential(sequence, PhiPsiList, funcType='AMBERPERIODIC', weight0=1, predDisorder=None):

	constraints = []

	for i, PhiPsi in zip(range(len(sequence)), PhiPsiList):

		if predDisorder is not None:
			weight = weight0 * (1-predDisorder[i])
		else:
			weight = weight0

		## the three elements in PhiPsi are the predicted two means (phi and psi) and two variances
		## for phi
		if i > 0:
			if funcType == 'AMBERPERIODIC':
				phi_mean = '%.4f' % (PhiPsi[0] + np.pi)
			else:
				phi_mean = '%.4f' % (PhiPsi[0])


			if funcType == 'CHARMM':
				phi_sig = '%.4f' % (weight * 2./(eps + PhiPsi[2]) )
			elif funcType == 'AMBERPERIODIC':
				phi_sig = '%.4f' % (weight /(eps + PhiPsi[2]) )
			else:
				phi_sig = '%.4f' % np.sqrt(PhiPsi[2]/(eps + weight) )


			resNum1, resNum2, resNum3, resNum4 = str(i), str(i+1), str(i+1), str(i+1)
			atomName1, atomName2, atomName3, atomName4 = 'C', 'N', 'CA', 'C'
			n_periodic = '1.0'

			if funcType in [ 'HARMONIC', 'CIRCULARHAMONIC']:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, phi_mean, phi_sig])
			else:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, phi_mean, n_periodic, phi_sig])
			constraints.append(line)

        ## for psi
		if i < len(sequence)-1 :
			if funcType == 'AMBERPERIODIC':
				psi_mean = '%.4f' % (PhiPsi[1] + np.pi)
			else:
				psi_mean = '%.4f' % (PhiPsi[1])

			if funcType == 'CHARMM':
				psi_sig = '%.4f' % (weight * 2./(eps + PhiPsi[3]) )
			elif funcType == 'AMBERPERIODIC':
				psi_sig = '%.4f' % (weight /(eps + PhiPsi[3]) )
			else:
				psi_sig = '%.4f' % np.sqrt(PhiPsi[3]/(eps + weight) )

			resNum1, resNum2, resNum3, resNum4 = str(i+1), str(i+1), str(i+1), str(i+2)
			atomName1, atomName2, atomName3, atomName4 = 'N', 'CA', 'C', 'N'
			n_periodic = '1.0'
			
			if funcType in [ 'HARMONIC', 'CIRCULARHAMONIC']:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, psi_mean, psi_sig])
			else:
				line = ' '.join(['Dihedral', atomName1, resNum1, atomName2, resNum2, atomName3, resNum3, atomName4, resNum4, funcType, psi_mean, n_periodic, psi_sig])
			constraints.append(line)

	return constraints


def parse_distance_pickle_file(args):

	f = tf.io.gfile.GFile(args.distance, 'rb')
	
	contact_dict = pickle.load(f, encoding='latin1')
	
	#num_res = len(contact_dict['sequence'])
	#print("Number of residues: %d" % (num_res))
	#print(contact_dict.keys())

	#dict_keys(['min_range', 'max_range', 'num_bins', 'domain', 'sequence', 'probs'])
	#print(contact_dict)
	#print(contact_dict)

	probs = contact_dict['probs']
	#print(probs.shape)
	#print(distCutoffs)

	potential = CalcPotentialByDFIRE(probs)
	pairConstraints = GenerateSplinePotential4Distance(contact_dict['sequence'],potential)
	WriteSplineConstraints(pairConstraints, savefile=args.out, savefolder4histfile=args.hist)


def parse_torsion_pickle_stat(filename):

	f = tf.io.gfile.GFile(filename, 'rb')
	
	contact_dict = pickle.load(f, encoding='latin1')

	print(contact_dict.keys())
	PhiPsiList = contact_dict['torsion_stats']
	sequence = contact_dict['sequence']

	print(PhiPsiList)
	
	PhiPsiConstraints = GeneratePhiPsiPotential(sequence, PhiPsiList)
	print(PhiPsiConstraints)

	
	with open('test', 'a') as fh:
		fh.write('\n'.join(PhiPsiConstraints) )
	




if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('--hist', '-h', help='directory containing histogram files', required=True)
	parser.add_argument('--distance', '-d', help='Distance probability pickle file', required=True)
	parser.add_argument('--torsion', '-t', help='Distance probability pickle file', required=True)		
	parser.add_argument('--out', '-o', help='Rosetta contraint file', required=True)

	args = parser.parse_args()

	parse_pickle_file(args)


	#parse_pickle_file('test_output/T0955/distogram/ensemble/T0955.pickle')
	#parse_torsion_pickle_stat('torsion.pickle')

