import pickle
import tensorflow as tf  
import numpy as np
from scipy.stats import vonmises
import os
import argparse

def CalcTorsionStats(torsion_prob):

	#glob_prob = 0
	phi_prob = []

	for i in range(36):
		row_prob = 0
		for j in range(36):
			row_prob += torsion_prob[i*36+j]

		#glob_prob += row_prob
		#print("i:%d row_prob: %f" % (i,row_prob))
		phi_prob.append(row_prob)

	#print(phi)
	#print("glob_prob: %f" % glob_prob)

	#glob_prob = 0
	psi_prob =[]

	for i in range(36):
		col_prob = 0
		for j in range(36):
			col_prob += torsion_prob[i+j*36]

		#glob_prob += col_prob
		#print("i:%d col_prob: %f" % (i,col_prob))
		psi_prob.append(col_prob)

	#print(psi)		
	#print("glob_prob: %f" % glob_prob)

	angle_bins = np.array(np.linspace(-1*np.pi, np.pi, num=36).tolist() ).astype(np.float32)
	#print(angle_bins)
	#print(len(angle_bins))

	phi_dist = []

	for i in range(36):
		l = [angle_bins[i]]*int(phi_prob[i]*1000)
		phi_dist.extend(l)

	psi_dist = []

	for i in range(36):
		l = [angle_bins[i]]*int(psi_prob[i]*1000)
		psi_dist.extend(l)


	#print(phi_dist)
	#r = vonmises.rvs(3, size=1000)
	#print(r)

	phi_kappa, phi_loc, phi_scale = vonmises.fit(phi_dist, fscale=1)
	phi_var = vonmises.var(phi_kappa,phi_loc,phi_scale)
	
	#print((kappa,loc,scale))
	#print(vonmises.mean(kappa,loc,scale))
	#print(vonmises.var(kappa,loc,scale))


	psi_kappa, psi_loc, psi_scale = vonmises.fit(psi_dist, fscale=1)
	psi_var = vonmises.var(psi_kappa,psi_loc,psi_scale)
	#print((kappa,loc,scale))
	#print(vonmises.mean(kappa,loc,scale))
	#print(vonmises.var(kappa,loc,scale))

	return (phi_loc,psi_loc,phi_var,psi_var)

def save_torsions_stats(filename, torsion_stats, sequence):
  """Save Torsions to a file as pickle of a dict."""
  #filename = os.path.join(torsions_dir, filebase + '.torsions')

  t_dict = dict(torsion_stats=torsion_stats, sequence=sequence)
  with tf.io.gfile.GFile(filename, 'w') as fh:
    pickle.dump(t_dict, fh, protocol=2)

def parse_torsion_pickle(args):

	f = tf.io.gfile.GFile(args.input, 'rb')
	
	contact_dict = pickle.load(f, encoding='latin1')
	#print(contact_dict.keys())

	torsion_probs = contact_dict['probs']
	sequence = contact_dict['sequence']

	#print(sequence)
	#print(torsion_probs.shape)
	#print(len(torsion_prob[0]))

	#check_prob = torsion_prob[0]
	PhiPsiList = []

	for i in range(torsion_probs.shape[0]):
		print("Calculating torsion stats for residue %d" %(i+1))
		phipsi_stat = CalcTorsionStats(torsion_probs[i])
		print(phipsi_stat)
		PhiPsiList.append(phipsi_stat)
		#PhiPsiList.append(CalcTorsionStats(torsion_probs[i]))

	#print(PhiPsiList)
	save_torsions_stats(args.out,PhiPsiList,sequence)


	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()

	parser.add_argument('--input', '-i', help='torsion pickle file', required=True)
	parser.add_argument('--out', '-o', help='output pickle file for torsion stats', required=True)

	args = parser.parse_args()

	#parse_torsion_pickle('test_output/T0955/torsion/0/torsions/T0955.torsions')
	parse_torsion_pickle(args)

