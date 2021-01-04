import matplotlib.pyplot as plt
import tensorflow as tf 
import argparse
import pickle
import numpy as np

def plot_distogram(args):

	f = tf.io.gfile.GFile(args.distance, 'rb')	
	contact_dict = pickle.load(f, encoding='latin1')

	dist_bins = np.array(np.linspace(2.0, 20.0, num=64).tolist() ).astype(np.float32)

	dist_probs = contact_dict['probs']
	idx = np.argmax(dist_probs[:,:,:],axis=2)

	dist_mode = [dist_bins[x] for i,x in enumerate(idx)]

	fig, axs = plt.subplots( nrows=1, ncols=1 ) 
	axs.imshow(dist_mode,cmap='viridis_r') 

	#plt.colorbar(dist_bins,orientation='vertical',ax=axs)
	#plt.suptitle('ProSPr '+domain_id+' Prediction', fontsize=14)
	#plt.title('Network '+network+', Stride '+stride, fontsize=10)	
	plt.xlabel('Residue i')
	plt.ylabel('Residue j')

	fig.savefig(args.out,bbox_inches='tight')


if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('--distance', '-d', help='Distance probability pickle file', required=True)
	parser.add_argument('--out', '-o', help='Rosetta contraint file', required=True)

	args = parser.parse_args()
	plot_distogram(args)


