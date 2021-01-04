# AlphaFoldv1 Pipeline

The code here was copied directly over from [AlphaFoldv1 github repository](https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13).
Copied over as I didn't want to fork the entire deepmind-research reposity. At the moment many have realized that the AlphaFoldv1 code publicly shared is 
very limited and does not include code for creating the input .trec files or more importantly the folding to produce the folded structure (i.e. PDB file).
The code only produces the distance probabilities (i.e. distograms) and torsion probabilities, which form the constraints used for folding. This repository
aims to implement the missing preparation of the .trec files and the folding.

## Requirements

## Step 1: Creation of Input trec files
TO DO - I would like to create my own implementation.  


Currently there are two other implementations  
https://github.com/Urinx/alphafold_pytorch  
https://github.com/dellacortelab/prospr  

However both these implementation are again aimed at only producing the distance probabilities (i.e. distograms) and torsion probabilities.

## Step 2: Creation of distance and torsion probabilities
See the AlphaFoldv1 [Instructions](https://github.com/deepmind/deepmind-research/blob/master/alphafold_casp13/README.md) 
particular the [run_eval.sh script](https://github.com/deepmind/deepmind-research/blob/master/alphafold_casp13/run_eval.sh)  

Note: This is limited to only the .trec files created for CASP13. To run this on a protein sequence of interest you will need to prepare your own input files detailed
in step 1.

## Step 3: Folding and creating PDB files
I decided to create bare minimum code that takes the distance and torsion constraints and uses it to fold the protein. This can be done by using the pickle files
produced by AlphaFoldv1 and converting them to Rosetta constraints that are then used to fold the protein using PyRosetta.

```
#Create torsion stat files
python make_torsion_stats.py -i test_output/T0955/torsion/0/torsions/T0955.torsions -o torsion.pickle

#Create Rosetta constraint file and histogram files
python make_constraints.py --hist test_hist --distance test_output/T0955/distogram/ensemble/T0955.pickle --torsion torsion.pickle --out test

#Fold using Rosetta
python fold.py --fasta 5W9F.fasta --constraints test --out test4.pdb
```

This code that produces the constraint and folding has been adapted from the [RaptorX-3DModelling](https://github.com/j3xugit/RaptorX-3DModeling)

## TO DO


## References and Readings







