# MSM-Design
Multi-Step Protein Design Process.
Involves Multivariate Linear Regression Model Training followed by Brute for Optimization of a conformational state or pair of conformational states
Requires Rosetta Protein Free Energy Scoring Software (script included to score sequences)
Optionally the Amber Force Field from OpenMM can be used (Requires OpenMM python API) a script to generate Amber Energies is also included
Other Requirements include PDB files that serve as the conformational state(s) of the protein being studied.
At least 2 PDB files are required for the package, although ideally one would want to use more. 

Design workflow is as follows:

Generate PDB files through MD simulations of the specified target protein
Optional - Run N_Body_Select.pair_select() from Tools.py to generate a list of close pairwise interactions for model training

Next select two of the following ways to train the models
Individual Model Training uses a list of sequences and energies that are generated with Rosetta or OpenMM/Amber before model training
Online Model Training scores sequences with Rosetta or OpenMM/Amber and simultaneously updates the model. (Require an energy scoring script which is provided)

Individual Model Training - Run SGD_Static.static_model() from Msm_Design.py - (This requires a path to the pdb file, a path to the sequence/energy file (this file should have sequences and energies separated by a comma. Each new seq/E pair should be put on a new line in the text file), a sequence alignment file (optional but useful and should follow the format of the provided sequence alignment file), the wild type sequence, and an optional list pairwise interactions (output of Run N_Body_Select.pair_select())). There are also 3 regression types included, regular linear regression, L1 regularization (Lasso) and L2 Regularization (Ridge). Output from this method is a list of the parameter weights from the linear regression model

Online Model Training - Run SGD_Online.online_model() from Msm_Design.py - (This requires an energy scoring function (python script provided either 'score_sequence_amber.py' or 'score_sequence_rosetta.py', a path to the pdb file, a sequence alignment file (optional but useful and should follow the format of the provided sequence alignment file), the wild type sequence, and an optional list pairwise interactions (output of Run N_Body_Select.pair_select())). There are also 3 regression types included, regular linear regression, L1 regularization (Lasso) and L2 Regularization (Ridge). Output from this method is a text file where each of the parameter weights is output on a new line.

Next we run Optimization on the models
Single State Optimization - Run Sequence_Optimization.run_opt() from Msm_Design.py
