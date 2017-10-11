# MSM-Design
Multi-Step Protein Design Process.
Involves Multivariate Linear Regression Model Training followed by Brute force Optimization of a conformational state or pair of conformational states.
Requires Rosetta Protein Free Energy Scoring Software (function included to score sequences in msm_design.msm_model).
Optionally the Amber Force Field from OpenMM can be used (Requires OpenMM python API) a function to generate Amber Energies is also included.
Other Requirements include PDB files that serve as the conformational state(s) of the protein being studied. Examples can be found in /data/pdb_files/. 
At least 2 PDB files are required for the package, although ideally one would want to use more. 

Design workflow is as follows:

Generate PDB files through MD simulations of the specified target protein suggested to use OpenMM and MdTraj for consistency
Naming convention should be StateX.pdb where X is a number corresponding to the specific state.
Optional - Run pair_select from msm_design.msm_model to generate a list of close pairwise interactions for model training. **Note Optimization currently does not support models trained with pair interactions.

Next select two of the following ways to train the models.
Individual Model Training uses a list of sequences and energies that are generated with Rosetta or OpenMM/Amber before model training.
Online Model Training scores sequences with Rosetta or OpenMM/Amber and simultaneously updates the model. (Require an energy scoring function which is provided and requires to be on linux).

**Again, Rosetta is Required to utilize this package, make sure you know the location of Fixbb.linuxgccrelease and the location of the Rosetta db in your linux directories
**Ideally the scoring portion could be re written using PyRosetta

Individual Model Training - Run SGD_Static.static_model() from msm_model - (This requires a path to the pdb file of interest, a path to the sequence/energy file (this file should have sequences and energies separated by a comma. Each new seq/E pair should be put on a new line in the text file), a sequence alignment file (optional but useful and should follow the format of the provided sequence alignment file in /data/other_data/), the wild type sequence, and an optional list pairwise interactions (output of pair_select from Msm_Design)). There are also 3 regression types included, regular linear regression, L1 regularization (Lasso) and L2 Regularization (Ridge). Output from this method is a list of the parameter weights from the linear regression model, optionally you can save the parameter weights as a text file. Ideally, naming convention should be wStateX.txt where X corresponds with the pdbfile number.

Online Model Training - Run SGD_Online.online_model() from msm_model - (This requires an energy scoring function (python definition in msm_model called either score_sequence_amber or score_sequence_rosetta, a path to the pdb file, a sequence alignment file (optional but useful and should follow the format of the provided sequence alignment file), the wild type sequence, and an optional list pairwise interactions (output of pair_select from msm_model.pair_select())). There are also 3 regression types included, regular linear regression, L1 regularization (Lasso) and L2 Regularization (Ridge). Output from this method is a text file where each of the parameter weights is output on a new line. Again, Ideally the name of the text file is automatically generated from the pdbfile name but naming convention should be wStateX.txt where X corresponds with the pdbfile number.

Next we run Optimization on the models
Single State Optimization - Run Sequence_Optimization.single_state_design() from msm_optimization - (This requires a directory filled with model output text files in the proper format. Each file should have a number in the file name corresponding to the conformational state that the model represents). The sequence alignment file as well as the wild type sequence are also needed. There are various options when it comes to input and output. As far as input one can change temperature, k (gas constant), the number of mutation pathways, the maximum number of allowed mutations, the number of intra pathway attempts, and of course the list of states to optimize. As far as output one can select to display plots as output or write them as pdfs, one can write an optimization text file which contains the probability, optimized state number, and best optimized sequence for all mutations allowed (ie if we allow 5 mutations then we will have a "best" sequence at 0, 1, 2, 3, 4, and 5 mutations per state selected).

Two State Optimization (Similar to the single state optimization) - Run Sequence_Optimization.two_state_design() from msm_optimization - (This requires a directory filled with model output text files in the proper format. Each file should have a number in the file name corresponding to the conformational state that the model represents). The sequence alignment file as well as the wild type sequence are also needed. There are various options when it comes to input and output. As far as input one can change temperature, k (gas constant), the number of mutation pathways, the maximum number of allowed mutations, the number of intra pathway attempts, and of course the list of pairs of states to optimize. **the input of pairs of states to optimize for MUST be a list of lists ie [[4,5],[1,7]].

As far as output one can select to display plots as output or write them to a pdf file, one can write an optimization text file which contains the probability, optimized state number, and best optimized sequence for all mutations allowed (ie if we allow 5 mutations then we will have a "best" sequence at 0, 1, 2, 3, 4, and 5 mutations).

We provide an example script that goes through training a model using the static method, and running an optimization with 10 states, optimizing for one of them. 

If you would like to use this package please do so and if you have any questions do not hesitate to contact the author
Trenth12@gmail.com
