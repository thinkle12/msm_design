# MSM-Design
Computational Protein Engineering Masters Thesis Files
Multi Step computational project

Starts by using 'Model_Generation_SGD_Ubiquitin_4positionmutations_Server.py' coupled with 'score_sequence_rosetta_submitnode.py'. Using these two scripts we score protein energies given random mutations (Need input PDB State file) and we generate a linear regression model using Stochastic Gradient Descent.
These two scripts depend on the input PDB state file. We use Ubiqutin PDB files, If pairs/triplets are enabled then 'Model_Generation_SGD_Ubiquitin_4positionmutations_Server.py' also depends on 'coordinate_pair_selection_ubiquitin.py' or 'triplet_selection.py'. 
Once we generate Models, we then use 'Optimization.py' to optimize for a given sequence using a Boltzman Distribution (Requires multiple states each with a linear regression model). 
Once we do this we generate a candidate sequence and then use 'rosetta_score_all_structures_one_conformation.py' to generate energies of all States of the candidate solution. 
Then we use 'Optimization_Check_Rosetta_Energies.py' to compare the model energies to the true rosetta energy values.
