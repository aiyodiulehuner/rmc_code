Monotone Retargeted Preference Completion
-----------------------------------------

Contact: Suriya Gunasekar (suriya@ttic.edu)

This is a reference implementation of the algorithm described in the paper:

Suriya Gunasekar, Oluwasanmi Koyejo, and Joydeep Ghosh, "Preference Completion from Partial Rankings" NIPS 2016.

It is being made available without any guarantees for non-commercial research use only. If you find this code useful in your research, please consider citing the paper.

------------------------------------------

smc.m: function implementing standard matrix completion
rmc_fixed_margin.m: function implementing the proposed monorone retargetted preference completion algorithm. 

movienlens.m and neurosynth.m scripts to run smc and rmc on movielens and neurosynth datasets (Warning: has some hardcoded numbers and paths)

Movielens/:  Folder with scripts to process data and results of Movielens dataset
Neurosynth/: Folder with scripts to process data and results of Neurosynth dataset


