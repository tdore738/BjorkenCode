# BjorkenCode
Bjorken code used during Travis Dore's PhD

Most of the code here was used for writing the paper: https://arxiv.org/abs/2207.04086
Hydro_DNMR_NoCS.py is an exception to this, but I think it's a bit easier to understand first so I made some comments in the file to try and help a new user parse through it a bit.

The rest of the files work roughly as follows:
- First, make a list of rho's you would like to run using GenRhoList.py. It's not necessary to use the script, just easier
-  Then, you'll also need an EoS list. I provided one here, without the EoS files. You'll need to construct whatever EoS files you need to run first. In general, you will probably need to change a lot of the directory locations.
-  The file 'txtRead_eosLoop.sh' is the main driver for running a bunch of hydro runs on the cluster. You need to run it with two arguments: the eos list and then the rho list. This script creates the output directories and then submits the job script. 
