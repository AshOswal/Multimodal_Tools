# Multimodal_Tools
 analysis tools

Tools for analysing data. Below are the subdirectories and the relevant files they contain. Please note that the scripts are written in MATLAB (2019a) and that they should work with any recent installation of SPM (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and FieldTrip
(http://www.fieldtriptoolbox.org/). 

1) Computational model contains code for generating figures pertaining to the computational model. PD_sim_master_coh will generate figure 7 which is the simulation of power and coherence. Run PD_sim_master to generate the subplots within Figures 6 and Supplementary Figure 6. Finally master_bif generates the numerical bifurcation plot in Figure 6B. 

2) Delays contains a script that simulates and estimates phase delays in the high and low beta frequency ranges. This sort of analysis has been used to compute delays in Figure 4C and also Supplementary Figure 3. See also Oswal et al., 2016 for a validation of this approach in similar datasets (https://pubmed.ncbi.nlm.nih.gov/27017189/). 

3) VoxelWiseCorrelation conatins an exemplar script detailing how nifti format images from two different modalities (in this case tractography and MEG) may be correlated at a voxel wise level. The folder contains example nifti images from tractography and MEG derived coherence. The master script to run is master_voxelwise_GLM.

4) Sharpness contains exemplar code for estimating time series sharpness and non-linear features. There is a separate README file in this directory which contains some example data and guidance on how to use scripts (read me DON and Sharpness). 

A Oswal
Please contact me if you have any further questions, feedback or comments
ashwini.oswal@ndcn.ox.ac.uk
