load coh_imagesG.mat         % coherence nifti images
load DMsub                   % Design matrix
load tract_imagesG.mat       % tract density images
correlate3Dimages(tract_imagesG,coh_imagesG,DMsub,'HB_',{'R_beta_mask.nii','L_beta_mask.nii'},10,1)
