# FCT

This spider will preprocess fMRI data as well as corresponding T1 data, extract mean time-courses of each predefined ROI and compute the correlation matrices between white matter ROIs and gray matter ROIs. Please see Gao’s publications [1, 2] for more details. The spider will also compute FALFF, ALFF and ReHo maps. 

This XNAT spider is currently designed for three databases (ADNI_23, BLSA and OASIS-3) which are proposed to be analyzed in white matter reanalysis project (PI: Dr. Gore and Dr. Landman). 

## Inputs: 

Configuration file (inptpara.mat)
includes: 
root_dir: root directory of the project
projID: ID of the related dataset be used
pathi: path for input directory
patho: path for output directory

## Output:

FCT tensors saved as 5D mat

# References

[1] Ding, Zhaohua, et al. "Visualizing functional pathways in the human brain using correlation tensors and magnetic resonance imaging." Magnetic resonance imaging 34.1 (2016): 8-17.

[2] Ding, Zhaohua, et al. "Spatio-temporal correlation tensors reveal functional structure in human brain." PloS one 8.12 (2013): e82107.
