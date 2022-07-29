# FCT

This spider will preprocess fMRI data as well as corresponding T1 data, and calculated functional correlational tensor (FCT) for each subject.

This spider is currently designed for three databases (ADNI_23, BLSA and OASIS-3). The pipeline is under progress. There will be updates afterwards.

## Inputs: 

Configuration file (inptpara.mat) which includes: 

root_dir: root directory of the project

projID: ID of the related dataset be used

pathi: path for input directory

patho: path for output directory

## Output:

FCT tensor saved as 5D mat

# References

[1] Ding, Zhaohua, et al. "Visualizing functional pathways in the human brain using correlation tensors and magnetic resonance imaging." Magnetic resonance imaging 34.1 (2016): 8-17.

[2] Ding, Zhaohua, et al. "Spatio-temporal correlation tensors reveal functional structure in human brain." PloS one 8.12 (2013): e82107.
