clc;
clear;

nii=load_nii('Template_T1_IXI555_MNI152.nii');
img=nii.img;
img_f=img(end:-1:1,:,:,:);
img_sy=(img+img_f)./2;
nii.img=img_sy;
save_nii(nii,'Template_T1_IXI555_MNI152_sym.nii');