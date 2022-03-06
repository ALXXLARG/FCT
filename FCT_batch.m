% FCT pipeline

%%%
function FCT_batch(inptpara)


root_dir = inptpara.rootdir;
out_dir = inptpara.outdir;
projID = inptpara.projID;
rfMRIpath = inptpara.rfMRI;
T1path = inptpara.T1;

%%% add paths
if ~isdeployed
    addpath(genpath([root_dir,'\add2path\']));
    addpath(genpath([root_dir '\myfunc\']))
end

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

list=dir([out_dir,'\preprocess\']);
if length(list)>2
    rmdir([out_dir,'\preprocess\*.*'],'s');
end


mkdir([out_dir,'\preprocess\FunImg\subject']);
mkdir([out_dir,'\preprocess\T1img\subject']);

%%% copy files to the preprocess folders
copyfile(rfMRIpath,[out_dir,'\preprocess\FunImg\subject']);
copyfile(T1path,[out_dir,'\preprocess\T1img\subject']);

load(which([root_dir,'\prepro_saved.mat']),'Cfg');
cd(out_dir);
Cfg.WorkingDir=[out_dir,'\preprocess'];
Cfg.DataProcessDir=[out_dir,'\preprocess'];
list=dir([out_dir,'\preprocess\FunImg\']);
Cfg.SubjectID={};
for i=3:length(list)
    Cfg.SubjectID{i-2,1}=list(i).name;
end
%######################## set parameters #######################%
list1 = dir([out_dir,'\preprocess\FunImg\']);
list1(1:2,:)=[];
tmp = dir(strcat(out_dir,'\preprocess\FunImg\',list1(1).name,'\*.nii.gz'));
info = niftiinfo(strcat(out_dir,'\preprocess\FunImg\',list1(1).name,'\',tmp.name));
timepoints = info.ImageSize(4);
TRvalue = info.PixelDimensions(4);
slicenum = info.ImageSize(3);
%%%%%%%%set slice order%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(projID, 'ADNI_23')
    sliceorder = [1:2:slicenum 2:2:slicenum];
    if mod(slicenum,2)==1
        refslice = slicenum;
    elseif mod(slicenum,2)==0
        refslice = slicenum-1;
    end
elseif strcmp(projID, 'BLSA')
    sliceorder = 1:slicenum;
    refslice = round(slice/2);
elseif strcmp(projID, 'OASIS-3')
    sliceorder = 1:slicenum;
    refslice = round(slice/2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cfg.TimePoints=timepoints; % set number of TRs
p_removed=0; %remove first number of TRs
Cfg.TR=TRvalue; % set TR
%%%%%%%%%%slice timing or not?
Cfg.IsSliceTiming = 1;
%%%%%%%%%%
Cfg.SliceTiming.SliceNumber=slicenum; % set number of slices
Cfg.SliceTiming.SliceOrder=sliceorder; % set slice order
Cfg.SliceTiming.ReferenceSlice=refslice; % set reference slice
%###############################################################%
% list=dir([path01,'\preprocess\FunImg\']);
% list(1:2)=[];
% if p_removed>1
%     for i=1:length(list)
%         tmplist=dir([path01,'\preprocess\FunImg\',list(i).name,'\*.nii']);
%         nii=load_untouch_nii([path01,'\preprocess\FunImg\',list(i).name,'\',tmplist(1).name]);
%         nii.hdr.dime.dim(5)=nii.hdr.dime.dim(5)-p_removed;
%         img=nii.img;
%         img(:,:,:,1:p_removed)=[];
%         nii.img=img;
%         save_untouch_nii(nii,[path01,'\preprocess\FunImg\',list(i).name,'\',tmplist(1).name]);
%     end
%     Cfg.TimePoints=Cfg.TimePoints-p_removed;
% end
save([out_dir '\prepro_saved.mat'], 'Cfg');
DPARSFA_run([out_dir '\prepro_saved.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% copy files to the data folder
copyfile([out_dir, '\preprocess\FunImgARCWD\'],[out_dir '\data\'])

%
flags = struct('mask', false, 'mean', false, 'interp', 1, 'which', 1, 'wrap', [0 0 0], 'prefix', 'r');
path_brod = [root_dir '\atlas\rBrodmann_YCG.nii'];
path_mask = [out_dir '\preprocess\Masks\AllResampled_BrainMask_05_91x109x91.nii'];
%%%%%%%%%%%
spm_reslice_quiet({path_brod path_mask},flags);%use spm function 'spm_reslice_quiet.m' in SCZ-WM-pipeline by Yurui and Dylan
copyfile([out_dir '\preprocess\Masks\rAllResampled_BrainMask_05_91x109x91.nii'],...
  [out_dir '\brainmask.nii']);

%%%%%% remove all the folders in preprocess to save space
cd(out_dir);
rmdir([out_dir '\preprocess\'], 's');
%%% load data
listdata = dir([out_dir '\data\']);
listdata(1:2,:)=[];
% saveDir = [path01 '\FCT\'];
%%%brain mask

%%%%%%%%%%%%%
for i = 1 : length(listdata)
%     reslice_nii([rootDir listdata(i).name '\Detrend_4DVolume.nii'],[rootDir listdata(i).name '\rest1.nii'],[2 2 2]);
    path1 = [[out_dir '\data\'] listdata(i).name '\Detrend_4DVolume.nii'];
    spm_reslice_quiet({path_brod path1},flags);
    %%%% save dir
    mkdir([out_dir '\results\',listdata(i).name,'\']);
    saveDir = [out_dir '\results\',listdata(i).name,'\'];
    %%% set up parameters for FTI functions
%     info=niftiinfo([outDir listdata(i).name '\rDetrend_4DVolume.nii']);
    nhood = 7;
    rpower = 2;
    TR = TRvalue;
    
    %%% Step 0: preprocess the fmri
%     GMnii  = load_nii(strcat(path01,'\brainmask.nii'));
    GM  = niftiread(strcat(out_dir,'\brainmask.nii'));
    GM = flip(GM,1);
%     fMRnii = load_nii(strcat(rootDir,listdata(i).name,'\rest1.nii'));
    fMR = niftiread(strcat([out_dir '\data\'],listdata(i).name,'\rDetrend_4DVolume.nii'));
    fMR = flip(fMR,1);
    xcell_mask = GM;
    xcell = fMR;
    xcell = double(xcell);
    
    % space smooth
    FWHM = 14;
    voxelsize = 2;
    
    for time = 1 : size(xcell,4)
        xcell(:,:,:,time) = imgaussfilt3(xcell(:,:,:,time),FWHM/voxelsize/2.355);
        %     xcell(:,:,:,time) = smooth3(xcell(:,:,:,time),'gaussian',[3 15 3],0.65);
    end
    
    % generate signal matrix
    [xcellmatrix, I, J, K, colindex,size01] = mypreprocessxcell(xcell, xcell_mask);
    xcellmatrix1 = xcellmatrix;
    % regress motion and physiological parameters
    
    % remove time points
    N_delete = 0;
    xcellmatrix1(1:N_delete,:) = [];
    size01_new = [size01(1),size01(2),size01(3),size01(4)-N_delete];
    
    Fs = 1/TR;
    xcellmatrix1 = xcellmatrix1-mean(xcellmatrix1);
    xcellmatrix1 = detrend(xcellmatrix1);
    xcellmatrix2 = mylpfilt(xcellmatrix1,TR);
    
    % normalize and back to 4D
    fMR1 = back2fourd(xcellmatrix2,I,J,K,size01_new);
    y0 = 1;
    %%% calculate FCT 5D
    fti5d= reconFTI5d_m4_sm_saveFTI5d(GM, double(fMR1), nhood,rpower,saveDir,y0);
    save(strcat(saveDir,'fti5d.mat'),'fti5d');

end
%%%%%%END!!

