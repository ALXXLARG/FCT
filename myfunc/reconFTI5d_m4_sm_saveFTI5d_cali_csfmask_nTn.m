% Function to reconstruct FT and save FT (method 1)
% inputs: 
%   T1_img - 3D T1 matrix      (orientation has to be RAS) 
%   GM_img - 3D GM mask matrix (orientation has to be RAS)
%   WM_img - 3D WM mask matrix (orientation has to be RAS)
%   fmri4d - 4D fMRI matrix    (orientation has to be RAS)
%   reso - 1x3 vector of spacial resolution in mm
%   nhood - size of cube neighborhood (3 or 5)
%   rpower - power of r
%   which_view - specify the view for FT visualization(1-axial, 2-coronal, 3 sagittal)
%   first_slice/last_slice/slice_interval - specify slices for FT visualization 
%   tensor_size - tensor size for FT visualization 
%   prefix/FW/saveDir - name strings and directory for saved fname 

function fti5d= reconFTI5d_m4_sm_saveFTI5d_cali_csfmask_nTn(T1_img, GM_img, WM_img, csfmask,fmri4d,nhood, rpower, prefix, FW, saveDir,calipara)

% remove NaN and Inf from fMRI data
fmri4d(isnan(fmri4d)) = 0;  
fmri4d(isinf(fmri4d)) = 0;  

% make a mask to specify where FT will be calculated
fmri_mask = (WM_img + GM_img) > 1; % fmri_mask = WM_img > 1;


rsize = size(fmri4d, 1);
csize = size(fmri4d, 2);
ssize = size(fmri4d, 3);
fprintf('fMRI data size: %d %d %d\n', rsize, csize, ssize);
fprintf('Anatomical size: %d %d %d\n', size(T1_img, 1), size(T1_img, 2), size(T1_img, 3));


% compute functional tensor
if nhood==3
    nhood_grid = -1:1;
    startingpoint = 2;
elseif nhood==5
    nhood_grid = -2:2;
    startingpoint = 3;
elseif nhood==7
    nhood_grid = -3:3;
    startingpoint = 4;
elseif nhood==9
    nhood_grid = -4:4;
    startingpoint = 5;
else
    disp('nhood has to be 3 or 5!');
end

fti5d = zeros(rsize, csize, ssize, 3, 3);
cc = zeros(rsize,csize,ssize);

for r=startingpoint:rsize-(startingpoint-1)
    disp(r);
    for c=startingpoint:csize-(startingpoint-1)
        for s=startingpoint:ssize-(startingpoint-1)
            
            if  ((fmri_mask(r, c, s)<1)||(csfmask(r, c, s)>0))  
                continue;  
            end
            
            % define directions for covariance in RAS system

            
            cijtotal = [];
            X = [];
            for dr=nhood_grid
                for dc=nhood_grid
                    for ds=nhood_grid
                        if  ((dr==0 && dc==0 && ds==0) || (fmri_mask(r+dr, c+dc, s+ds)<1)...
                                ||(csfmask(r+dr, c+dc, s+ds)>0))
                            continue;  
                        end
                        
                           vec1 = [dr dc ds];
%                         vec1 = 10*vec1/norm(vec1);
%                         
%                         dr1= vec1(1);
%                         dc1= vec1(2);
%                         ds1 = vec1(3);


                        
                        
%                         dr1 = dr/sqrt(dr^2+dc^2+ds^2);
%                         dc1 = dc/sqrt(dr^2+dc^2+ds^2);
%                         ds1 = ds/sqrt(dr^2+dc^2+ds^2);
                        
%                         Apara = [dr1^2 2*dr1*dc1 2*dr1*ds1 dc1^2 2*dc1*ds1 ds1^2];
                        Apara = [dr^2 2*dr*dc 2*dr*ds dc^2 2*dc*ds ds^2];
                        
                        X = [X;Apara];
                        loc_tmp=corr(squeeze(fmri4d(r, c, s, :)), squeeze(fmri4d(r+dr, c+dc, s+ds, :)));
                        %%%%%
                        %%%%%calibration process: modify the corr_update
                        x00 = abs(dr);
                        y00 = abs(dc);
                        z00 = abs(ds);
                        
%                         RR = sqrt((x00^2+y00^2+z00^2)/(x00^2+y00^2+(1/(calipara^2))*z00^2));
                        RR = sqrt((x00^2+y00^2+z00^2)/(x00^2+(1/(calipara^2))*y00^2+z00^2));
                        loc_tmp = RR*loc_tmp;

                        %%%%%
%                         loc_tmp1 = loc_tmp/norm(vec1);
%                         loc_tmp1 = loc_tmp*(norm(vec1))^2;
                        %%%%%
                           loc_tmp1 = loc_tmp;
                        %%%%%
                        
                       cijtotal = [cijtotal; abs((loc_tmp1))^rpower];
                       %cijtotal = [cijtotal; abs((loc_tmp1))^rpower*(loc_tmp1>=0)];
                        
                        
                    end
                end
            end
            
             b= regress(cijtotal,X);
%%%%%%%%%%%%%%%%%%
%            b = fgls(X,double(cijtotal),'intercept',false);
            %%%%%%%%%%%%%
            if length(b)~=6
                error('CHECK!');
            end
            
            if ~isreal(b)
                error('CHECK REAL!')
            end
            
            ftensor = [b(1) b(2) b(3);b(2) b(4) b(5); b(3) b(5) b(6)];
            
            %%%
            %%%
            
            % plot FTI
            ftensor(isnan(ftensor)) = 0;  % remove NaN
            ftensor(isinf(ftensor)) = 0;  % remove Inf
            
            %%% extra 
            if isempty(ftensor)==1
                continue;
            end
            %%%
            
            %%%%
            
            %%%%
            
            fti5d(r, c, s, :, :) = ftensor;
        end % sfor
    end % cfor
end % rfor

save([saveDir 'FTI5d_' prefix '_' FW '_m4_r' num2str(rpower) '_nhood' num2str(nhood) '.mat'], 'fti5d');

return;

