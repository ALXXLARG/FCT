function out = cat_long_biascorr(job)
%cat_long_biascorr. Longitudinal bias correction by average segmentation.
% 
%  job      .. SPM job structure
%   .images .. cell of realigend images
%   .p0     .. cell with the avg p0 label map 
%   .str    .. strength of correction (0=soft, 0.5=default, 1=strong corr.)
%   .prefix .. filename prefix (default 'm')
%   .fs     .. filter size in mm (2=strong, 4=default, 8=soft correction)
%              (if it is empty then it is defined by job.str)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_long_biascorr.m 1944 2022-02-07 17:49:38Z dahnke $

% ToDo:
%  * Integration and Test of WMHs.
%  * Strong correction (str = 0.75) caused GM overestimation and adaptation  
%    of filter thresholds and smoothin is required. 
%  * Iterative correction with low to high filter size

% RD202010: First tests showed clear improvements of the timepoints but the
%           whole pipeline seems to be less affected.
%           Hence, corrections are maybe more relevant for plasticity
%           studies or in case of artifacts.
%           Strong correction (str = 0.75) caused GM overestimation,
%           whereas low correction (str = 0.25) seemed to be much better.


  def.images    = {};
  def.segment   = {};
  def.str       = 0.5; 
  def.prefix    = 'm';
  def.hardcore  = 0;                          % direct correction based on the T1 maps > not realy working
  def.LASstr    = 0.5;                        % use (local) intensity normalization    > 
  job           = cat_io_checkinopt(job,def);
  job.str       = max(0.1,min(1,job.str)); 
  
  if ~isfield(job,'fs') ||  isempty(job.fs)  
    job.fs = 2^(4 * (1 - job.str));
  end
  %%
  if job.LASstr > 0
    fprintf('  Biasstr = %0.2f, filtersize = %0.2f, LASstr = %0.2f\n',job.str,job.fs,job.LASstr);
  else
    fprintf('  Biasstr = %0.2f, filtersize = %0.2f \n',job.str,job.fs);
  end
  
  if job.hardcore
    %% 202201 intensity based correction > remove if not further needed in 202301
    %  This is just a temporary test if stronger changes based on the original 
    %  images intens can improve the preprocessing, especially in case of 
    %  different protocols. 
    %  This can be seen as some very strong bias correction that changes
    %  also contrast differences (that why we cannot smooth stonger),
    %  whereas zero smooting would simply replace the image by the average.
    %  Although, there are some visual improvements, the is much to simple
    %  and an intensity normalization and more smoothing are probably
    %  necessary to stabilize the effect. 
    [pp,ff,ee] = spm_fileparts(job.segment{1});
    if strcmp(pp(end-2:end),'mri'), pp = pp(1:end-3); end
    Vavg   = spm_vol( fullfile(pp,[ff(3:end) ee]) ); 
    Yavg   = spm_read_vols(Vavg); 
  end
  
  % load segment and estimate voxel size
  Vp0    = spm_vol(job.segment{1}); 
  Yp0    = spm_read_vols(Vp0); 
  vx_vol = sqrt(sum(Vp0.mat(1:3,1:3).^2));

  % apply bias correction for all images
  for ii = 1:numel(job.images)
    Vo = spm_vol(job.images{ii}); 
    Yo = single(spm_read_vols(Vo)); 
  
    %% avoid edges (artifacts) in the original image
    Yg   = cat_vol_grad(Yo)  ./ (1+Yo);    % normalized gradient map of the original image  
    Ygp0 = cat_vol_grad(Yp0) ./ (1+Yp0);   % normalized gradient map of the average segmentation map
    Ydiv = cat_vol_div(Yo)   ./ Yo; 

    
    %% hard intensity correction (this is too simple)
    %Ywa = cat_vol_smooth3X(Yavg ./ Yo,0); % correct addtive to adapt avg to tp 
    if job.hardcore
      fprintf('Hard intensity correction based on the average:\n'); 
      Ybg = cat_vol_smooth3X( Yo ./ cat_stat_kmeans(Yo(round(Yp0)==3)) <0.5 & Yp0==0 , 4) ; 
      Yi  = max(0.7,max(eps,Yo) ./ max(eps,Yavg)).*(1-Ybg) + Ybg;
      Yi  = cat_vol_median3(Yi); 
      Yw  = cat_vol_smooth3X(Yi, 2 ); 

      % quantify difference
      Ygow = cat_vol_grad(Yo./Yw)  ./ (1+(Yo./Yw));
      fprintf('  Yg overlap after hard corection: %0.4f\n',cat_stat_nanmean(abs(Ygow(Yp0(:)>0) - Yg(Yp0(:)>0)))); 

      Yo = Yo ./ Yw; 
    end
    
    
    %% intensity normalization 
    T3th  = [ min(Yo(:)) min(Yo(Yp0(:)>0.5 & Yg(:)<0.2))+eps ...
              cat_stat_kmeans(Yo(round(Yp0)==1 & Yg<0.2)) ...
              cat_stat_kmeans(Yo(round(Yp0)==2 & Yg<0.2)) ...
              cat_stat_kmeans(Yo(round(Yp0)==3 & Yg<0.2)) ...
              max(Yo(Yp0>0  & Yg<0.2))];
    T3thx = [ 0 min(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_kmeans(Yo(round(Yp0)==1 & Yg<0.2))+eps ...
              1 2 3 ...
              max(Yo(Yp0(:)>0.5 & Yg(:)<0.2))/cat_stat_kmeans(Yo(round(Yp0)==3  & Yg<0.2))*3 ];

    % intensity scaling
    Ym = Yo; 
    for i=2:numel(T3th)
      M = Yo>T3th(i-1) & Yo<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Yo(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Yo>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Yo(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    

    
    %% avoid regions that strongly changed between the time points
    Yd = single(abs(Ym - Yp0) .* (Yp0>0));
    Yd = cat_vol_median3(Yd,Yp0>0,Yp0>0,0.1); 
    Yd = cat_vol_smooth3X(Yd,0.5); 
    Yd = Yd .* Yg; 
    

    %% extract save segment
    gstr = 2^( 2 * (job.str - 0.5) ); 
    Ywm = abs(Ydiv)<0.2*mean(vx_vol) & Yd<0.2*mean(vx_vol) & Yg<0.6*mean(vx_vol) * gstr & Ygp0<0.6*mean(vx_vol) & Yp0>2.5; % & cat_vol_morph(Yp0>2,'e',1,vx_vol);
    Ygm = abs(Ydiv)<0.2*mean(vx_vol) & Yd<0.2*mean(vx_vol) & Yg<0.6*mean(vx_vol) * gstr & Ygp0<0.6*mean(vx_vol) & Yp0>1.5 & Yp0<2.5; % & cat_vol_morph(Yp0>1,'e') & cat_vol_morph(Yp0<3,'e'); 
    Ycm = abs(Ydiv)<0.2*mean(vx_vol) & Yd<0.2*mean(vx_vol) & Yg<0.6*mean(vx_vol) * gstr & Ygp0<1.2*mean(vx_vol) & Yp0>0.5 & Yp0<1.5 & cat_vol_morph(Yp0>0,'e',3,vx_vol) & cat_vol_morph(Yp0<2,'de',3,vx_vol);
    Ycm(smooth3(Ycm)<0.3) = 0; 
    Ybg = Ym<0.5 & single(cat_vol_morph(Yp0==0,'de',40,vx_vol)); 
    %clear Yg; 
    
    
    %% get values from segment
    Yi = Ybg + Yo .* Ywm ./ cat_stat_kmeans(Yo(Ywm(:)));
    Yi = Yi  + Yo .* Ygm ./ cat_stat_kmeans(Yo(Ygm(:)));
    Yi = Yi  + Yo .* Ycm ./ cat_stat_kmeans(Yo(Ycm(:)));
    if 1
      if T3th(4)<T3th(5) % T1
        Yi = cat_vol_localstat(Yi,Yi~=0,1,3,round(2 ./ mean(vx_vol)));
      else
        Yi = cat_vol_localstat(Yi,Yi~=0,1,2,round(2 ./ mean(vx_vol)));
      end
    end
    % remove outlier 
    Yi2 = Yi; Yi2( smooth3(Yi2)<0.8 ) = 0;  
    Yi2 = cat_vol_localstat(Yi2,Yi2~=0,1,1,round(2 ./ mean(vx_vol)));
    Yw  = cat_vol_approx(Yi2,'nh',vx_vol,4); 
    Yi( abs(Yi-Yw)>0.2 * max(0.5,job.str) ) = 0; 
    % remove background
    Yi(Ybg) = 0; 
    Yi = cat_vol_median3(Yi,Yi~=1 & Yi~=0,Yi~=1 & Yi~=0);
    Yi = cat_vol_localstat(Yi,Yi~=0,round(3/mean(vx_vol)),1,round(10 ./ mean(vx_vol)));
    % approximate bias field
    Yw = cat_vol_approx(Yi,'nh',vx_vol,job.fs); % overcorrection of subcortical structures
    Yw = Yw ./ cat_stat_kmeans(Yw(Ywm(:)));
    
    
    
    %% local intensity correction based on the average segmentation 
    %  The priciple idea is to do some further corrections in case of
    %  changed protocols (indicated by the user). The hope is that this
    %  reduce the bias like local differences related to contrast
    %  differences. 
    %  For test the ADNI 1.5 and 3.0 Tesla scans with 100 days differences 
    %  and normal longitudinal scans with 1 and 2 years differences are used, 
    %  where the correction should result in a more likely aging pattern in 
    %  the 100-days rescans and not worse (quite similar) changes within the 
    %  normal longitudinal images. 
    if job.LASstr
      
      % apply bias correction 
      Ysrc = Yo ./ Yw; 

if 0      
      clsdef = { 
     ... p0-  g   d   k op
        0.01  0.1 0.1 1;
        1.05  0.1 0.1 1; ... C
       ... 1.33  0.1 0.1 1; 
       ... 1.66  0.1 0.1 1; 
        1.95  0.1 0.1 1; 
        2.05  0.1 0.1 1; ... G
       ... 2.33  0.1 0.1 1; 
       ... 2.66  0.1 0.1 1; 
        2.95  0.1 0.1 1; 
        3.15  0.1 0.1 1; ... W
        };
      clsdef = { 
     ... p0-  g   d   k op
        0.50  0.1 0.1 1;  %C
        1.25  1.1 1.5 1; % G
        1.75  1.1 1.1 1; % C
        2.25  1.1 1.1 1; % G
        2.75  1.3 1.1 1; % W
        3.10  0.3 0.5 1; % W
        };
      for ci = 2:size(clsdef,1)
        Ymsk = Yp0>=clsdef{ci-1,1} & Yp0<clsdef{ci,1} & Yg<clsdef{ci,2} & abs(Ydiv)<clsdef{ci,3};
        
        [srcm,strs,strn] = cat_stat_kmeans(Ysrc(Ymsk(:)),clsdef{ci,4});
        [avgm,avgs,avgn] = cat_stat_kmeans(Yavg(Ymsk(:)),clsdef{ci,4});
        
        Tth.T3th(ci)  = srcm(strn==max(strn));
        Tth.T3thx(ci) = avgm(avgn==max(avgn));
      end  
     
      if 0
        T3thx = interp1(Tth.T3thx,1:0.1:numel(Tth.T3thx),'spline');
        T3th  = interp1(Tth.T3th ,1:0.1:numel(Tth.T3th) ,'spline');
      else
        T3thx = Tth.T3thx;
        T3th  = Tth.T3th;
      end
      
      [T3th,si] = sort(T3th);
      T3thx     = T3thx(si);
      
    %  T3th  = interp1(T3thx,1:0.5:numel(Tth.T3thx),'spline');
    %  T3thx = interp1(T3th ,1:0.5:numel(Tth.T3th) ,'spline');
     
      
      isc=1; Ym = Ysrc; Ysrc = Ysrc * 0; % * 1000; Tth.T3th = Tth.T3th * 1000; 
      for i = 2:numel(T3th)
        M = Ym>T3th(i-1) & Ym<=T3th(i);
        Ysrc(M(:)) = T3thx(i-1) + (Ym(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
      end
      M  = Ym>=T3th(end); 
      Ysrc(M(:)) = numel(T3th)/isc/6 + (Ym(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
      %%
   Yww = smooth3(Ysrc)./smooth3(Yavg) .* (Yp0>0) .* (Yg<0.3 | abs(Yd)>0.1); Yww = cat_vol_median3(Yww,Yp0>0); 
   %Yww = cat_vol_localstat(Yww,Yp0>0,1,1,10); 
   Yww = cat_vol_approx(Yww,1);
   Yww = Yww ./ ( cat_stat_nanmean(Ysrc(:)./Yavg(:)./Yww(:) .* (Yp0(:)>0)) ./ cat_stat_nanmean(Ysrc(:)./Yavg(:) .* (Yp0(:)>0)));
   %   Ysrc2 = Ysrc / 1000; 
  
 %  Yn = (Ysrc ./ Yww  -  Yavg); 
 %  cat_sanlm(Yn,1,3); 
   
  %%
      Ym = cat_main_gintnormi(Ysrc,Tth);
end

      %%
      
      % estimate threshholds
      T3th2 = zeros(1,3); 
      for ti=1:3, T3th2(ti) = cat_stat_kmeans(Ysrc(round(Yp0)==ti & Yg<0.3)); end
    
        
      
      
      %% the most important segment is the GM 
      Ygm2    = Ygp0<0.6*mean(vx_vol) & Yp0>1.5 & Yp0<2.5 & Yg<0.6*mean(vx_vol); 
      Ygm2(smooth3(Ygm2)<0.3) = 0; % remove small dots
      Yi      = Ysrc .* Ygm2;
      Yi      = cat_vol_median3(Yi,Yi~=1 & Yi~=0,Yi~=1 & Yi~=0);
      Yi      = cat_vol_localstat(Yi,Yi~=0,1,1,round(20 ./ mean(vx_vol)));
      % first approximation to remove local outlier
      Ylab{1} = cat_vol_approx(Yi,'nh',vx_vol,job.fs * 4 / job.LASstr); 
      Yi( abs( log( Yi ./ Ylab{1} )) > 0.2 ) = 0; % arbitrary value between 0.10 (remove more) and 0.25 (remove less)  
      % final approximation 
      Ylab{1} = cat_vol_approx(Yi,'nh',vx_vol,job.fs * 4 / job.LASstr); 
      % for all other segments we just use the global values
      Ylab{2} = T3th2(3);
      Ylab{3} = T3th2(1); 
      Ylab{6} = min(Yo); % ####### inoptimal#

      % intensity normalization 
      Yml = zeros(size(Ysrc)); 
      Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc - Ylab{2}) ./ max(eps,Ylab{2} - Ylab{1})) );
      Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc - Ylab{1}) ./ max(eps,Ylab{2} - Ylab{1})) );
      Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc - Ylab{3}) ./ max(eps,Ylab{1} - Ylab{3})) );
      Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc - Ylab{6}) ./ max(eps,Ylab{3} - Ylab{6})) );
      Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
    end
    
    %% create some final measurements 
    %  CJV ? 
    %  RMSE ? 
    %  COV ?
    %  hist overlap? > rmse? 
    
    
    %% write corrected output
    out.bc{ii} = spm_file(Vo.fname,'prefix',job.prefix); 
    Vw = Vo; Vw.fname = out.bc{ii};
    if job.LASstr
      spm_write_vol(Vw,Yml);
    else
      spm_write_vol(Vw,Yo ./ Yw);
    end
  end
end