function [rel_rms, abs_rms] = y_FD_Jenkinson(RealignmentParameterFile,ReferenceImage)
% function [rel_rms, abs_rms] = y_FD_Jenkinson(RealignmentParameterFile,ReferenceImage)
% Calculate FD Jenkinson (relative RMS) and absolute RMS based on SPM's realignment parameters
% Reference: Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841.
%            Jenkinson, M. 1999. Measuring transformation error by RMS deviation. Internal Technical Report TR99MJ1, FMRIB Centre, University of Oxford. Available at www.fmrib.ox.ac.uk/analysis/techrep for downloading.
% Input:
% 	RealignmentParameterFile  -   The realignment parameter file for a given participant generated by SPM. E.g., rp***.txt
%   ReferenceImage            -   The reference image for realignment (usually the first time point (one-pass) or the mean image after an initial motion correction (two-pass))
% Output:
%	rel_rms      -   relative RMS (FD Jenkinson)
%	abs_rms      -   absolute RMS
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120930.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


rmax = 80.0; %The default radius (as in FSL) of a sphere represents the brain

RP=load(RealignmentParameterFile);
nTimePoint=size(RP,1);
sinq1=sin(RP(:,4));
sinq2=sin(RP(:,5));
sinq3=sin(RP(:,6));
cosq1=cos(RP(:,4));
cosq2=cos(RP(:,5));
cosq3=cos(RP(:,6));

[RefData, RefHead] = y_Read(ReferenceImage,1);
center = RefHead.mat*([0.5*(RefHead.dim(1));0.5*(RefHead.dim(2));0.5*(RefHead.dim(3));1]);
center = center(1:3); %Get the coordinate for the center

abs_rms = zeros(nTimePoint,1);
for t=1:nTimePoint

    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    MA1=eye(4);
    MA2=(M_RigidBodyTransform);
    
    M = MA1*inv(MA2) - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    abs_rms(t) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
end


rel_rms = zeros(nTimePoint-1,1);
for t=2:nTimePoint
    M1=[1       0        0     0;...
        0    cosq1(t)  sinq1(t)  0;...
        0    -sinq1(t) cosq1(t)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t)  0    sinq2(t)     0;...
        0        1       0        0;...
        -sinq2(t) 0    cosq2(t)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t)   sinq3(t)   0     0;...
        -sinq3(t)  cosq3(t)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t,1);...
        0    1     0     RP(t,2);...
        0    0     1     RP(t,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform=MT*M1*M2*M3;
    
    
    M1=[1       0        0     0;...
        0    cosq1(t-1)  sinq1(t-1)  0;...
        0    -sinq1(t-1) cosq1(t-1)  0;...
        0       0        0     1;];
    
    M2=[cosq2(t-1)  0    sinq2(t-1)     0;...
        0        1       0        0;...
        -sinq2(t-1) 0    cosq2(t-1)     0;...
        0       0        0        1;];
    
    M3=[cosq3(t-1)   sinq3(t-1)   0     0;...
        -sinq3(t-1)  cosq3(t-1)   0     0;...
        0           0       1     0;...
        0           0       0     1;];
    
    MT=[1    0     0     RP(t-1,1);...
        0    1     0     RP(t-1,2);...
        0    0     1     RP(t-1,3);...
        0    0     0     1;];
    
    M_RigidBodyTransform_1=MT*M1*M2*M3;
    
    MA1=(M_RigidBodyTransform_1);
    MA2=(M_RigidBodyTransform);
    
    M = MA1*inv(MA2) - eye(4);
    
    A = M(1:3,1:3);
    
    T = M(1:3,4);
    
    rel_rms(t-1) = sqrt(rmax*rmax/5*trace(A'*A) + (T+A*center)'*(T+A*center));
    
end

rel_rms=[0;rel_rms]; %The FD_Jenkinson at time point t means the movement from time point t-1 to time point t. (Put the FD_Jenkinson for the first time point to "0".)



