function H = plottensortotal1(T1_img,tensor,fmri_mask,which_view,k,FA)
rsize = size(T1_img, 1);
csize = size(T1_img, 2);
ssize = size(T1_img, 3);   
H = figure; 
    set(H, 'Position', [20, 200, rsize*8, csize*8]); set(H, 'Visible', 'On');
    cslice = round(csize/2);  rslice = round(rsize/2);  sslice=round(ssize/2);
    radd = 0.5; cadd = 0.5; sadd = 0.5;

rfor = 4:rsize-3;
cfor = 4:csize-3;
sfor = 4:ssize-3;
    switch which_view
        case 1 % axial slice
           sslice = k; sfor = k; az = 90; el = -90; sadd = -1;   
        case 2 % coronal slice
           cslice = k; cfor = k; az = 90; el = 0;   cadd =  1;   
        case 3 % sagital slice
           rslice = k; rfor = k; az = 0; el = 0;    radd = -1;  
        otherwise        
    end
    h = slice(single(T1_img), cslice, rslice, sslice);
    set(h, 'edgecolor', 'none'); colormap gray; 
    xlim([1 csize]); ylim([1 rsize]); zlim([1 ssize]); axis equal; 
    hold on
    for r = rfor
        for c= cfor
            for s = sfor
            if  (fmri_mask(r, c, s)<1)
                continue;  
            end
                ftensor = squeeze(tensor(r,c,s,:,:));
                ftensor = ftensor/(ftensor(1,1)+ftensor(2,2)+ftensor(3,3)+realmin);
                plottensor3_1(r-radd, c+cadd, s+sadd, ftensor, 12, 0,FA(r,c,s));
%                 drawnow;
            end
        end
    end


end