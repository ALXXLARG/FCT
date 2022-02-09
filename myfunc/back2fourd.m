function img1 = back2fourd(xcellmatrix,I,J,K,size01)

% mark as the identifier of whether it is needed to conduct the
% normalization for all the voxels

if size01(4)~=size(xcellmatrix,1)
    error('wrong size!');
end

img1 = single(zeros(size01));



for i = 1 : size(xcellmatrix,2)
    tmp = (xcellmatrix(:,i));
    if std(tmp)>0
        d = std(tmp);
        tmp = tmp - mean(tmp);
        tmp = tmp/d;
    end
    img1(I(i),J(i),K(i),:) = tmp;
end

end