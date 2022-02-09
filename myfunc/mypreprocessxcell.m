function [xcellmatrix, I, J, K, colindex,size01] = mypreprocessxcell( xcell , xcellmask)


size01 = size(xcell);

[xrows, xcols, xplanes] = size(xcellmask);
colindex = zeros(xrows, xcols, xplanes);
xcellmatrix = [];
I = [];
J = [];
K = [];
colcounter = 0;
for i=1:xrows
    for j=1:xcols
        for k=1:xplanes
            if xcellmask(i,j,k) > 0
                tmp = squeeze(xcell(i,j,k,:));
                tmp = tmp - mean(tmp);
                if norm(tmp) > 0
                   colcounter = colcounter + 1;
                   colindex(i,j,k) = colcounter;
                   xcellmatrix = [xcellmatrix tmp];
                   I = [I i];
                   J = [J j];
                   K = [K k];
                end
            end
        end
    end
end

end
