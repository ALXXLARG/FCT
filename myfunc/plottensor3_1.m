function [xx, yy, zz] = plottensor3_1(varargin)

% Usage: plottensor(r_center, c_center, s_center, tensor, number_of_grids,
% figure_type). If figure_type = 1, it draws a peanut plot, otherwise a
% ellipsoid.

% Parse possible Axes input
% error(nargchk(6, 6, nargin));
narginchk(6,7);
[cax, args, nargs] = axescheck(varargin{:});

[rcenter, ccenter, scenter, T, n, fig_type,FAvalue] = deal(args{1:7});
[c, r, s] = sphere(n);


if  fig_type == 1   % peanut shape
    L = r.*(T(1, 1)*r + T(1, 2)*c + T(1, 3)*s) + ...
        c.*(T(2, 1)*r + T(2, 2)*c + T(2, 3)*s) + ...
        s.*(T(3, 1)*r + T(3, 2)*c + T(3, 3)*s);
    Lr = r.*L;
    Lc = c.*L;
    Ls = s.*L;
else   % ellipsoid
    Lr = T(1, 1)*r + T(1, 2)*c + T(1, 3)*s;
    Lc = T(2, 1)*r + T(2, 2)*c + T(2, 3)*s;
    Ls = T(3, 1)*r + T(3, 2)*c + T(3, 3)*s;
end


% % colormapping 1
% red_m  = abs(Lr);
% green_m= abs(Lc);
% blue_m = abs(Ls);
% clr_3d = cat(3, red_m, green_m, blue_m);
% clr_max = max(clr_3d,[],3);
% if (max(clr_max(:))>1)
%     clr_3d_norm= clr_3d./sqrt(red_m.^2+green_m.^2+blue_m.^2);
%     red_m_norm = clr_3d_norm(:,:,1);
%     green_m_norm= clr_3d_norm(:,:,2);
%     blue_m_norm = clr_3d_norm(:,:,3);
%     red_m(clr_max>1)=red_m_norm(clr_max>1);
%     green_m(clr_max>1)=green_m_norm(clr_max>1);
%     blue_m(clr_max>1)=blue_m_norm(clr_max>1);
%     clr_3d = cat(3, red_m, green_m, blue_m);
% end

% % colormapping 2
% [V,D]=eig(T); 
% v1=abs(V(:,3)); v1_n=v1/norm(v1);
% sz=size(L); 
% clr_3d = cat(3, ones(sz)*v1_n(1), ones(sz)*v1_n(2), ones(sz)*v1_n(3));

% colormapping 3
[V,D]=eig(T); 

[D_sorted,ind] = sort(diag(D),'descend');
V_new = V(:,ind);
%
v1=abs(V_new(:,1));
v1_n=v1/norm(v1);
%
e1=D_sorted(1); e2=D_sorted(2); e3=D_sorted(3); 
%FA=sqrt(sqrt(((e1-e2)^2+(e1-e3)^2+(e2-e3)^2)/2/(e1^2+e2^2+e3^2)));
sz=size(Lr); 
clr_3d = FAvalue*cat(3, ones(sz)*v1_n(1), ones(sz)*v1_n(2), ones(sz)*v1_n(3));

% clr_3d = (clr_3d>1) + (clr_3d<=1).*clr_3d;


surf(ccenter + Lc, rcenter + Lr, scenter + Ls, clr_3d, 'EdgeColor', 'none');
axis equal;
% shading(h, 'interp');
lighting phong;
hold on;
