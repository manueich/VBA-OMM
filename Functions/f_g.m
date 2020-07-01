function [gx,dG_dX,dG_dPhi] = f_g(Xt,Phi,ut,inG)
% Direct Maping g(x) = x

gx = Xt(1);
dG_dX = [1;0];
dG_dPhi = [];

end

