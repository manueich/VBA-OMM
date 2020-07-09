function [gx,dG_dX,dG_dPhi] = f_g(Xt,Phi,ut,inG)
% Observation function y(t) = G(t)

gx = Xt(1);
dG_dX = [1;0];
dG_dPhi = [];

end

