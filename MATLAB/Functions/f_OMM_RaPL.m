function [fx,dF_dX,dF_dTheta] = f_OMM_RaPL(X,th,u,inF)
% OMM with piecewise-linear GA function
% -----------------------------
nf = 3;

%% -----------------------
% Toolbox
if nargin == 4
    
dt = inF.dt;
V = inF.V;
Gb = inF.Gb;
A = inF.A;
al = inF.alpha;
tb = inF.tb;

nb = length(inF.tb)-2;

t = u(1);
I = u(2);
Rap = u(3);

Ra = RaPL(t,tb,th,A,al);

fx      = [-X(2)*X(1) - exp(th(1))*(X(1)-Gb) + (Ra+Rap)/V; 
          -exp(th(2))*(X(2) - exp(th(3))*I)];

fx      = dt.*fx + X;

dx1dx1 = -exp(th(1))-X(2);                  dx1dx2 = -X(1);
dx2dx1 = 0;                                 dx2dx2 = -exp(th(2));
 
J       = [dx1dx1, dx1dx2
           dx2dx1, dx2dx2];

dF_dX   = dt*J' + eye(2);

dF_dTheta = [-exp(th(1))*(X(1)-Gb),                     0
             0,                                         -exp(th(2))*(X(2) - exp(th(3))*I)
             0,                                         exp(th(2)+th(3))*I
             zeros(nb,1)                                zeros(nb,1)];
         

dF_dTheta = get_dRa(dF_dTheta,t,th,tb,A,al,V);
dF_dTheta = dt*dF_dTheta;

end 

%% --------------------
% Get Ra
if nargin == 3

t = X;
inF = u;

fx = RaPL(t,inF.tb,th,inF.A,inF.alpha);

end

end


function [ y,A ] = RaPL( t,tb,th,A,al )
%PIECEWISE_LINEAR Summary of this function goes here
%   Detailed explanation goes here

k = exp(th(4:end));
s = (-2*A*al + al*(-k(2)*tb(2) + k(1)*tb(3) - k(3)*tb(3) + k(2)*tb(4) - k(4)*tb(4) + k(3)*tb(5) ...
    -k(5)*tb(5) + k(4)*tb(6) + k(5)*tb(7)) + k(6)*(2 - al*tb(7) + al*tb(8)))/(al*(tb(6)-tb(8)));

k = [k(1:5);s;k(6:end)];

A = 0;
k = [0;k];

for i=2:length(tb)
    if t>=tb(i-1) && t<=tb(i)
        y = k(i-1)+(k(i)-k(i-1))/(tb(i)-tb(i-1))*(t-tb(i-1));
    end
    A = A+ 0.5*(tb(i)-tb(i-1))*(k(i)+k(i-1));
end

if t>tb(end)
    y = k(end)*exp(-(t-tb(end))*al);
end

A = A+k(end)/al;

end

function dF_dTheta = get_dRa(dF_dTheta,t,th,tb,A,al,V)

nf = 3;
nb = length(tb)-2; 

k = exp(th(4:end));
 
if t>=tb(1) && t<=tb(2)
    dF_dTheta(nf+1,1) = k(1)*t/tb(2)/V;
end
 
if t>=tb(2) && t<=tb(3)
    dF_dTheta(nf+1:nb+nf,1) = [k(1)*(tb(3)-t)/(tb(3)-tb(2))/V;
        k(2)*(t-tb(2))/(tb(3)-tb(2))/V; 0; 0; 0; 0];
end

if t>=tb(3) && t<=tb(4)
    dF_dTheta(nf+1:nb+nf,1) = [0; k(2)*(tb(4)-t)/(tb(4)-tb(3))/V;
        k(3)*(t-tb(3))/(tb(4)-tb(3))/V; 0; 0; 0];
end

if t>=tb(4) && t<=tb(5)
    dF_dTheta(nf+1:nb+nf,1) = [0; 0; k(3)*(tb(5)-t)/(tb(5)-tb(4))/V;
        k(4)*(t-tb(4))/(tb(5)-tb(4))/V; 0; 0];
end

if t>=tb(5) && t<=tb(6)
    dF_dTheta(nf+1:nb+nf,1) = [0; 0; 0; k(4)*(tb(6)-t)/(tb(6)-tb(5))/V;
        k(5)*(t-tb(5))/(tb(6)-tb(5))/V; 0];
end
      
if t>=tb(6) && t<=tb(7)
     dF_dTheta(nf+1:nb+nf,1) = [k(1)*tb(3)*(t-tb(6))/(tb(7)-tb(6))/(tb(6)-tb(8));
         k(2)*(tb(2)-tb(4))*(t-tb(6))/(tb(6)-tb(7))/(tb(6)-tb(8));
         k(3)*(tb(3)-tb(5))*(t-tb(6))/(tb(6)-tb(7))/(tb(6)-tb(8));
         -k(4)*(tb(6)-tb(4))*(t-tb(6))/(tb(6)-tb(7))/(tb(6)-tb(8));
         k(5)*(-tb(5)*tb(6)+t*(tb(5)+tb(6)-tb(7)-tb(8))+tb(7)*tb(8))/(tb(6)-tb(7))/(tb(6)-tb(8));
         k(6)*(t-tb(6))*(2+al*(tb(8)-tb(7)))/(tb(7)-tb(6))/(tb(6)-tb(8))/al]/V;
end

if t>=tb(7) && t<=tb(8)
     dF_dTheta(nf+1:nb+nf,1) = [k(1)*tb(3)*(tb(8)-t)/(tb(8)-tb(7))/(tb(6)-tb(8));
         k(2)*(tb(2)-tb(4))*(t-tb(8))/(tb(6)-tb(8))/(tb(8)-tb(7));
         k(3)*(tb(3)-tb(5))*(t-tb(8))/(tb(6)-tb(8))/(tb(8)-tb(7));
         k(4)*(tb(4)-tb(6))*(t-tb(8))/(tb(6)-tb(8))/(tb(8)-tb(7));
         k(5)*(tb(5)-tb(7))*(t-tb(8))/(tb(7)-tb(8))/(tb(8)-tb(6));
         k(6)*(t*(-2+al*(tb(6)+tb(7)-2*tb(8))) + 2*tb(8) + al*(-tb(6)*tb(7)+tb(8)^2))/(tb(6)-tb(8))/(tb(8)-tb(7))/al]/V;
end

end

                       