function [fx,dF_dX,dF_dTheta] = f_OMM_RaLN(X,th,u,inF)
% OMM with log-normal GA function
%-----------------------------
nf = 3;

%% -----------------------
% Toolbox
if nargin == 4
    
dt = inF.dt;
V = inF.V;
Gb = inF.Gb;
A = inF.A;

nb = 5;

t = u(1);
I = u(2);
Rap = u(3);

Ra = RaLN(t,th,A);

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
         

dF_dTheta = get_dRa(dF_dTheta,t,th,A,V);
dF_dTheta = dt*dF_dTheta;

end 

%% --------------------
% Get Ra
if nargin == 3

t = X;
inF = u;

[fx,dF_dX,dF_dTheta] = RaLN(t,th,inF.A);

end

end

function [y,f1,f2] = RaLN(t,th,A)

nf = 3;

if t==0
    y = 0;
    f1 = 0;
    f2 = 0;
else
    Rh = 1/(1+exp(-th(nf+5)));
    f1 = (1-Rh)*A./(t*sqrt(exp(th(nf+2))*pi)).*exp(-(log(t/exp(th(nf+1)))-exp(th(nf+2))/2).^2/exp(th(nf+2)));
    f2 = A*Rh./(t*sqrt(exp(th(nf+4))*pi)).*exp(-(log(t/exp(th(nf+3)))-exp(th(nf+4))/2).^2/exp(th(nf+4)));

    y = f1 + f2;
end

end

function dF_dTheta = get_dRa(dF_dTheta,t,th,A,V)

nf = 3;

if t==0
    dF_dTheta(nf+1:nf+5) = zeros(5,1);
else
    Rh = 1/(1+exp(-th(nf+5)));
    fu1 = (1-Rh)*A./(t*sqrt(exp(th(nf+2))*pi)).*exp(-(log(t/exp(th(nf+1)))-exp(th(nf+2))/2).^2/exp(th(nf+2)));
    fu2 = A*Rh./(t*sqrt(exp(th(nf+4))*pi)).*exp(-(log(t/exp(th(nf+3)))-exp(th(nf+4))/2).^2/exp(th(nf+4)));

    dF_dTheta(nf+1,1) = 2*fu1*exp(-th(nf+2))*(-exp(th(nf+2))/2+log(t/exp(th(nf+1))))/V;
    dF_dTheta(nf+2,1) = fu1*(-exp(th(nf+2))/2+log(t/exp(th(nf+1)))+exp(-th(nf+2))*(-exp(th(nf+2))/2+log(t/exp(th(nf+1)))).^2)/V ...
        - fu1/2/V;
    dF_dTheta(nf+3,1) = 2*fu2*exp(-th(nf+4))*(-exp(th(nf+4))/2+log(t/exp(th(nf+3))))/V;
    dF_dTheta(nf+4,1) = fu2*(-exp(th(nf+4))/2+log(t/exp(th(nf+3)))+exp(-th(nf+4))*(-exp(th(nf+4))/2+log(t/exp(th(nf+3)))).^2)/V ...
        - fu2/2/V;
    dF_dTheta(nf+5,1) = fu2*exp(-th(nf+5))/(1+exp(-th(nf+5)))/V ...
        - A*Rh./(t*sqrt(exp(th(nf+2))*pi)).*exp(-(log(t/exp(th(nf+1)))-exp(th(nf+2))/2).^2/exp(th(nf+2)))*exp(-th(nf+5))/(1+exp(-th(nf+5)))/V;
end

end
                       