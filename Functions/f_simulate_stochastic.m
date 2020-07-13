function [G,X,Ra] = f_simulate_stochastic(p,u,options,n,x0)
%F_SIMULATE_STOCAHSTIC Summary of this function goes here
%   Detailed explanation goes here

theta = p.muTheta;
sigTh = p.SigmaTheta;
fname = options.f_fname;

TH = mvnrnd(theta,sigTh,n);

if x0
    X0(:,1) = randn(n,1)*0.1 + p.muX0(1);
else
    X0(:,1) = ones(n,1)*p.muX0(1);
end

X0(:,2) = ones(n,1)*p.muX0(2);

fprintf(1,'Stochastic Simulation Progress: %3d%%\n',0);

parfor i=1:n
    th = TH(i,:)';
    [Z,Ra(i,:)] = solve_ode(fname,th,u,options,X0(i,:));
    G(i,:) = Z(1,:);
    X(i,:) = Z(2,:);
    fprintf(1,'\b\b\b\b%3.0f%%',i/n*100);
end
fprintf('\n');

end


function [X,Ra] = solve_ode(fname,th,u,options,X0)

X(:,1) = X0;

for i=1:size(u,2)
    [X(:,i+1)] = feval(fname,X(:,i),th,u(:,i),options.inF);
    Ra(i) = feval(fname,u(1,i),th,options.inF) + u(3,i);
end
Ra(end+1) = feval(fname,u(1,end),th,options.inF) + u(3,end);

end