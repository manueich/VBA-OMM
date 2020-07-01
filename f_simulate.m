function [X,Ra,SigX,SigRa] = f_simulate(p,u,options,fname)
%F_SIMULATE Summary of this function goes here
%   Detailed explanation goes here

% Preallocate variables
SigX = cell(length(u),1);
dfdx0 = cell(length(u),1);
dfdth = cell(length(u),1);

% Set initial conditions
SigRa{1} = 0;
dfdx0{1} = eye(length(p.muX0));
dfdth{1} = zeros(length(p.muTheta),length(p.muX0));
X(:,1) = p.muX0;
if options.updateX0 
    SigX{1} = p.SigmaX0;
else
    SigX{1} = zeros(length(p.muX0));
end

% Get Covariace matrices 
    % - GA function parameters
SigThRa = p.SigmaTheta(4:end,4:end);
    % - Parameters + initial conditions
try     % Posterior
    SigThX0 = p.out.suffStat.ODE_posterior.SigmaPhi;        
catch   % Prior
    if options.updateX0 
        SigThX0 = diag([diag(p.SigmaTheta); diag(p.SigmaX0)]);
    else
        SigThX0 = p.SigmaTheta;
    end    
end
             
for i=1:size(u,2)-1
    [X(:,i+1),Jt,Ht] = feval(fname,X(:,i),p.muTheta,u(:,i),options.inF);
    
    if ~isempty(Jt)
        dfdth{i+1} = dfdth{i}*Jt + Ht;
        dfdx0{i+1} = dfdx0{i}*Jt;

        if options.updateX0 
            D = [dfdth{i+1}; dfdx0{i+1}];
        else
            D = [dfdth{i+1}];
        end
        
        % Calculate SigX
        SigX{i+1} = D'*SigThX0*D;
        
        % Calculate SigRa
        Ht = Ht(4:end,1)/options.inF.dt*options.inF.V;
        SigRa{i+1} = Ht'*SigThRa*Ht;

    else
        SigX{i+1} = SigX{i};
    end   
    
    Ra(i) = feval(fname,u(1,i),p.muTheta,options.inF) + u(3,i);
end

SigX = sqrt(VBA_getVar(SigX,length(SigX)));
SigRa = sqrt(VBA_getVar(SigRa,length(SigRa)));

Ra(end+1) = feval(fname,u(1,end),p.muTheta,options.inF) + u(3,end);

end

