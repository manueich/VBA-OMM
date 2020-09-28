function [out] = VBA_OMM_G(dat,priors,const,opt)
%VBA_OMM Inversion of the oral minimal model (OMM) of glucose dynamics using 
% variational Bayesian analysis
%
% This function inverts the OMM in the followig form
%       dG/dt = - G*X - p1*(G-Gb) + (Ra+Rap)/V      G(0) = G0
%       dX/dt = - p2* (X - SI*(I-Ib))               X(0) = X0
%
% INPUT
%   - dat: Structure containg the data
%       - t: Time of sampling points in min. The first datapoint must be at t=0
%       - G: Glucose data at time points t in mmol/L
%       - I: Insulin data at time poinst t in arbitrary insulin unit (IU),
%       e.g. mU/L
%   - opt: Structure specifing inversion options
%       - GA_fun: Either 'RaPL' or 'RaLN', using either the piecewise
%       linear or log-normally based fucntions representing Ra
%       - tb: Time of breakpoints of RaPL in min. First breakpoint must be 
%       0 and last breakpoint must conincide with last datapoint.
%       - alpha: Exponential decay rate of RaPL after last breakpoint in
%       1/min 
%       - full_disp: 0 or 1 specifying whether full inversion figures are
%       displayed
%   - priors: Strcuture spefiying the priors. A coefficient of variation
%   (CV) of zero means the parameter is fixed and not updated during
%   inversion
%       - p1: 1x2 vector specifying median and CV of log-normal distribution
%       in 1/min and %
%       - p2: 1x2 vector specifying median and CV of log-normal distribution
%       in 1/min and %
%       - SI: 1x2 vector specifying median and CV of log-normal distribution
%       in 1/min per IU and %
%       - k: Mx2 matrix specifying median and CV of Ra function parameters.
%       CV in %.
%       For RaPL: M=length(tb)-2. Heights of Ra at breakpoints tb in 
%       mmol/kg/min, starting  with tb(2). The height at the penultimate 
%       breakpoint is calculated from the total AUC of RA and therefore 
%       not specified. 
%       For RaLN: M=5. Representing in order T1 in min, W1 no unit, 
%       T2 in min, W2 no unit and RH no unit. RH isrestricted to (0,1) 
%       and is not log-normally distributed.
%   - const: Structure specifying model constants
%       - dt: ODE integration step size in min
%       - A: AUC of Ra in mmol/kg
%       - V: Glucose distribution volume in L/kg
%       - X0: Initial condition of X, i.e. X(0) in 1/min
%       - G0: Initial condition of G, i.e. G(0) in mmol/L
%       - Gb: Basal level of G in mmol/L
%       - Ib: Basal level of I in IU
%       - measCV: measurement uncertainty CV of glucose assay in %
%       - Rap: Persiting absoption from a previous meal. An empty input 
%       means no persisting absoption. If present, Rap has to be a row 
%       vector coinciding with the integration time points on the
%       grid t(1):dt:t(end)
%
% OUTPUT out: Structure containing inversion results and input
%   - priors, options, data, const: see INPUT
%   - posterior: Structure containing posterior parameter distributions
%   analogous to prior structure
%   - Model_Output: Structure containing model output
%       - t: ODE integration time grid in min
%       - G, SigG: Model inferred glucose and uncertainty (SD) on t 
%       in mmol/L
%       - X, SigX: Model inferred state X and uncertainty (SD) on t 
%       in 1/min
%       - Ra, SigRa: Model inferred glucose appearance and uncertainty (SD) 
%       on t in in mmol/kg/min
%       - Rap: Persisting appearance for possible consecutive meal from
%       t(end) to 2*t(end)
%   - Performance: Structure containing model performance metrics
%       - FreeEnergy: Variational free energy, i.e. lower bound on log
%       model evidence
%       - R2: Coefficient of determination between data and model output
%       - AIC: Akaike Information Criterion
%       - BIC: Bayesian Information Criterion
%       - LL: Log Likelihood
%       - RMSE: Root mean squared error between data and model output
%       - wres: Weighted residuals between data and model output
%   - VB_Toolbox: Raw output of the VB toolbox. See description of function
%   VBA_NLStateSpaceModel() for details.
%
%   M. Eichenlaub 28/06/2020

%% Check Input
if nargin ~= 4
    disp('ERROR: Wrong number of input argumnets'); out = []; return
end

% Check dat
try
    t = dat.t(:);
    G = dat.G(:);
    I = dat.I(:);    
    if t(1)~=0 || length(t)~=length(G) || length(t)~=length(I)
        disp('ERROR: Dat structure is flawed'); out = []; return
    end
catch
    disp('ERROR: dat structure is flawed'); out = []; return
end

% Check const
try
    A = const.A;    
    V = const.V;
    dt = const.dt;
    if isempty(const.Rap)
        Rap = zeros(1,t(end)/dt+1);
    else
        Rap = const.Rap(:)';
    end
    X0 = const.X0;
    measCV = const.measCV;
    Gb = const.Gb;
    G0 = const.G0;
    Ib = const.Ib;
    
    if length(Rap)~=length(t(1):dt:t(end))
        disp('ERROR: Rap is the wrong length'); out = []; return
    end
catch
    disp('ERROR: const structure is flawed'); out = []; return
end

% Check Priors
try
    p_p1 = priors.p1;
    p_p2 = priors.p2;
    p_SI = priors.SI;
    p_k = priors.k;
catch
    disp('ERROR: Priors are flawed'); out = []; return
end

% Check opt
try
    full_disp = opt.full_disp;
    GA_fun = opt.GA_fun;
    
    switch GA_fun
        case 'RaPL'
            fname = @f_OMM_RaPL;
            tb = opt.tb;
            alpha = opt.alpha;
            p_k = priors.k;
            if tb(1)~=t(1); disp('ERROR: First breakpoint of RaPL must coincide with first datapoint, i.e t=0'); 
                out = []; return; end
            if tb(end)~=t(end); disp('ERROR: Last breakpoint of RaPL must coincide with last datapoint'); 
                out = []; return; end
            if size(p_k,1)+2~=length(tb) 
                disp('ERROR: Number of Breakpoints in RaPL does not match number of prior parameters');
                 out = []; return; end            
        case 'RaLN'
            fname = @f_OMM_RaLN;
            if size(p_k,1)~=5
                disp('ERROR: Wrong number of prior parameters for RaLN');
                 out = []; return; end
        otherwise
            disp('ERROR: GA function type flawed');
                 out = []; return; 
    end
catch
    disp('ERROR: opt structure flawed');  out = []; return
end

%% Construct model options

% Construct u
ti = t(1):dt:t(end);
Ii = interp1(t,I,ti);
u(1,:) = ti(1:end-1);
u(2,:) = Ii(1:end-1)-Ib;
u(3,:) = Rap(1:end-1);

% Interpolate Data and construct yout
Gdat = interp1(t,G,t(1):t(end));
if isfield(dat,'yout')
    yout = dat.yout;
    t = t(yout==0);
    G = G(yout==0);
    I = I(yout==0);
else
    yout = ones(1,length(Gdat)); yout(t+1) = 0;
end

% Construct inF
options.isYout = yout(2:end);
options.inF.V = V;
options.inF.Gb = Gb;
options.inF.G0 = G0;
options.inF.X0 = X0;
options.inF.dt = dt;
options.inF.A = A;

% Set options
options.decim = 1/(dt);
options.updateX0 = 0;
options.updateHP = 0;
options.TolFun = 1e-4;
options.GnTolFun = 1e-4;
options.backwardLag = 30;   
options.MaxIter = 100;  
options.MaxIterInit = 50;
options.microU = 1;
options.checkGrads = 0;
options.f_fname = fname;

% Construct priors
    % - System parameters
pr.muTheta(1) = log(p_p1(1));
pr.muTheta(2) = log(p_p2(1));
pr.muTheta(3) = log(p_SI(1));

if p_p1(2)==0; sth(1)=0; else sth(1) = log((p_p1(2)/100)^2+1); end
if p_p2(2)==0; sth(2)=0; else sth(2) = log((p_p2(2)/100)^2+1); end
if p_SI(2)==0; sth(3)=0; else sth(3) = log((p_SI(2)/100)^2+1); end

    % - Input parameters
switch GA_fun
    case 'RaPL'
        for i=1:size(p_k,1)
            pr.muTheta(3+i) = log(p_k(i,1));
            sth(3+i) = log((p_k(i,2)/100)^2+1);
        end
        options.inF.alpha = alpha;
        options.inF.tb = tb;
    case 'RaLN'
        for i=1:4
            pr.muTheta(3+i) = log(p_k(i,1));
            sth(3+i) = log((p_k(i,2)/100)^2+1);
        end
        [pr.muTheta(8),sth(8)] = f_logistic_mapping(p_k(5,1),p_k(5,2),1);
end

pr.muTheta = pr.muTheta';
pr.SigmaTheta = diag(sth);

    % - Initial Conditions
pr.muX0 = [G0;X0];
sX0(1) = (measCV*G0/100)^2; sX0(2)=0;
pr.SigmaX0 = diag(sX0);

    % - Measurement error
[pr.a_sigma,pr.b_sigma] = VBA_Create_NoisePrior(mean(Gdat)*measCV/100,mean(Gdat)*measCV/100*0.1);
for i=2:length(Gdat)
    pr.iQy{i-1,1} = 1/(Gdat(i)/mean(Gdat));
end

options.priors = pr;

% Construct dim
dim.n_theta = length(options.priors.muTheta);
dim.n_phi   = 0;
dim.n       = 2;

%% Simulate priors

if full_disp
    disp('Simulating model with PRIOR values ...');
    [X,Ra,SigX,SigRa] = f_simulate(options.priors,u,options,fname);

    col = lines;
    figure('Units','Normalized','Position',[0.05 0.5 0.4 0.3]),
    sgtitle('Model ouput and data from PRIOR values');
    subplot(121); hold on; box on;
    yyaxis left
    plot(t,G,'o'); ylabel('Glucose [mmol/L]');
    f_plotUncTimeSeries(ti,X(1,:),col(1,:),1.5,SigX(1,:));
    yyaxis right
    plot(t,I,'-o','Color',col(2,:)); ylabel('Insulin [IU]');
    xlabel('Time [min]');

    subplot(122); hold on; box on;
    f_plotUncTimeSeries(ti,Ra,col(1,:),1.5,SigRa);
    ylabel('Glucose Appearance [mmol/kg/min]')
    xlabel('Time [min]');
end


%% Model inversion

if full_disp
    options.verbose = 1;
    options.DisplayWin = 1;
else
    options.verbose = 0;
    options.DisplayWin = 0;
end

disp('Model Inversion ...');
[posterior,out_TB] = VBA_NLStateSpaceModel(Gdat(2:end),u,fname,@f_g,dim,options);
disp('DONE ...');

%% Wrapup results

out.priors = priors;
out.options = opt;
out.const = const;
out.data = dat;
out.VB_Toolbox.posterior = posterior;
out.VB_Toolbox.out = out_TB;

% Posterior
out.posterior.p1 = [exp(posterior.muTheta(1)) sqrt(exp(posterior.SigmaTheta(1,1))-1)*100];
out.posterior.p2 = [exp(posterior.muTheta(2)) sqrt(exp(posterior.SigmaTheta(2,2))-1)*100];
out.posterior.SI = [exp(posterior.muTheta(3)) sqrt(exp(posterior.SigmaTheta(3,3))-1)*100];
switch GA_fun
    case 'RaPL'
        for i=1:length(posterior.muTheta(4:end))
            out.posterior.k(i,:) = [exp(posterior.muTheta(3+i)) sqrt(exp(posterior.SigmaTheta(3+i,3+i))-1)*100];
        end
    case 'RaLN'
        for i=1:4
            out.posterior.k(i,:) = [exp(posterior.muTheta(3+i)) sqrt(exp(posterior.SigmaTheta(3+i,3+i))-1)*100];
        end
        [out.posterior.k(5,1),out.posterior.k(5,2)] = ...
            f_logistic_mapping(posterior.muTheta(8),posterior.SigmaTheta(8,8),2);
end

% Rap for possible next meal
u_temp = u;
u_temp(1,:) = u_temp(1,:)+ti(end);
[~,Ra,~,~] = f_simulate(posterior,u_temp,options,fname);
out.Model_Output.Rap = Ra - [u_temp(3,:) u_temp(3,end)];

% Model output
[X,Ra,SigX,SigRa] = f_simulate(posterior,u,options,fname);
out.Model_Output.t = u(1,:);
out.Model_Output.G = X(1,:);
out.Model_Output.SigG = SigX(1,:);
out.Model_Output.X = X(2,:);
out.Model_Output.SigX = SigX(2,:);
out.Model_Output.Ra = Ra;
out.Model_Output.SigRa = SigRa;

% Model performance
out.Performance.FreeEnergy = out.VB_Toolbox.out.F;
out.Performance.R2 = out.VB_Toolbox.out.fit.R2;
out.Performance.AIC = out.VB_Toolbox.out.fit.AIC;
out.Performance.LL = out.VB_Toolbox.out.fit.LL;
out.Performance.BIC = out.VB_Toolbox.out.fit.BIC;
    % = - Weighted residulas + RMSE
for i=1:length(t); idx = find(ti==t(i)); 
    wres(i)=(G(i)-out.Model_Output.G(idx))/(const.measCV/100*G(i)); 
    RMSE(i)=(G(i)-out.Model_Output.G(idx))^2;
end
out.Performance.wres = wres;
out.Performance.RMSE = sqrt(mean(RMSE));

%% Simulate Posterior

if full_disp
    disp('Simulating model with POSTERIOR values ...');
    
    col = lines;
    figure('Units','Normalized','Position',[0.05 0.1 0.4 0.3]),
    sgtitle('Model ouput and data from POSTERIOR values');
    subplot(121); hold on; box on;
    yyaxis left
    plot(t,G,'o'); ylabel('Glucose [mmol/L]');
    f_plotUncTimeSeries(ti,X(1,:),col(1,:),1.5,SigX(1,:));
    yyaxis right
    plot(t,I,'-o','Color',col(2,:)); ylabel('Insulin [IU]');
    xlabel('Time [min]');

    subplot(122); hold on; box on;
    f_plotUncTimeSeries(ti,Ra,col(1,:),1.5,SigRa);
    ylabel('Glucose Appearance [mmol/kg/min]')
    xlabel('Time [min]');
end


end

