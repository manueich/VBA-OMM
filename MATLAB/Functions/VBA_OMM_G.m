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
%       - displayWin: 0 or 1 specifying whether inversion results are
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
%       For RaPL: M=6. Heights of Ra at breakpoints tb in 
%       mmol/kg/min, starting  with tb(2). The height at the penultimate 
%       breakpoint is calculated from the total AUC of RA and therefore 
%       not specified. NB: the toolbox currrently only supports RaPL with 
%       exactly 8 breakpoints. The times of these breakpoints can however 
%       be chosen freely. Please contact the developers if a different 
%       number of breakpoints is required.
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
%   analogous to prior structure with additional matrix of posterior
%   correlation coefficients
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
    displayWin = opt.displayWin;
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
u(1,:) = ti(1:end);
u(2,:) = Ii(1:end)-Ib;
u(3,:) = Rap(1:end);

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
options.verbose = 0;
options.DisplayWin = 0;

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

%% Model inversion

if displayWin; disp('Model Inversion ...'); end

[posterior,out_TB] = VBA_NLStateSpaceModel(Gdat(2:end),u,fname,@f_g,dim,options);
if isempty(posterior)
    disp('Model Inversion FAILED');
    out = []; 
    return
end

if displayWin; disp('DONE ...'); end

%% Wrapup results

out_TB.u = u;
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
        out.posterior.Correlation = f_get_correlation(posterior.muTheta,...
            posterior.SigmaTheta,1);
    case 'RaLN'
        for i=1:4
            out.posterior.k(i,:) = [exp(posterior.muTheta(3+i)) sqrt(exp(posterior.SigmaTheta(3+i,3+i))-1)*100];
        end
        [out.posterior.k(5,1),out.posterior.k(5,2)] = ...
            f_logistic_mapping(posterior.muTheta(8),posterior.SigmaTheta(8,8),2);
        out.posterior.Correlation = f_get_correlation(posterior.muTheta,...
            posterior.SigmaTheta,2);
end

% Rap for possible next meal
u_temp = u;
u_temp(1,:) = u_temp(1,:)+ti(end);
[~,Ra,~,~] = f_simulate(posterior,u_temp,options,fname);
out.Model_Output.Rap = Ra - u_temp(3,:);

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

%% Plot Figure with Results

if displayWin
    disp('Creating Results Figure ...')
    col = lines;
    % Set up Figure Window
    pos0 = get(0,'ScreenSize');
    pos = [0.51*pos0(3),0.3*pos0(4),0.45*pos0(3),0.5*pos0(4)];
    fig = figure('Name','OMM Results','Position',pos,'menubar','none');
    tg = uitabgroup(fig);
    tb5 = uitab(tg,'Title','Summary');
    tb1 = uitab(tg,'Title','Data');
    tb2 = uitab(tg,'Title','Priors');
    tb3 = uitab(tg,'Title','Posterior');
    tb4 = uitab(tg,'Title','Parameters');
    
    % Summary Tab ------------------------------------------------
    
    str{1} = ['Inversion results of the oral minimal model using ' GA_fun];
    str{end+1} = '';    
    str{end+1} = sprintf('Elaspsed time %.1f seconds',out.VB_Toolbox.out.dt);
    str{end+1} = '';
    str{end+1} = 'Model fit criteria';
    str{end+1} = sprintf('   - R2:  %2.1f %%',out.Performance.R2*100);
    str{end+1} = sprintf('   - RMSE:  %.2f mmol/L',out.Performance.RMSE);
    str{end+1} = '';
    str{end+1} = 'Parameter estimation results [median +/- CV in %]';
    str{end+1} = sprintf('   - Glucose effectiveness p1 [1E-3 1/min]:  %2.2f +/- %2.1f %%',...
        out.posterior.p1(1)*1e3,out.posterior.p1(2));
    str{end+1} = sprintf('   - Decay parameter p2 [1E-3 1/min]:  %2.2f +/- %2.1f %%',...
        out.posterior.p2(1)*1e3,out.posterior.p2(2));
    str{end+1} = sprintf('   - Insulin sensitivity SI [1E-4 1/min per IU]:  %2.2f +/- %2.1f %%',...
        out.posterior.SI(1)*1e4,out.posterior.SI(2));
    
    uicontrol('Parent',tb5,...
        'Style','text',...
        'units','normalized',...
        'position',[0.1,0,0.8,0.85],...
        'fontsize',11,...
        'HorizontalAlignment','left',...
        'string',str)        
   
    % Data Tab ------------------------------------------------
    axes('Parent',tb1)
    sgtitle('Provided DATA');
    subplot(2,2,1), hold on; box on; 
    title('Glucose Data')
    errorbar(t,G,G*const.measCV/100,'o--','Color',col(1,:),'MarkerFaceColor',col(1,:),...
        'MarkerSize',3);
    line([t(1) t(end)],[Gb Gb],'LineStyle','-.','Color',col(1,:))
    ylim([floor(min(G)-1) ceil(max(G)+1)])
    ylabel('Glucose [mmol/L]');
    xlabel('Time [min]');    
    
    subplot(2,2,2), hold on; box on; 
    title('Insulin Data')
    plot(t,I,'o--','Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerSize',3);
    line([t(1) t(end)],[Ib Ib],'LineStyle','-.','Color',col(2,:))
    ylim([0 floor(max(I)+5)])
    ylabel('Insulin [IU]');
    xlabel('Time [min]');    
        
    subplot(2,2,3), hold on; box on; 
    title('Persisting Absorption Rap')
    plot(ti,Rap,'-','Color',col(1,:)); 
    ylabel('GA [mmol/kg/min]');
    xlabel('Time [min]');
    
    % Prior Tab ------------------------------------------------
    axes('Parent',tb2)
    sgtitle('Model ouput and data from PRIOR values');
    [X,Ra,SigX,SigRa] = f_simulate(options.priors,u,options,fname);
    
    subplot(2,2,1), hold on; box on; 
    title('Glucose G(t)')
    plot(t,G,'o','Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerSize',4);
    line([t(1) t(end)],[Gb Gb],'LineStyle','-.','Color',col(1,:))
    f_plotUncTimeSeries(ti,X(1,:),col(1,:),1.5,SigX(1,:));
    ylabel('Glucose [mmol/L]');
    xlabel('Time [min]');
    
    subplot(2,2,2), hold on; box on; 
    title('Active Insulin X(t)')
    f_plotUncTimeSeries(ti,X(2,:),col(2,:),1.5,SigX(2,:));
    ylabel('X [1/min]');
    xlabel('Time [min]');
    
    subplot(2,2,3), hold on; box on; 
    title('Glucose Absorption Ra(t)')
    f_plotUncTimeSeries(ti,Ra,col(1,:),1.5,SigRa);
    ylabel('GA [mmol/kg/min]')
    xlabel('Time [min]');
    
    % Posterior Tab ------------------------------------------------
    axes('Parent',tb3)
    sgtitle('Model ouput and data from POSTERIOR values');
    [X,Ra,SigX,SigRa] = f_simulate(posterior,u,options,fname);
    
    subplot(2,2,1), hold on; box on; 
    title('Glucose G(t)')
    plot(t,G,'o','Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerSize',4);
    line([t(1) t(end)],[Gb Gb],'LineStyle','-.','Color',col(1,:))
    f_plotUncTimeSeries(ti,X(1,:),col(1,:),1.5,SigX(1,:));
    ylabel('Glucose [mmol/L]');
    xlabel('Time [min]');
    
    subplot(2,2,2), hold on; box on; 
    title('Active Insulin X(t)')
    f_plotUncTimeSeries(ti,X(2,:),col(2,:),1.5,SigX(2,:));
    ylabel('X [1/min]');
    xlabel('Time [min]');
    
    subplot(2,2,3), hold on; box on; 
    title('Glucose Absorption Ra(t)')
    f_plotUncTimeSeries(ti,Ra,col(1,:),1.5,SigRa);
    ylabel('GA [mmol/kg/min]')
    xlabel('Time [min]');
    
    subplot(2,2,4), hold on; box on; 
    title('Weighted Residuals')
    plot(t,wres,'-o','Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerSize',4)
    line([t(1) t(end)],[1 1],'LineStyle','-.','Color','k')
    line([t(1) t(end)],[-1 -1],'LineStyle','-.','Color','k')
    line([t(1) t(end)],[0 0],'LineStyle','-','Color','k')
    xlabel('Time [min]');
    
    % Parameters Tab ------------------------------------------------
    axes('Parent',tb4)
    sgtitle('Posterior parameter distributions');
    subplot(2,2,1), hold on; box on; 
    title('Glucose effectiveness p1')
    p(1) = plot_dist(options.priors.muTheta(1),sqrt(options.priors.SigmaTheta(1,1)),...
        priors.p1(1),priors.p1(2),col(1,:));
    p(2) = plot_dist(posterior.muTheta(1),sqrt(posterior.SigmaTheta(1,1)),...
        out.posterior.p1(1),out.posterior.p1(2),col(2,:));
    legend(p,{'Prior Distr.','Posterior Distr.'},'Location','best')
    xlabel('p_{1} [1/min]'); ylabel('Probability Density')
    
    subplot(2,2,2), hold on; box on; 
    title('Decay parameter p2')
    p(1) = plot_dist(options.priors.muTheta(2),sqrt(options.priors.SigmaTheta(2,2)),...
        priors.p2(1),priors.p2(2),col(1,:));
    p(2) = plot_dist(posterior.muTheta(2),sqrt(posterior.SigmaTheta(2,2)),...
        out.posterior.p2(1),out.posterior.p2(2),col(2,:));
    legend(p,{'Prior Distr.','Posterior Distr.'},'Location','best')
    xlabel('p_{2} [1/min]'); ylabel('Probability Density')
    
    subplot(2,2,3), hold on; box on; 
    title('Insulin Sensitivity SI')
    p(1) = plot_dist(options.priors.muTheta(3),sqrt(options.priors.SigmaTheta(3,3)),...
        priors.SI(1),priors.SI(2),col(1,:));
    p(2) = plot_dist(posterior.muTheta(3),sqrt(posterior.SigmaTheta(3,3)),...
        out.posterior.SI(1),out.posterior.SI(2),col(2,:));
    legend(p,{'Prior Distr.','Posterior Distr.'},'Location','best')
    xlabel('S_{I} [1/min per IU]'); ylabel('Probability Density')
    
    subplot(2,2,4), hold on; box on; 
    title('GA function parameters')
    p(1)=plot_errorbars([1:size(priors.k,1)]-0.1,priors.k,col(1,:),GA_fun);
    p(2)=plot_errorbars([1:size(out.posterior.k,1)]+0.1,...
        out.posterior.k,col(2,:),GA_fun);
    legend(p,{'Prior Distr.','Posterior Distr.'},'Location','best') 
    
    switch GA_fun
        case 'RaPL'
            ylabel('GA [mmol/kg/min]');
            for i=1:size(priors.k,1); labels{i} = ['k' num2str(i)]; end
            labels{i} = ['k' num2str(i+1)];
            set(gca,'XTick',1:1:size(priors.k,1));
            set(gca,'XTickLabels',labels);
        case 'RaLN'
            set(gca,'XTick',1:1:size(priors.k,1));
            set(gca,'XTickLabels',...
                {'     T_1\newline[10^{2}min]','W_1','     T_2\newline[10^{2}min]','W_2','R_H'});
    end

end

end

function [p] = plot_dist(mu,sig,med,CV,col)   
    xp = linspace(0,med+4*(med*CV/100),1000);
    p = plot(xp,lognpdf(xp,mu,sig),'LineWidth',1.1,'Color',col);
    line([med med],[0 lognpdf(med,mu,sig)],'LineStyle','--','Color',col)   
end
   

function [p] = plot_errorbars(pos,k,col,GA_fun)   
    
    switch GA_fun
        case 'RaPL'
            for i=1:length(pos)
                mu = k(i,1);
                sig = exp(sqrt(log((k(i,2)/100)^2+1)));
                p = errorbar(pos(i),mu,mu-mu/sig,mu*sig-mu,'o','Color',col,...
                    'MarkerFaceColor',col,'MarkerSize',5,'LineWidth',1.2);
            end
        case 'RaLN'
            for i=1:length(pos)-1
                sc = 1; if i==1 || i==3; sc = 0.01; end
                mu = k(i,1)*sc;
                sig = exp(sqrt(log((k(i,2)/100)^2+1)));
                p = errorbar(pos(i),mu,mu-mu/sig,mu*sig-mu,'o','Color',col,...
                    'MarkerFaceColor',col,'MarkerSize',5,'LineWidth',1.2);
            end
            [m_o,s_o] = f_logistic_mapping(k(end,1),k(end,2),1);
            lo = f_sig(m_o-sqrt(s_o));
            hi = f_sig(m_o+sqrt(s_o));  
            mu = k(end,1);
            errorbar(pos(end),mu,mu-lo,hi-mu,'o','Color',col,...
                'MarkerFaceColor',col,'MarkerSize',5,'LineWidth',1.2);      
            
    end 
end


function s=f_sig(x)
s = 1./(1+exp(-x));
end
    


