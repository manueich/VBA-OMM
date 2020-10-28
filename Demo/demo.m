%% Demo on VBA-OMM
%   All folders and subfolders of the VBA toolbox as well as the Functions
% folder must be in the MATLAB path

clc
clear
close all

% Load demo glucose and insulin data from Dalla Man et. (2002): "The oral 
% glucose minimal model: estimation of insulin sensitivity from a meal test"
% 1st column: time of measurement points
% 2nd column: glucose data
% 3rd column: insulin data

load('demo_dat');

% Construt the data struture
dat.t = demo_dat(:,1);                              % Time
dat.G = demo_dat(:,2);                              % Glucose 
dat.I = demo_dat(:,3);                              % Insulin

% Construct the constants structure
const.A = 6;                                    % AUC of glucose appearance function in mmol/kg
const.V = 0.145;                                % volume of glucose distribution in L/kg
const.dt = 0.1;                                 % integration time step in min
const.Rap = [];                                 % persisting absoption in mmol/kg/min. An empty array means no absorption
% const.Rap = 0.01*...                            % Example of exponetial decay from height 0.01 and rate 0.017     
%      exp(-0.017*[dat.t(1):const.dt:dat.t(end)]);     
const.X0 = 0;                                   % initial condition of state X
const.measCV = 2;                               % CV glucose assay in %
const.Gb = dat.G(1);                            % basal level of glucose
const.G0 = dat.G(1);                            % initial condition of glucose
const.Ib = dat.I(1);                            % basal level of insulin

% Construct inversion options
opt.GA_fun = 'RaPL';                        % glucose appearance function either 'RaPL' or 'RaLN'
opt.tb = [0 10 30 60 90 120 180 300];       % breakpoints of RaPL'
opt.alpha = 0.017;                          % decay rate of RaPL'
opt.displayWin = 1;                         % display Figure with inversion results

% Priors
    % - System Parameters [median CV]
priors.p1 = [0.025 25];
priors.p2 = [0.012 40];
priors.SI = [7.1e-4 100];

    % - Input function Parameters
switch opt.GA_fun
    case 'RaPL'
        priors.k = [3.2e-3*const.A 50;     % k1 at tb(2)
                    7.3e-3*const.A 50;     % k2 at tb(3)
                    5.4e-3*const.A 50;     % k3 at tb(4)
                    5.1e-3*const.A 50;     % k4 at tb(5)
                    3.7e-3*const.A 50;     % k5 at tb(6)
                    1.8e-3*const.A 50];    % k7 at tb(8)
    case 'RaLN'
        priors.k = [30 30;      % T1
                    0.5 30;     % W1
                    100 30;     % T2
                    0.5 30;     % W2
                    0.7 30];    % RH 
end

% Call inversion routine
out = VBA_OMM_G(dat,priors,const,opt);

