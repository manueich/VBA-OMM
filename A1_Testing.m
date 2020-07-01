clc
clear
close all

%%

load('dat_ex');

dat.t = dat_ex(:,1);
dat.G = dat_ex(:,2);
dat.I = dat_ex(:,3);

const.A = 6;
const.V = 0.145;
const.dt = 0.1;
const.Rap = zeros(1,length(dat.t(1):const.dt:dat.t(end)));
const.X0 = 0;
const.measCV = 2;
const.Gb = dat.G(1);
const.G0 = dat.G(1);
const.Ib = dat.I(1);

opt.GA_fun = 'RaPL';
opt.tb = [0 10 30 60 90 120 180 300];
opt.alpha = 0.017;
opt.full_disp = 1;

priors.p1 = [0.025 25];
priors.p2 = [0.012 40];
priors.SI = [7.1e-4 100];
priors.k = [3.2e-3*const.A 50;     % k1 at tb(2)
            7.3e-3*const.A 50;     % k2 at tb(3)
            5.4e-3*const.A 50;     % k3 at tb(4)
            5.1e-3*const.A 50;     % k4 at tb(5)
            3.7e-3*const.A 50;     % k5 at tb(6)
            1.8e-3*const.A 50];    % k7 at tb(8)
        
% priors.k = [30 30;      % T1
%             0.5 30;     % W1
%             100 30;     % T2
%             0.5 30;     % W2
%             0.7 30];    % RH     
        
out = VBA_OMM_G(dat,priors,const,opt);

