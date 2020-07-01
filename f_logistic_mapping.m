function [m_o,s_o] = f_logistic_mapping(m_i,s_i,mod)
%F_LOGISTIC_MAPPING Converts median and CV for logistic mapping
% Mod 1: Logistic to Normal
%   m_i: Median of Logistic
%   s_i: CV of Logistic
%   
%   m_o: Mean of Normal
%   s_o: Var of Normal

% Mod 2: Normal to Logistic
%   m_i: Mean of Normal
%   s_i: Var of Normal
%   
%   m_o: Median of Logistic
%   s_o: CV of Logistic

if mod==1       % Logistic to Normal
    [m_o,s_o] = Log2Norm(m_i,s_i);    
else            % Normal to logistic
    [m_o,s_o] = Norm2Log(m_i,s_i);
end

end

function [m_o,s_o] = Log2Norm(m_i,s_i)

m_o = -log(1/m_i-1);
sig = 0.01:0.01:10;
for i=1:length(sig)
    [~,CV(i)] = Norm2Log(m_o,sig(i));
end

[~,idx] = min(abs(s_i-CV));
s_o = sig(idx);

end

function [m_o,s_o] = Norm2Log(m_i,s_i)
a = 3/pi^2;

m_o = f_logistic(m_i);
mea = f_logistic(m_i/sqrt(1+a*sqrt(s_i)));
vari = mea*(1-mea)*(1-1/sqrt(1+a*sqrt(s_i))); 
s_o = sqrt(vari)/mea*100;
    
end

function s = f_logistic(x)
s = 1./(1+exp(-x));
end

