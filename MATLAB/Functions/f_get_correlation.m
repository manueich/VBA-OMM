function [C] = f_get_correlation(mu,Sigma,CASE)
%f_get_correlation Calculates the posterior correlation matrix
%   Takes the normally distributed posterior parameter distribution and
%   calculates the posterior corraltion matrix based on their log-normally
%   distributed substitutes
%   IN:
%       - mu: Posterior vector of mean (normally distributed)
%       - Sigma: Posterior covariance matrix (normally distributed)
%       - CASE: 1: all parameters have been transformed to log-normal (RaPL)       
%               2: all parameters have been transformed to log-normal
%               except for the last one which has been logistically trasnformed (RaLN) 
%   OUT:
%       - C: Posterior correlation matrix

switch CASE
    case 1      % all Lognormal
       for i=1:size(Sigma,1)
            for j=1:size(Sigma,2)
                % Formula from Zhang et al. (2015): Inferences on correlation 
                % coefficients of bivariate log-normal distributions, IN:
                % Journal of Applied statistics, 42(3)
                C(i,j) = (exp(Sigma(i,j))-1)/(sqrt((exp(Sigma(i,i))-1)*(exp(Sigma(j,j))-1)));
            end
        end
    case 2      % Last logictic
        % Runs a Monte-Carlo Simulation with 1e6 samples from the
        % multivariate normal distribution and calculates the sample
        % correlation after transformation
        N=1e6;
        try
            R = mvnrnd(mu,Sigma,N);
            R = [exp(R(:,1:end-1)) 1./(1+exp(-R(:,end)))];
            C = corrcoef(R);
        catch
            for i=1:size(Sigma,1)
                for j=1:size(Sigma,2)
                    C(i,j) = (exp(Sigma(i,j))-1)/(sqrt((exp(Sigma(i,i))-1)*(exp(Sigma(j,j))-1)));
                end
            end
        end
end

end

