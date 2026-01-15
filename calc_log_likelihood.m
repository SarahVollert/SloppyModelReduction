function loglikelihood = calc_log_likelihood(theta, dataset, observed_output,model_output)
% This function calculates the Gaussian log likelihood of a model, for a
% given input dataset (dataset), output dataset (observed_output) and
% parameter values (theta). %Note that theta must include sigma as the 
% final parameter for this Gaussian log-likelihood. \

[n_data,~] = size(dataset);

%Gaussian log-likelihood
loglikelihood = -n_data*log(theta(end))-sum((model_output(theta,dataset)-observed_output').^2/(2*theta(end)^2));
end