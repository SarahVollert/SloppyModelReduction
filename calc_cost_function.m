function cost = calc_cost_function(theta, dataset, observed_output, sigma, model_output)
% This function calculates the cost function for the
% evaluation of the Hessian for the analysis of model sloppiness. 
[n_data,~] = size(dataset);

%Gaussian cost function
cost = n_data*log(sqrt(2*pi)*sigma)+(2*sigma^2)^-1*sum((model_output(theta,dataset) - observed_output').^2);

end