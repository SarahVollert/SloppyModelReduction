function [prediction, loglikelihood] = calc_prediction_and_likelihood(parameter_sample, dataset, observed_output,model_output)
% This function returns the model predictions for each of the data points,
% and the log-likelihood, for each parameter sample passed to it. The model
% is passed as model_output. 


%get number of datapoints
n_data = length(dataset);

%set up calcification and negative log likelihood storage
prediction = zeros(size(parameter_sample,1),n_data);
loglikelihood = zeros(size(parameter_sample,1),1);

%calculate calcification rate and negative log likelihood for each posterior sample
parfor k=1:size(parameter_sample,1)
    prediction(k,:) = model_output(parameter_sample(k,:),dataset);
    loglikelihood(k) = calc_log_likelihood(parameter_sample(k,:),dataset, observed_output,model_output);
end

end