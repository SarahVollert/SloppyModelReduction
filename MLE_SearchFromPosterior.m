%% MLE_SearchFromPosterior.m
function [MLE_sample,MLE_loglikelihood] = MLE_SearchFromPosterior(MLE_options, posterior_sample, ...
                                    posterior_loglikelihood, Priors, dataset, observed_output,model_output)

% This function searches for the MLE multiple times (as specified by 
% MLE_options.n_MLEs) using the fmincon search function. Each of these
% searches are initialised at the posterior samples with the highest 
% likelihoods. 


%sort posterior samples by likelihood
[~,sort_index] = sort(posterior_loglikelihood,'descend');
sorted_posterior_sample = posterior_sample(sort_index,:);
   
%maximise the log likelihood
MLEs = zeros(MLE_options.n_MLEs,size(posterior_sample,2));
LLs = zeros(MLE_options.n_MLEs,1);

%create a function to calculate the negative log-likelihood
calc_nll = @(parameter_vals)calc_nll_with_error_check(parameter_vals, dataset, observed_output,model_output);
%Note: due to errors within the pH algorithm for the calcification model, a
%try-catch error check is built into this function. 

%find each of the MLEs
parfor i=1:MLE_options.n_MLEs
    [MLE,nll] = fmincon(calc_nll,sorted_posterior_sample(i,:),[],[],[],[],Priors.prior_lowers,Priors.prior_uppers,[]);
    MLEs(i,:) = MLE;
    LLs(i) = -nll;
end

%keep only unique MLEs
[MLE_sample,index2] = unique(MLEs,'rows','stable');
MLE_loglikelihood = LLs(index2);

%sort the results from best to worst
[MLE_loglikelihood, index3] = sort(MLE_loglikelihood,'descend');
MLE_sample = MLE_sample(index3,:);
                      
end

