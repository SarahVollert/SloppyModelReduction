function [S_Ls, eigenvalues,eigenparameters] = SensitivityMatrices_likelihoodInformedSubspace(dataset, posterior_samples, prior_samples, delta, observed_output, model_output)
% This function runs multiple analyses of model sloppiness using the
% LIS approach, each evaluated for a different posterior sample.

%extract
[~,n_params, n_sets] = size(posterior_sample_multiple);

%initialise
eigenvalues = zeros(n_params,n_sets);
eigenparameters = zeros(n_params,n_params,n_sets);
S_Ls = zeros(n_params,n_params,n_sets);

%run sloppy analysis on each set of samples
for i=1:n_sets
    %calculate LIS matrix
    S_Ls(:,:,i) = SensitivityMatrix_likelihoodInformedSubspace(dataset, posterior_samples(:,:,i), prior_samples(:,:,i), delta, observed_output, model_output);
    
    %analyse sloppiness
    [eigenvalues(:,i), eigenparameters(:,:,i)] = AnalyseSloppiness(S_Ls(:,:,i));
end

end