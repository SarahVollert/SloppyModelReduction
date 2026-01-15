function [S_Ps,eigenvalues,eigenparameters] = SensitivityMatrices_posteriorCovariance(posterior_sample)
% This function runs multiple analyses of model sloppiness using the
% posterior covariance method, each evaluated at an indepenently generated
% (using SMC) posterior sample.

%set up a sensitivity matrix for each posterior sample
n_params = size(posterior_sample,2);
n_samples = size(posterior_sample,3);
S_Ps = zeros(n_params-1,n_params-1,n_samples);
eigenvalues = zeros(n_params-1,n_samples); %each column is the normalised 
                                        %eigenvalues from a different
                                        %posterior sample
eigenparameters = zeros(n_params-1,n_params-1,n_samples); %the rows of each 
                                        % matrix are the normalised 
                                        % eigenparameters such that each 
                                        % matrix is an analysis of model
                                        % sloppiness for an posterior
                                        % sample.

% run the analysis for each sample                                     
for i=1:n_samples
    %get the sensitivity matrix
    S_Ps(:,:,i) = SensitivityMatrix_posteriorCovariance(posterior_sample(:,:,i));

    %get the eigenvalues and eigenparameters
    [eigenvalues(:,i), eigenparameters(:,:,i)] = AnalyseSloppiness(S_Ps(:,:,i));
end

end