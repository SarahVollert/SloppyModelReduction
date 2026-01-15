function [S_Hs,eigenvalues,eigenparameters] = SensitivityMatrices_Hessian(dataset, ...
                                    MLE_sample, delta, observed_output, model_output)
% This function runs multiple analyses of model sloppiness using the
% Hessian matrix approach, each evaluated at different maximum likelihood
% estimates. 

%set up a sensitivity matrix for each MLE
n_params = size(MLE_sample,2);
n_mles = size(MLE_sample,1);
S_Hs = zeros(n_params-1,n_params-1,n_mles);
eigenvalues = zeros(n_params-1,n_mles); %each column is the normalised 
                                        %eigenvalues from a different MLE
eigenparameters = zeros(n_params-1,n_params-1,n_mles); %the rows of each 
                                        % matrix are the normalised 
                                        % eigenparameters such that each 
                                        % matrix is an analysis of model
                                        % sloppiness for an MLE.

% run the analysis for each MLE                                     
parfor i=1:size(MLE_sample,1)
    %get the sensitivity matrix
    S_Hs(:,:,i) = SensitivityMatrix_Hessian(dataset, MLE_sample(i,:), delta, observed_output, model_output);

    %get the eigenvalues and eigenparameters
    [eigenvalues(:,i), eigenparameters(:,:,i)] = AnalyseSloppiness(S_Hs(:,:,i));
end

end

