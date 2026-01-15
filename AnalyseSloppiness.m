function [normalised_eigenvalues, eigenparameters] = AnalyseSloppiness(S)
% This function does an analysis of model sloppiness for a given
% sensitivity matrix (S). Here we do a spectral decomposition, and output 
% the normalised eigenvalues, the normalised eigenvector values for each
% eigenparameter.

%initialise
n_params = length(S);

%eigendecomposition
[eigenvectors,eigenvalues_matrix,~] = svd(S); 
eigenvalues = diag(eigenvalues_matrix);

%Sort by eigenvalues (largest to smallest)
[sorted_eigenvalues,index_eigenvalues] = sortrows(eigenvalues,'descend');
sorted_eigenvectors = eigenvectors(:,index_eigenvalues); 

%normalise eigenvalues by the leading eigenvalue (now each eigenvalue is
%between 0 and 1, where 1 indicates the highest relative influence of an 
%eigenparameter)
normalised_eigenvalues = sorted_eigenvalues/(max(eigenvalues));

%transform each eigenvector to be between -1 and 1 (now the contribution of
%each eigenvector value is a clear ratio of the parameters)
normalised_eigenvectors = zeros(n_params);
for i=1:n_params
    if (max(abs(sorted_eigenvectors(:,i))) == max(sorted_eigenvectors(:,i))) %if positive
        normalised_eigenvectors(:,i) = sorted_eigenvectors(:,i)/max(abs(sorted_eigenvectors(:,i)));
    else %negative leading value
        normalised_eigenvectors(:,i) = -sorted_eigenvectors(:,i)/max(abs(sorted_eigenvectors(:,i)));
    end
end

%transpose such that eigenparameters are the rows
eigenparameters = transpose(normalised_eigenvectors);

end