function S_P = SensitivityMatrix_posteriorCovariance(posterior_sample)
% This function calculates the posterior covariance method sensitivity
% matrix, using a posterior sample. 

%find sensitivity matrix
S_P = inv(cov(log(posterior_sample(:,1:end-1))));

end