%% Main.m
% This is the main file which can run each component of the paper. From
% this document, you can run:
%(1) a Bayesian inference model calibration (SMC for uniform prior distributions) 
%(2) a maximum likelihood estimate model calibration 
%(3) a model sloppiness calculation
%(4) a comparison between two models 

% Please note that while this code is somewhat general, it was written
% specifically for application to the coral calcification model proposed by
% Galli & Solidoro (2018). If you are using this code for a different
% model, please use caution when modifying this code. 

%% (0) Select the model that you wish to analyse
%change the end of the function name to switch between calcification
%models:
% (a) _OriginalModel
% (b) _NoBATModel
% (c) _NoTranscellularModel
% (d) _NoBATNoTranscellularModel

define_uniform_priors = @define_uniform_priors_OriginalModel;
define_parameter_units = @define_parameter_units_OriginalModel;
model_output = @ model_output_OriginalModel;

%% (1) Bayesian inference for model calibration 

%import dataset
load('dataset.mat')
observed_output = dataset(:,6);

%define uniform priors
[Priors.prior_uppers, Priors.prior_lowers] = define_uniform_priors();
[n_parameters,~] = size(Priors.prior_uppers);
units = define_parameter_units(n_parameters); 

%%%%% Please select whether you want to generate the posterior sample:
%%%%% (TRUE) run the SMC algorithm to generate the sample, or 
%%%%% (FALSE) load the posterior sample shown in the manuscript
%%%%% Note: it can be computationally intensive to run the SMC algorithm.
Run_SMC_algorithm_option = false;

if Run_SMC_algorithm_option == true %Run the SMC algorithm
    
    %set other simulation conditions
    SMC_options.n_particles = 5000;                 %number of samples
    SMC_options.ESS = SMC_options.n_particles/2;    %minimum ESS of samples
    SMC_options.unique_particles_required = 0.95;   %portion of unique samples required
    SMC_options.MCMCmult = 0.99;                    %tuning parameter for speeding sampling
    SMC_options.max_gamma = 1;                      %target distribution temperature
    
    %run SMC algorithm (caution - computationally expensive for large samples)
    [prior_sample, posterior_sample] = SMC(dataset, observed_output, ...
        Priors, SMC_options,model_output);
    %NOTE: This is a generalised SMC algorithm, however, it was implemented
    %specifically for the coral calcification model. This SMC algorithm can
    %only work given uniform prior distributions (because of a hardcoded 
    %transform) and contains error checks each time the model is run, due to a 
    %pH convergence issue specific to the coral calcification model. 

else %load the results
    load('SMC_sample.mat');
    fprintf('Posterior sample from manuscript has been loaded\n')
end

%calculate the posterior calcification rates for this posterior sample
fprintf('Calculating posterior predictions, this may take a while for large samples...\n')
[posterior_prediction, posterior_loglikelihood] = calc_prediction_and_likelihood(posterior_sample, dataset, observed_output,model_output);
fprintf('...posterior predictions calculated\n')

%run Bayesian inference plots for the calcification model
Plots_PosteriorDensity(prior_sample, posterior_sample, Priors, units)
Plots_GoodnessOfFit(posterior_prediction, observed_output)


%% (2) MLE for model calibration

%%%%% Please select whether you want to generate the MLEs:
%%%%% (TRUE) run the MLE search algorithm, or 
%%%%% (FALSE) load the local maxima shown in the manuscript
%%%%% Note: it can be computationally intensive to search for many MLEs.
Run_MLE_algorithm_option = false;

if Run_MLE_algorithm_option == true %run MLE search algorithm

    %set up options for MLEs
    MLE_options.n_MLEs = 50;  %how many times we search for the MLE
    fprintf('Running the MLE search algorithm %d times\n',MLE_options.n_MLEs)
    
    %run the mle search through fmincon
    [MLE_sample,~] = MLE_SearchFromPosterior(MLE_options, posterior_sample, ...
                                        posterior_loglikelihood, Priors, dataset, observed_output,model_output);
    fprintf('MLE search algorithm has finished\n')
                                    
else %load the local maxima from manuscript
    load('MLE_sample.mat');
    fprintf('Local maxima from manuscript have been loaded\n')
end

%calculate the models predictions of the dataset for each MLE
[MLE_prediction, MLE_loglikelihood] = calc_prediction_and_likelihood(MLE_sample, dataset, observed_output, model_output);

%sort from best to worst fits
[MLE_loglikelihood, MLE_order] = sort(MLE_loglikelihood,'descend');
MLE_sample = MLE_sample(MLE_order,:);
MLE_prediction = MLE_prediction(MLE_order,:);

%number of MLEs to highlight for analysis
n_highlight = 5;
if n_highlight>size(MLE_sample,1)
    fprintf('Warning: Local maxima to be highlighted exceeds those available\n')
    fprintf('         Highlighting all local maxima instead\n')
    n_highlight = size(MLE_sample,1);
end

%run MLE in comparison to SMC plots for the calcification model
Plots_PosteriorDensityMLE(prior_sample, posterior_sample, Priors, units, MLE_sample,n_highlight)
Plots_GoodnessOfFitMLE(posterior_prediction, observed_output, MLE_prediction,n_highlight) 


%% (3) Analysis of Model Sloppiness 
%% (3a) A single sensitivity matrix at a time

%To calculate one sensitivity matrix of each type:

%initialise
delta = 1e-02;              %step size for finite differencing
mle = MLE_sample(1,:);      %selected MLE for evaluating Hessian

%calculate each sensitivity matrix
S_H = SensitivityMatrix_Hessian(dataset, mle, delta, observed_output, model_output);
S_P = SensitivityMatrix_posteriorCovariance(posterior_sample);
fprintf('Calculating LIS matrix, this may take a while for large posterior samples...\n')
S_L = SensitivityMatrix_likelihoodInformedSubspace(dataset, posterior_sample, prior_sample, delta, observed_output, model_output);
fprintf('...LIS matrix calculated\n')

%analyse the model sloppiness of one sensitivity matrix:
[eigenvalues_H, eigenparameters_H] = AnalyseSloppiness(S_H);
[eigenvalues_P, eigenparameters_P] = AnalyseSloppiness(S_P);
[eigenvalues_L, eigenparameters_L] = AnalyseSloppiness(S_L);

%compare each approach
Plots_eigenvalues_methodComparison(eigenvalues_H, eigenvalues_P, eigenvalues_L)

%% (3b) multiple sensitivity matrices

% Hessian evaluated at MLEs
%%%%% Please select whether you want to calculate the MLE Hessian:
%%%%% (TRUE) calculate the Hessian matrix at each of the MLEs in the MLE_sample, or 
%%%%% (FALSE) load the Hessian matrices shown in the manuscript
%%%%% Note: it can be computationally intensive to calculate for many MLEs.
Run_MLE_Hessian_calculation_option = false;
if Run_MLE_Hessian_calculation_option == true %calculate Hessian at MLEs in MLE_sample
    [S_H_multiple,eigenvalues_H_multiple,eigenparameters_H_multiple] = ...
        SensitivityMatrices_Hessian(dataset, MLE_sample, delta, observed_output, model_output);
else %load the Hessian matrices from the manuscript
    load('Hessian_MLE_analyses.mat')
end
        
% Posterior covariance method
%load multiple independent SMC samples used in the manuscript
load('SMC_sample_multiple.mat') %currently coral calcification posteriors
[S_P_multiple,eigenvalues_P_multiple,eigenparameters_P_multiple] = ...
    SensitivityMatrices_posteriorCovariance(posterior_sample_multiple);

% Likelihood informed subspace
%%%%% Please select whether you want to calculate the LIS matrices:
%%%%% (TRUE) calculate the LIS matrices for each posterior_sample, or 
%%%%% (FALSE) load the LIS matrices shown in the manuscript
%%%%% Note: it can be computationally intensive to calculate LIS matrices
%%%%% for large posterior samples. 
Run_LIS_calculation_option = false;
if Run_LIS_calculation_option == true %calculate Hessian at MLEs in MLE_sample
    [S_L_multiple,eigenvalues_L_multiple,eigenparameters_L_multiple] = ...
        SensitivityMatrices_likelihoodInformedSubspace(dataset, posterior_sample_multiple, ...
        prior_sample_multiple, delta, observed_output, model_output);
else %load the Hessian matrices from the manuscript
    load('LIS_analyses.mat')
end

% Plots
Plots_eigenvalues_methodComparison(eigenvalues_H_multiple, eigenvalues_P_multiple, eigenvalues_L_multiple)
%Note: this multiple analyses of model sloppiness can be added here where
%each column of the matrix is a different analysis. 

%% (4) Model Comparison                                                                              
%Note that this section was set up to compare the original and refined
%coral calcification models. As such, the results were saved under
%Original and Refined structures. This section is set up to import the
%results that were presented in the manuscript, however, these results
%could be generated by running the SMC algorithm for each of the models
%considered. 

%load the results for all models
load('SMC_sample_allModels.mat')
% Each model (Original, NoBAT, NoTranscellular, NoBATNoTranscellular) is 
% set up as a structure which contrains the information on priors (Priors),
% the posterior sample (posterior_sample), the prior sample (prior_sample)
% and the calculated posterior predictions (posterior_prediction). 

%Plot useful model comparisons
Plots_ModelComparison_PosteriorDensity(Original, NoBATNoTranscellular)
Plots_ModelComparison_GoodnessOfFit(Original, NoBAT, NoTranscellular, NoBATNoTranscellular, observed_output)

%compare eigenvalues of sloppy analysis between models
S_original = SensitivityMatrix_posteriorCovariance(Original.posterior_sample);
S_reduced = SensitivityMatrix_posteriorCovariance(NoBATNoTranscellular.posterior_sample);
[eigenvalues_original, eigenparameters_original] = AnalyseSloppiness(S_original);
[eigenvalues_reduced, eigenparameters_reduced] = AnalyseSloppiness(S_reduced);
Plots_eigenvalues_modelComparison(eigenvalues_original, eigenvalues_reduced)
