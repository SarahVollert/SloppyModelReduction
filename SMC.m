function [prior_sample, posterior_sample] = SMC(dataset, observed_output, ...
    Priors, SMC_options,model_output)
%This function runs a Sequential Monte Carlo algorithm to generate a
%posterior sample.

%NOTE: This is a generalised SMC algorithm, however, it was implemented
%specifically for the coral calcification model. This SMC algorithm can
%only work given uniform prior distributions (because of a hardcoded 
%transform) and contains error checks each time the model is run, due to a 
%pH convergence issue specific to the coral calcification model. 

%INITIALISATION
n_particles = SMC_options.n_particles;
ESS = SMC_options.ESS;
unique_particles_required = SMC_options.unique_particles_required;
MCMCmult = SMC_options.MCMCmult; 
prior_lowers = Priors.prior_lowers;
prior_uppers = Priors.prior_uppers;
n_params = length(prior_lowers);
max_gamma = SMC_options.max_gamma;




%GENERATE PRIOR SAMPLE
theta = draw_Priors(n_particles,Priors); %particle values
particle_w = ones(n_particles,1)/n_particles; %particle weightings
log_likelihood = zeros(n_particles,1); %log likelihood of each particle

%remove particles which cause errors from the prior sample
parfor j=1:n_particles %for each particle
    % set a flag to see if an error has occurred
    flag = false; %F = no error, T = error
    
    while flag == false %while there are no errors
        flag = true;
        
        try %to run the model
            model_output(theta(j,:),dataset);
            
        catch %any errors running the model
            flag = false; %there has been an error - flag it
            
            %replace problematic particle
            new_theta = draw_Priors(1,Priors);
            theta(j,:) =  new_theta(1,:);
        end
    end %continue once there are no errors
    
    %calculate particle likelihoods
    log_likelihood(j) = calc_log_likelihood(theta(j,:), dataset, observed_output,model_output);
end
prior_sample = theta;




%RUN THE SEQUENCE TO POSTERIOR
%while the temperature is < desired temperature
gamma_t = 0; %intial temperature
while gamma_t < max_gamma
    
    %Select next temperature in sequence
    ess_posterior = calc_ess(max_gamma,gamma_t,particle_w,log_likelihood); %check ess for final posterior
    if ess_posterior > ESS %if ESS at posterior is acceptable
        newgamma = max_gamma; %move to posterior
    else
        %find the temperature at which the ESS becomes unacceptable
        newgamma = fzero(@(newgamma) calc_ess(newgamma,gamma_t,particle_w,log_likelihood)-ESS,[gamma_t max_gamma]);
    end
    
    fprintf('*** likelihood annealing temperature is %f***\n',newgamma);
    



    %REWEIGHT
    %reweight
    log_particle_w = log(particle_w) +(newgamma - gamma_t)*log_likelihood;
    
    %numerically stabilise
    log_particle_w = log_particle_w - max(log_particle_w);
    particle_w = exp(log_particle_w);
    
    %normalise weights
    particle_w = particle_w/sum(particle_w);
    



    %RESAMPLING
    %apply transform (uniform priors only)
    theta_tilde = zeros(n_particles,n_params);
    for k=1:n_params
        theta_tilde(:,k)=log((theta(:,k)-prior_lowers(k))./(prior_uppers(k)-theta(:,k)));
    end
    
    %calc cov for random walk MH-MCMC step
    cov_rw = cov(theta_tilde);
    
    %duplicate good particles
    r = randsample(1:n_particles,n_particles,'true',particle_w);
    theta_tilde = theta_tilde(r,:);
    log_likelihood = log_likelihood(r);
    particle_w = ones(n_particles,1)/n_particles;
    



    %MUTATION

    %keep track of acceptances
    accept = zeros(n_particles,1);
    accepted_ratio = 0;
    mutation = 0;
    
    % while there are too many duplicates, keep applying MCMC
    while accepted_ratio < unique_particles_required %for each MCMC step
        
        parfor j=1:n_particles %for each particle
            reject = false; %flag for parameters that give errors when modelling
                            
            %propose new values
            theta_tilde_prop = mvnrnd(theta_tilde(j,:),cov_rw*MCMCmult^mutation);
            
            %transform back for likelihood calc (uniform priors only)
            theta_prop = zeros(1,n_params);
            for k=1:n_params
                theta_prop(k) = (prior_lowers(k)+prior_uppers(k)*exp(theta_tilde_prop(k)))./(1+exp(theta_tilde_prop(k)));
            end
            
            %test to see if this parameter combination produces errors in the model
            try
                model_output(theta_prop,dataset);
            catch %if it throws an error - then set reject flag to true
                reject = true;
            end
            
            if reject == false %only try this parameter value if it will not throw an error
                
                %calculate log likelihood for proposed value
                log_like_prop = calc_log_likelihood(theta_prop,dataset,observed_output,model_output);
                
                %calculate priors
                log_prior_curr = sum(theta_tilde(j,:)-2*log(exp(theta_tilde(j,:))+1));
                log_prior_prop = sum(theta_tilde_prop-2*log(exp(theta_tilde_prop)+1));
                
                %MH accept or reject particle
                mh = (exp(newgamma*(log_like_prop - log_likelihood(j)) + log_prior_prop - log_prior_curr));
                if rand < mh
                    %accept particle and update loglikelihood
                    theta_tilde(j,:) = theta_tilde_prop;
                    log_likelihood(j) = log_like_prop;
                    
                    %change acceptance flag
                    accept(j) = 1;
                end
            end %finish error check
            
        end %for each particle
        
        accepted_ratio = sum(accept)/n_particles;
        mutation = mutation+1;
        fprintf('MCMC Step %d: Percentage of accepted particles is %f\n',mutation,accepted_ratio);
        
    end  %for each MCMC step
    
    %convert back to theta for reweighting
    for k=1:n_params
        theta(:,k) = (prior_lowers(k)+prior_uppers(k)*exp(theta_tilde(:,k)))./(1+exp(theta_tilde(:,k)));
    end
    
    %output the results (useful if running this on HPC)
    save('Results.mat')
    
    %update the temperature
    gamma_t = newgamma;
    
end %temperature step end

posterior_sample = theta;
fprintf('SMC algorithm complete \n')

end