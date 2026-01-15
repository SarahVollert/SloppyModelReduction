function S_L = SensitivityMatrix_likelihoodInformedSubspace(dataset, posterior_sample, prior_sample, delta, observed_output, model_output)
% This function calculates the LIS matrix evaluated for a posterior
% sample. Here the Hessian is evaluated using central finite differencing, 
% with a step size of delta. 

%initialise variables
n_params = size(posterior_sample,2)-1;

% Sloppy analysis samples
posterior_sample_sloppy = posterior_sample;
prior_sample_sloppy = log(prior_sample(:,1:n_params));

%find Cholesky factor L
Omega = cov(prior_sample_sloppy);
L = chol(Omega);

%initialise matrices
S = zeros(n_params,n_params);
err_count = 0;

%for each particle in the sample - approximate the Hessian
parfor m=1:size(posterior_sample_sloppy,1)
    particle_error=false;
    sigma = posterior_sample(m,n_params+1);
    
    %calculate hessian matrix
    H = zeros(n_params,n_params);
    for i=1:n_params
        for j=1:i
            % Calculate parameter values
            theta_i_up_j_up = posterior_sample_sloppy(m,:);
            theta_i_up_j_up(i) = theta_i_up_j_up(i) + delta*posterior_sample_sloppy(m,i)/2;
            theta_i_up_j_up(j) = theta_i_up_j_up(j) + delta*posterior_sample_sloppy(m,j)/2;
            
            theta_i_up_j_dw = posterior_sample_sloppy(m,:);
            theta_i_up_j_dw(i) = theta_i_up_j_dw(i) + delta*posterior_sample_sloppy(m,i)/2;
            theta_i_up_j_dw(j) = theta_i_up_j_dw(j) - delta*posterior_sample_sloppy(m,j)/2;
            
            theta_i_dw_j_up = posterior_sample_sloppy(m,:);
            theta_i_dw_j_up(i) = theta_i_dw_j_up(i) - delta*posterior_sample_sloppy(m,i)/2;
            theta_i_dw_j_up(j) = theta_i_dw_j_up(j) + delta*posterior_sample_sloppy(m,j)/2;
            
            theta_i_dw_j_dw = posterior_sample_sloppy(m,:);
            theta_i_dw_j_dw(i) = theta_i_dw_j_dw(i) - delta*posterior_sample_sloppy(m,i)/2;
            theta_i_dw_j_dw(j) = theta_i_dw_j_dw(j) - delta*posterior_sample_sloppy(m,j)/2;
            
            % Evaluate the cost function
            try
                model_output(theta_i_up_j_up,dataset);
                model_output(theta_i_dw_j_up,dataset);
                model_output(theta_i_up_j_dw,dataset);
                model_output(theta_i_dw_j_dw,dataset);
            catch
                particle_error=true
                err_count = err_count +1;
            end
            
            if particle_error == true
                break
            end

            %calculate the cost function (likelihood) for parameter changes
            Cost_i_up_j_up = calc_cost_function(theta_i_up_j_up, dataset, observed_output, sigma, model_output)
            Cost_i_up_j_dw = calc_cost_function(theta_i_up_j_dw, dataset, observed_output, sigma, model_output)
            Cost_i_dw_j_up = calc_cost_function(theta_i_dw_j_up, dataset, observed_output, sigma, model_output)
            Cost_i_dw_j_dw = calc_cost_function(theta_i_dw_j_dw, dataset, observed_output, sigma, model_output)

            % Calculate the Hessian entries
            H(i,j) = (Cost_i_up_j_up - Cost_i_up_j_dw - Cost_i_dw_j_up + Cost_i_dw_j_dw) / delta^2;
        end
        
        if particle_error == true %included because of coral calcification pH error
            break
        end
        
    end
    
    %if there is no error caused by particle
    if particle_error == false
        %complete the hessian
        H_matrix = H'+triu(H',1)';
        
        %calculate sum-term
        S = S + transpose(L)*H_matrix*L;
    end
    
    fprintf('completed particle no. %d \n',m)
end

%calculate the LIS matrix by dividing by the number of particles used
S_L = S/(length(posterior_sample_sloppy)-err_count);

end