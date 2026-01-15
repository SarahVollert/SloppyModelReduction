function S_H = SensitivityMatrix_Hessian(dataset, MLE, delta, observed_output, model_output)
% This function calculates the Hessian matrix evaluated at a single set
% of parameter values (specified as MLE). Here the Hessian is evaluated
% using central finite differencing, with a step size of delta. 

%extract information
[~,n_params] = size(MLE);
sigma = MLE(n_params);

%initialise Hessian matrix
Hessian = zeros(n_params-1);

%for each element of the Hessian
for i=1:n_params-1   %don't calculate for sigma!
    parfor j=1:i     %Hessian is symmetric (only need to calculate half)
        
        % Calculate parameter values for central differencing
        theta_i_up_j_up = MLE;
        theta_i_up_j_up(i) = theta_i_up_j_up(i) + delta*MLE(i)/2;
        theta_i_up_j_up(j) = theta_i_up_j_up(j) + delta*MLE(j)/2;
        
        theta_i_up_j_dw = MLE;
        theta_i_up_j_dw(i) = theta_i_up_j_dw(i) + delta*MLE(i)/2;
        theta_i_up_j_dw(j) = theta_i_up_j_dw(j) - delta*MLE(j)/2;
        
        theta_i_dw_j_up = MLE;
        theta_i_dw_j_up(i) = theta_i_dw_j_up(i) - delta*MLE(i)/2;
        theta_i_dw_j_up(j) = theta_i_dw_j_up(j) + delta*MLE(j)/2;
        
        theta_i_dw_j_dw = MLE;
        theta_i_dw_j_dw(i) = theta_i_dw_j_dw(i) - delta*MLE(i)/2;
        theta_i_dw_j_dw(j) = theta_i_dw_j_dw(j) - delta*MLE(j)/2;
        
        %test for errors (implemented due to pH issue for calcification
        %model)
        try
            model_output(theta_i_up_j_up,dataset);
            model_output(theta_i_dw_j_up,dataset);
            model_output(theta_i_up_j_dw,dataset);
            model_output(theta_i_dw_j_dw,dataset);
        catch
            fprintf('pH issue when trying to calculate the calcification')
        end
        
        %calculate the cost function (likelihood) for parameter changes
        Cost_i_up_j_up = calc_cost_function(theta_i_up_j_up, dataset, observed_output, sigma, model_output)
        Cost_i_up_j_dw = calc_cost_function(theta_i_up_j_dw, dataset, observed_output, sigma, model_output)
        Cost_i_dw_j_up = calc_cost_function(theta_i_dw_j_up, dataset, observed_output, sigma, model_output)
        Cost_i_dw_j_dw = calc_cost_function(theta_i_dw_j_dw, dataset, observed_output, sigma, model_output)

        % Calculate the Hessian entries
        Hessian(i,j) = (Cost_i_up_j_up - Cost_i_up_j_dw - Cost_i_dw_j_up + Cost_i_dw_j_dw) / delta^2;
    end
end

% Copy the symmetric elements to complete sensitivity matrix
S_H = Hessian'+triu(Hessian',1)';

end