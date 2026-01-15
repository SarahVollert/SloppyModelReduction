function neg_loglike = calc_nll_with_error_check(parameter_vals, dataset, observed_output,model_output)
%This function calculates the negative loglikelihood given the parameter
%values and available data. 
%Note: due to errors within the pH algorithm for the calcification model, a
%try-catch error check is built into this function. 


%set an error flag to test paramter estimates
error = false;

%try catch in case tested theta causes an error
try
    model_output(parameter_vals,dataset);
catch
    error = true;
    neg_loglike = NaN;
end

%if no error, calculate the nll
if error==false
    neg_loglike = -calc_log_likelihood(parameter_vals, dataset, observed_output,model_output);
end

end