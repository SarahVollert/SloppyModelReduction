%% define_uniform_priors_NoBATNoTranscellularModel.m
function [Prior_uppers, Prior_lowers] = define_uniform_priors_NoBATNoTranscellularModel

% This function is used to define all of the uniform prior bounds for a
% SMC sampling algorithm. This function returns a structure for the upper
% bounds and lower bounds separately. 
% 
% The function is currently written to define the
% prior bounds for the coral calcification model proposed by Galli &
% Solidoro (2018) where the BAT and transcellular pathway mechanisms are excluded.  

%% Uniform prior bounds for coral calcification model:

%Define the prior upper bounds
k_pp_upper = 0.1;             %physiological
s_upper = 0.1;                %physiological
alpha_upper = 1;              %fixed
beta_upper = 1;               %fixed
v_H_c_upper = 250;            
E0_c_upper = 12e6;      
k1f_c_upper =140e-6;          
k2f_c_upper =500e-3;          
k3f_c_upper =800;             
k1b_c_upper =8;               
k2b_c_upper =500;             
k3b_c_upper =10e-8;                  
sigma_upper = 50;               

Prior_uppers = [k_pp_upper; s_upper;
    alpha_upper; beta_upper; v_H_c_upper; E0_c_upper;
    k1f_c_upper; k2f_c_upper; k3f_c_upper; k1b_c_upper;
    k2b_c_upper; k3b_c_upper; sigma_upper];

%Define the prior lower bounds
k_pp_lower = 0;               %fixed
s_lower =0;                   %fixed
alpha_lower = 0;              %fixed
beta_lower = 0;               %fixed
v_H_c_lower =0;               %fixed
E0_c_lower = 0;               %fixed
k1f_c_lower =0;               %fixed
k2f_c_lower =0;               %fixed
k3f_c_lower =0;               %fixed
k1b_c_lower =0;               %fixed
k2b_c_lower =0;               %fixed
k3b_c_lower =0;               %fixed
sigma_lower = 0;              %fixed

Prior_lowers = [k_pp_lower; s_lower;
    alpha_lower; beta_lower; v_H_c_lower; E0_c_lower;
    k1f_c_lower; k2f_c_lower; k3f_c_lower; k1b_c_lower;
    k2b_c_lower; k3b_c_lower; sigma_lower];

end