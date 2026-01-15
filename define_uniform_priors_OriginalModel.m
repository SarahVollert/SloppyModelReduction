%% define_uniform_priors_OriginalModel.m
function [Prior_uppers, Prior_lowers] = define_uniform_priors_OriginalModel

% This function is used to define all of the uniform prior bounds for a
% SMC sampling algorithm. This function returns a structure for the upper
% bounds and lower bounds separately. 
% 
% The function is currently written to define the
% prior bounds for the coral calcification model proposed by Galli &
% Solidoro (2018). 

%% Uniform prior bounds for coral calcification model:

%Define the prior upper bounds
k_CO2_upper = 0.1;            %physiological
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
E0_b_upper =1500;             
k1f_b_upper =50e-6;           
k2f_b_upper =0.01;            
k3f_b_upper =0.01;           
k1b_b_upper =20e-5;           
k2b_b_upper =10e-4;           
k3b_b_upper =35e-9;           
sigma_upper = 50;               

Prior_uppers = [k_CO2_upper; k_pp_upper; s_upper;
    alpha_upper; beta_upper; v_H_c_upper; E0_c_upper;
    k1f_c_upper; k2f_c_upper; k3f_c_upper; k1b_c_upper;
    k2b_c_upper; k3b_c_upper; E0_b_upper; k1f_b_upper;
    k2f_b_upper; k3f_b_upper; k1b_b_upper; k2b_b_upper;
    k3b_b_upper; sigma_upper];

%Define the prior lower bounds
k_CO2_lower = 0;              %fixed
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
E0_b_lower =0;                %fixed
k1f_b_lower =0;               %fixed
k2f_b_lower =0;               %fixed
k3f_b_lower =0;               %fixed
k1b_b_lower =0;               %fixed
k2b_b_lower =0;               %fixed
k3b_b_lower =0;               %fixed
sigma_lower = 0;              %fixed

Prior_lowers = [k_CO2_lower; k_pp_lower; s_lower;
    alpha_lower; beta_lower; v_H_c_lower; E0_c_lower;
    k1f_c_lower; k2f_c_lower; k3f_c_lower; k1b_c_lower;
    k2b_c_lower; k3b_c_lower; E0_b_lower; k1f_b_lower;
    k2f_b_lower; k3f_b_lower; k1b_b_lower; k2b_b_lower;
    k3b_b_lower; sigma_lower];

end