%% model_output_NoTranscellularModel.m
% This function defines the model. Here, we input the value of the
% parameters, and the input data, and the function returns the model's
% output. 

% The below code defines the coral calcification model proposed by Galli &
% Solidoro (2018) with the transcellular pathway mechanism removed. 
% This function takes in parameter values for the
% unknown model parameters and sigma (noise parameter), and the
% environmental conditions for the coral (dataset). This function then
% outputs the calcification predictions from the model, given these
% parameter values and environmental data. 

function model_predictions = model_output_NoTranscellularModel(parameter_vals,dataset)


%% Calcification model

%getting parameter values from the inputs:
p.k_pp = parameter_vals(1); p.s = parameter_vals(2); 
p.alpha=parameter_vals(3); p.beta=parameter_vals(4); 
p.v_H_c=parameter_vals(5); p.E0_c=parameter_vals(6); 
p.k1f_c=parameter_vals(7); p.k2f_c=parameter_vals(8); 
p.k3f_c=parameter_vals(9); p.k1b_c=parameter_vals(10); 
p.k2b_c=parameter_vals(11); p.k3b_c=parameter_vals(12); 
p.E0_b=parameter_vals(13); p.k1f_b=parameter_vals(14); 
p.k2f_b=parameter_vals(15); p.k3f_b=parameter_vals(16); 
p.k1b_b=parameter_vals(17); p.k2b_b=parameter_vals(18); 
p.k3b_b=parameter_vals(19); p.sigma=parameter_vals(20);


%known parameters (as given by Galli & Solidoro, 2018)
p.Ca_sw = 10.4e3;             % taken from https://www.iaea.org/ocean-acidification/download/OA%20guide/OA-Guide_Ch1_20100518.pdf
p.h_coel = 3000/ (10^4);      %[cm] coelenteron thickness
p.h_ecm = 5/ (10^4);          %[cm] calcyfying medium thickness
p.k_p = (1.1)/ (10^3);        %[umol cm kg^-1 s^-1]!!!!!  aragonite precipitation rate constant
p.n_p = 1.63;                 %[-] reaction order for precipitation
p.k_d = (27)/ (10^3);         %[umol cm kg^-1 s^-1]  aragonite dissolution rate constant
p.n_d = 2.5;                  %[-] reaction order for dissolution
p.H_O2 = 473;                 %[kJ mol^-1] Energy equivalent of Carbs
p.delta_G_ATP = 30.5;         %[kJ/mol]
p.Sal = 38.1;                 %[psu] Salinity

%extracting environmental conditions from the input dataset
Temp_data = dataset(:,1); % deg C
TA_sw_data = dataset(:,2); % umol kg^-1
DIC_sw_data = dataset(:,3); % umol kg^-1
Pg_data = dataset(:,4); % nmol cm^-2 h^-1
R_data = -dataset(:,5); % nmol cm^-2 h^-1

%initialising 
[n_data,~]=size(dataset);               %number of data points
model_predictions = zeros(1,n_data);    %predicted calcification rates
tspan = [0 1500];             %run up to 2000 seconds - assume equilibrium
y0 = 1500*ones(6,1);                    % An arbitrary initial guess

%predict the calcification rate for each datapoint
for n = 1:n_data

    %extract environmental conditions for this data point
    p.DIC_sw = DIC_sw_data(n);                      
    p.TA_sw = TA_sw_data(n);                       
    p.Pg = Pg_data(n) / 3600; % umol cm kg^-1 s^-1  converted from nmol cm^-2 h^-1
    p.R  = R_data(n) / 3600;  % umol cm kg^-1 s^-1  converted from nmol cm^-2 h^-1
    p.Temp = Temp_data(n);
    
    %solve system of ODEs
    [~,yF] = ode15s(@(t,y) Galli_ODEs(t,y,p),tspan,y0);
    
    %take calcification rates as steady state solution
    model_predictions(n) = calculate_calcification(yF(end,:),p)*3600;
    % convert units from umol cm kg^-1 s^-1 to nmol cm^-2 h^-1.
end

end


%% Galli_ODEs function
function [dydt] = Galli_ODEs(~,y,p)
% This function defines the system of ODEs as proposed by Galli & Solidoro
% (2018). This function takes in parameter values and environmental
% conditions using the structure p. These ODEs model the DIC, TA and Ca
% concentration within the coelenteron and ECM, using the variable y. 

% 1. Organise inputs
DIC_coel = y(1);
TA_coel = y(2);
Ca_coel = y(3);
DIC_ecm = y(4);
TA_ecm = y(5);
Ca_ecm = y(6);

% 2. Model
% 2.1 Calculating carbon equilibrium

K_ar = Aragonite_constant(p.Temp,p.Sal);
% K_ar is the solubility product of aragonite calculated from temperature
% and salinity according to Zeebe and Wolf-Gladrow (2001).

[pH_coel,CO2_coel,HCO3_coel,~] = pH_algorithm(DIC_coel,TA_coel,p.Temp,p.Sal);
Hplus_coel = 10^(-pH_coel + 6); % +6 required to convert mol/kg to umol/kg.
[pH_ecm,CO2_ecm,HCO3_ecm,CO3_ecm] = pH_algorithm(DIC_ecm,TA_ecm,p.Temp,p.Sal);
Hplus_ecm = 10^(-pH_ecm + 6); % +6 required to convert mol/kg to umol/kg

% H+, OH-, CO2, HCO3 and CO3 are calculated in each compartment from DIC,
% TA, temperature and salinity, assuming that chemical equilibrium is
% reached at each time step, according to the equilibrium relations in
% Zeebe and Wolf-Gladrow (2001). For the sake of simplicity, we consider
% just the carbonates' contribution to alkalinity and neglect, e.g.
% borates.

% 2.2 Calculating all fluxes

% Bicarbonate anion transport (BAT), directed from coelenteron to ECM.
den_b = p.k2f_b*p.k3f_b + p.k2f_b*p.k3b_b*HCO3_ecm + p.k1f_b*p.k2f_b*HCO3_coel ...
        + p.k1f_b*p.k3f_b*HCO3_coel + p.k2b_b*p.k3b_b*HCO3_ecm + p.k1b_b*p.k3b_b*HCO3_ecm ...
        + p.k1f_b*p.k2b_b*HCO3_coel + p.k1b_b*p.k2b_b + p.k1b_b*p.k3f_b;                      
    
J_BAT = p.E0_b * (p.k1f_b*p.k2f_b*p.k3f_b*HCO3_coel - p.k1b_b*p.k2b_b*p.k3b_b*HCO3_ecm)/den_b;

% Passive transport, directed from seawater to coelenteron.
J_sw_coel_DIC = p.s * (p.DIC_sw - DIC_coel);                                                  
J_sw_coel_TA = p.s * (p.TA_sw - TA_coel);                                                
J_sw_coel_Ca = p.s * (p.Ca_sw - Ca_coel);                                                   

% Passive transport, directed from coelenteron to ECM. 
J_coel_ecm_DIC = p.k_pp * (DIC_coel - DIC_ecm);                                               
J_coel_ecm_TA = p.k_pp * (TA_coel - TA_ecm);                                                
J_coel_ecm_Ca = p.k_pp * (Ca_coel - Ca_ecm);                                                   

% ATP production fluxes associated with photosynthesis and respiration.
J_Pg_ATP = p.H_O2 * p.Pg / p.delta_G_ATP;                                                    
J_R_ATP = p.H_O2 * p.R / p.delta_G_ATP;                                                      

% ATP flux used to run active transports, due to photosynthesis and
% respiration.
J_ATPCa_ATPase = p.alpha*J_Pg_ATP + p.beta*J_R_ATP;                                          

% Active transport due to Ca-ATPase pump, from coelenteron to ECM.
den_c = p.k2f_c*p.k3f_c*J_ATPCa_ATPase + p.k2f_c*p.k3b_c*(p.v_H_c*Hplus_coel)^2*J_ATPCa_ATPase ...
       + p.k1f_c*p.k2f_c*(p.v_H_c*Hplus_ecm)^2*J_ATPCa_ATPase + p.k1f_c*p.k3f_c*(p.v_H_c*Hplus_ecm)^2 ...
       + p.k2b_c*p.k3b_c*(p.v_H_c*Hplus_coel)^2 + p.k1b_c*p.k3b_c*(p.v_H_c*Hplus_coel)^2 ...
       + p.k1f_c*p.k2b_c*(p.v_H_c*Hplus_ecm)^2 + p.k1b_c*p.k2b_c + p.k1b_c*p.k3f_c;           
   
J_CaATPase = p.E0_c * (p.k1f_c*p.k2f_c*p.k3f_c*J_ATPCa_ATPase * (p.v_H_c*Hplus_ecm)^2 ...
                    - p.k1b_c*p.k2b_c*p.k3b_c*(p.v_H_c*Hplus_coel)^2 ) / den_c;               

% Aragonite precipitation and dissolution kinetics
Omega = Ca_ecm * CO3_ecm / K_ar;

if Omega >= 1
    J_CaCO3 = p.k_p*(Omega - 1)^p.n_p;                                                       
else
    J_CaCO3 = p.k_d*(1 - Omega)^p.n_d;                                                      
end

% 3. Final ODEs
dDIC_coel_dt = (-J_BAT + J_sw_coel_DIC - J_coel_ecm_DIC - p.Pg + p.R)/p.h_coel; 
dTA_coel_dt = (-J_BAT - 2*J_CaATPase + J_sw_coel_TA - J_coel_ecm_TA)/p.h_coel;          
dCa_coel_dt = (-J_CaATPase + J_sw_coel_Ca - J_coel_ecm_Ca)/p.h_coel;                   
dDIC_ecm_dt = (J_BAT + J_coel_ecm_DIC  - J_CaCO3)/p.h_ecm;                  
dTA_ecm_dt = (J_BAT + 2*J_CaATPase + J_coel_ecm_TA - 2*J_CaCO3)/p.h_ecm;           
dCa_ecm_dt = (J_CaATPase + J_coel_ecm_Ca - J_CaCO3)/p.h_ecm;                            

% 4. Organise outputs
dydt = [dDIC_coel_dt;
        dTA_coel_dt;
        dCa_coel_dt;
        dDIC_ecm_dt;
        dTA_ecm_dt;
        dCa_ecm_dt];

end

%% Aragonite_constant function
function K_ar = Aragonite_constant(Temp,S)

% This function takes in the temperature and the salinity and uses it to
% caclulate the aragonite constant for these environmental conditions. 

% This is from Mucci (1983), ref'd in Appendix A of Zeebe & Wolf-Gladrow
% (2001) Appendix A10. 


T = Temp + 273.15; %convert to kelvin
log_K_arg = -171.945 - 0.077993*T + 2903.293/T ...
            + 71.595*log10(T) ...
            + (-0.068393 + 0.0017276*T + 88.135/T) * sqrt(S) ...
            - 0.10018*S + 0.0059415* S^1.5;

K_ar = 10^(log_K_arg+12);
% The additional exponent +12 makes K_arg in units of umol^2 kg^-2.
end

%% pH_algorithm function
function [pH_final,CO2_final,HCO3_minus_final,CO3_2minus_final] = pH_algorithm(DIC,TA,Temp,Sal)
% This algorithm takes in the chemical concentrations (DIC, TA) and
% environmental conditions (temperature and salinity) and uses it to
% calculate dependent species concentrations and pH level. 

% Algorithm based on Verspagen et al. 2014 Ecology Letters 17:951-960, Appendix S1


% set up inputs
DIC = DIC/1e6;          % Convert to mol/L.
TA = TA/1e6;            % Convert to mol/L.
pH_tolerance = 1e-6;
T = Temp + 273.15;      % Convert to Kelvins
S = Sal;                % g/L or psu

% pKw, from Millero 1995, Geochimica et Cosmochimica Acta 59: 661-677, 
% equation (63).
pKw = (148.9802 - 13847.26/T - 23.6521*log(T) ...
    + (-5.977 + 118.67/T + 1.0495*log(T))*sqrt(S) ...
    - 0.01615*S) / -log(10);

% pK1 and pK2, from Millero et al. 2006, Marine Chemistry 100: 80-94,
% equations in the abstract.
pK10 = -126.34048 + 6320.813/T + 19.568224*log(T);
pK20 = -90.18333  + 5143.692/T + 14.613358*log(T);
A1 = 13.4191*sqrt(S) + 0.0331*S - 5.33e-5 * S^2; 
B1 = -530.123*sqrt(S) - 6.103*S;
C1 = -2.06950*sqrt(S);
A2 = 21.0894*sqrt(S) + 0.1248*S - 3.687e-4 * S^2;
B2 = -772.483*sqrt(S) - 20.051*S;
C2 = -3.3336*sqrt(S);
pK1star = A1 + B1/T + C1*log(T);
pK2star = A2 + B2/T + C2*log(T);
pK1 = pK10 + pK1star;
pK2 = pK20 + pK2star;

% Model parameters
K_w = 10^(-pKw);
K_1 = 10^(-pK1);
K_2 = 10^(-pK2);

pH = 7; % First guess

%initialise algorithm variables
pH_value_found = false; %flag for finding pH value
m=0; %number of iterations
Iter_Tol = 200; %how many search iterations are allowed before non-convergence declared

%Search for pH value
while pH_value_found == false
    
    H_plus = 10^(-pH); % mol/L
    OH_minus = K_w/H_plus; % mol/L
    
    alpha_c = H_plus^2 + K_1*H_plus + K_1*K_2;
    
    alpha_0 = H_plus^2/alpha_c;
    alpha_1 = K_1*H_plus/alpha_c;
    alpha_2 = K_1*K_2/alpha_c;
    
    CO2 = alpha_0 * DIC;
    HCO3_minus = alpha_1 * DIC;
    CO3_2minus = alpha_2 * DIC;
    
    Buffer_Capacity = (    H_plus + OH_minus + ...
        (alpha_1 * (alpha_0+alpha_2) + 4*alpha_0*alpha_2) * DIC) * log(10);
    
    ALK_test = HCO3_minus + 2 * CO3_2minus + OH_minus - H_plus;
    
    Delta_pH = (TA - ALK_test)/Buffer_Capacity;
    
    pH = pH + Delta_pH; % New estimate of pH
    
    if abs(Delta_pH) < pH_tolerance %If the change in pH is minimal (below our tolerance)
        break
    end
    
    m=m+1;
    if m>Iter_Tol %If we've searched too long
        error('The pH was set to -5 because the algorithm didnt converge')

    end
    if pH<=0 || pH>=14  %If pH values are invalid
        error('the pH algorithm is not converging properly')
    end
end

%return the outputs once pH value found
pH_final = pH;
CO2_final = CO2*1e6;
HCO3_minus_final = HCO3_minus*1e6;
CO3_2minus_final = CO3_2minus*1e6;
end

%% calculate_calcification function
function J_CaCO3 = calculate_calcification(y,p)
% This function calculates the calcification rate, given the environmental
% conditions and the solved ODE system. 

DIC_ecm = y(4);
TA_ecm = y(5);
Ca_ecm = y(6);

K_ar = Aragonite_constant(p.Temp,p.Sal);
% K_ar is the solubility product of aragonite calculated from temperature
% and salinity according to Zeebe and Wolf-Gladrow (2001).

[~,~,~,CO3_ecm] = pH_algorithm(DIC_ecm,TA_ecm,p.Temp,p.Sal);
   
% Aragonite precipitation and dissolution kinetics
Omega = Ca_ecm * CO3_ecm / K_ar;

if Omega >= 1
    J_CaCO3 = p.k_p*(Omega - 1)^p.n_p;       %Calcification precipitation                                           
else
    J_CaCO3 = p.k_d*(1 - Omega)^p.n_d;       %Calcification dissolution                      
end

end

%%

