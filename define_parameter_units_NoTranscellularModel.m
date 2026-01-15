function units = define_parameter_units_NoTranscellularModel(n_params)
% This function is used to specify the units for each of the model
% parameters, for plotting in figures. This model does not include the
% passive transcullular pathway for CO2. 

% The below units are for the parameters of the coral calcification model
% proposed by Galli & Solidoro (2018) excluding the transcellular pathway
% mechanism. 

units = cell(n_params,1);
units(1) ={'$k_{pp}$ (cm~s$^{-1}$)'}; 
units(2) ={'$s$ (cm~s$^{-1}$)'};
units(3) ={'$\alpha$ (unitless)'}; 
units(4) ={'$\beta$ (unitless)'}; 
units(5) ={'$v_{H_{c}}$ (cm~s$^{-1}$)'}; 
units(6) ={'$E0_{c}$ ($\mu$mol~cm$^{-2}$)'};
units(7) ={'$k1f_c$ (cm$^{4}$~s~$\mu$mol$^{-2}$)'}; 
units(8) ={'$k2f_c$ (cm$^{2}$~$\mu$mol$^{-1}$)'}; 
units(9) ={'$k3f_c$ (s$^{-1}$)'}; 
units(10) ={'$k1b_c$ (s$^{-1}$)'};
units(11) ={'$k2b_c$ (s$^{-1}$)'}; 
units(12) ={'$k3b_c$ (cm$^{4}$~s~$\mu$mol$^{-2}$)'}; 
units(13) ={'$E0_b$ ($\mu$mol~cm$^{-2}$)'}; 
units(14) ={'$k1f_b$ (cm$^{3}$~s$^{-1}~\mu$mol$^{-1}$)'};
units(15) ={'$k2f_b$ (s$^{-1}$)'}; 
units(16) ={'$k3f_b$ (s$^{-1}$)'}; 
units(17) ={'$k1b_b$ (s$^{-1}$)'}; 
units(18) ={'$k2b_b$ (s$^{-1}$)'};
units(19) ={'$k3b_b$ (cm$^{3}$~s$^{-1}~\mu$mol$^{-1}$)'}; 
units(20) ={'$\sigma$ ($\mu$mol~cm$^{-2}$~h$^{-1}$)'};

end