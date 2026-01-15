function units = define_parameter_units_NoBATModel(n_params)
% This function is used to specify the units for each of the model
% parameters, for plotting in figures. This model does not include the BAT
% pump mechanism. 

% The below units are for the parameters of the coral calcification model
% proposed by Galli & Solidoro (2018) excluding the BAT mechanism.

units = cell(n_params,1);
units(1) ={'$k_{CO_2}$ (cm~s$^{-1}$)'}; 
units(2) ={'$k_{pp}$ (cm~s$^{-1}$)'}; 
units(3) ={'$s$ (cm~s$^{-1}$)'};
units(4) ={'$\alpha$ (unitless)'}; 
units(5) ={'$\beta$ (unitless)'}; 
units(6) ={'$v_{H_{c}}$ (cm~s$^{-1}$)'}; 
units(7) ={'$E0_{c}$ ($\mu$mol~cm$^{-2}$)'};
units(8) ={'$k1f_c$ (cm$^{4}$~s~$\mu$mol$^{-2}$)'}; 
units(9) ={'$k2f_c$ (cm$^{2}$~$\mu$mol$^{-1}$)'}; 
units(10) ={'$k3f_c$ (s$^{-1}$)'}; 
units(11) ={'$k1b_c$ (s$^{-1}$)'};
units(12) ={'$k2b_c$ (s$^{-1}$)'}; 
units(13) ={'$k3b_c$ (cm$^{4}$~s~$\mu$mol$^{-2}$)'}; 
units(14) ={'$\sigma$ ($\mu$mol~cm$^{-2}$~h$^{-1}$)'};

end