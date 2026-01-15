function Plots_eigenvalues_modelComparison(eigenvalues_original, eigenvalues_reduced)
% This function plots the eigenvalues from an analysis of model sloppiness
% for both the original model and reduced model (without BAT or 
% transcellular mechanisms). 

%extract info
[n_params_original,~] = size(eigenvalues_original);
[n_params_reduced,~] = size(eigenvalues_reduced);
n_params = max([n_params_original,n_params_reduced]);

%prepare figure
figure
hold on
x_ticks = linspace(1,n_params,n_params);
xlim([x_ticks(1)-0.5 x_ticks(end)+0.5]) %align central on tick

%add legend
plot(x_ticks,eigenvalues_original,'.','Color','#a6a6a6','MarkerSize',15)
plot(x_ticks(1:n_params_reduced),eigenvalues_reduced,'.','Color','#ba5eba','MarkerSize',15)
legend('Original model','Reduced model','Location','NorthEast')

%formatting
set(gca,'YScale','log')
xlabel('Eigenvalue index ($j$)','Interpreter','Latex')
ylabel('Relative Eigenvalue Size ($\lambda_j/\lambda_1$)','Interpreter','Latex')
xticks(x_ticks)

end