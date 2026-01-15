function Plots_eigenvalues_methodComparison(eigenvalues_H, eigenvalues_P, eigenvalues_L)
% This function plots the eigenvalues from an analysis of model sloppiness
% using each of the approaches mentioned in the manuscript (a Hessian
% evaluated at a maximum likelihood estimate, the posterior covariance
% method and the likelihood informed subspace method). 

% Note: this multiple analyses of model sloppiness can be added here where
% each column of the matrix is a different analysis. 

%note the number of each type of analyses
[n_params,n_analyses_H] = size(eigenvalues_H);
[~,n_analyses_P] = size(eigenvalues_P);
[~,n_analyses_L] = size(eigenvalues_L);

%prepare figure
figure
hold on
x_ticks = linspace(1,n_params,n_params);
xlim([x_ticks(1)-0.5 x_ticks(end)+0.5])

%plot each MLE 
for i=1:n_analyses_H
    p1 = plot(x_ticks,eigenvalues_H(:,i),'-o','Color','#FDD0A3','MarkerSize',5);
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    p1.MarkerEdgeColor = [0.8500 0.3250 0.0980];
end

%plot each PCA 
for i=1:n_analyses_P
    p2 = plot(x_ticks,eigenvalues_P(:,i),'-s','Color','#D8F3BC','MarkerSize',7);
    p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    p2.MarkerEdgeColor = [0.4660 0.6740 0.1880];
end

%plot each LIS 
for i=1:n_analyses_L
    p3 = plot(x_ticks,eigenvalues_L(:,i),'-v','Color','#B0D8FF','MarkerSize',5);
    p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    p3.MarkerEdgeColor = [0 0.4470 0.7410];
end

%formatting
set(gca,'YScale','log')
xlabel('Eigenvalue index ($j$)','Interpreter','Latex')
ylabel('Relative Eigenvalue Size ($\lambda_j/\lambda_1$)','Interpreter','Latex')
xticks(x_ticks)

%add legend
legend([p1 p2 p3],...
    {'Hessian approach (evaluated at likelihood maxima)','Posterior covariance approach','LIS matrix approach'},...
    'Location','SouthWest')


end