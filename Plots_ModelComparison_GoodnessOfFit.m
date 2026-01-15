function Plots_ModelComparison_GoodnessOfFit(Original,NoBat, NoTrans, NoBatNoTrans, observed_output)
CL = get(gca,'colororder');

%% Plot No Trans No BAT model against original model
%set up plot
figure
pbaspect([1 1 1])
xlabel('Calcification Observations (nmol cm^{-2} h^{-1})')
ylabel('Model Calcification (nmol cm^{-2} h^{-1})')
hold on

% Plot reduced model (NoBatNoTrans)
%sort posterior calcficiations 
sorted_calcification = sort(NoBatNoTrans.posterior_prediction);
median_calcification = sorted_calcification(round(length(NoBatNoTrans.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(NoBatNoTrans.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(NoBatNoTrans.posterior_prediction)*0.975),:)';
%plot
median_plot1 = plot(observed_output,median_calcification,'.','Color','#ba5eba','MarkerSize',15);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.','Color','#ba5eba','LineWidth',1.5)

% Plot original model 
%sort posterior calcficiations
sorted_calcification = sort(Original.posterior_prediction);
median_calcification = sorted_calcification(round(length(Original.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(Original.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(Original.posterior_prediction)*0.975),:)';
%original
median_plot2 = plot(observed_output,median_calcification,'.k','MarkerSize',10);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.k','LineWidth',0.75)

%plot x=y line
plot([0 300],[0 300],'Color',[0.5 0.5 0.5])

% add legend
median_plot1.Annotation.LegendInformation.IconDisplayStyle = 'off';
median_plot2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reduced model','Original model','Location','northwest')


%% Plot original and No BAT
%set up plot
figure
pbaspect([1 1 1])
xlabel('Calcification Observations (nmol cm^{-2} h^{-1})')
ylabel('Model Calcification (nmol cm^{-2} h^{-1})')
hold on
 
% Plot reduced model (NoBat)
%sort posterior calcficiations
sorted_calcification = sort(NoBat.posterior_prediction);
median_calcification = sorted_calcification(round(length(NoBat.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(NoBat.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(NoBat.posterior_prediction)*0.975),:)';
%plot 
median_plot1 = plot(observed_output,median_calcification,'.','Color',CL(2,:),'MarkerSize',15);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.','Color',CL(2,:),'LineWidth',1.5)

% Plot original model
%sort posterior calcficiations
sorted_calcification = sort(Original.posterior_prediction);
median_calcification = sorted_calcification(round(length(Original.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(Original.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(Original.posterior_prediction)*0.975),:)';
%plot
median_plot2 = plot(observed_output,median_calcification,'.k','MarkerSize',10);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.k','LineWidth',0.75)

%plot x=y line
plot([0 300],[0 300],'Color',[0.5 0.5 0.5])

% add legend
median_plot1.Annotation.LegendInformation.IconDisplayStyle = 'off';
median_plot2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('No BAT mechanism model','Original model','Location','northwest')


%% Plot No trans model to original
%set up plot
figure
pbaspect([1 1 1])
xlabel('Calcification Observations (nmol cm^{-2} h^{-1})')
ylabel('Model Calcification (nmol cm^{-2} h^{-1})')
hold on
 
% Plot reduced model (NoTrans)
%sort posterior predictions
sorted_calcification = sort(NoTrans.posterior_prediction);
median_calcification = sorted_calcification(round(length(NoTrans.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(NoTrans.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(NoTrans.posterior_prediction)*0.975),:)';
%plot
median_plot1 = plot(observed_output,median_calcification,'.','Color',CL(1,:),'MarkerSize',15);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.','Color',CL(1,:),'LineWidth',1.5)

% Plot original model
%sort posterior calcficiations
sorted_calcification = sort(Original.posterior_prediction);
median_calcification = sorted_calcification(round(length(Original.posterior_prediction)*0.5),:)';
CI95_lower = sorted_calcification(ceil(length(Original.posterior_prediction)*0.025),:)';
CI95_upper = sorted_calcification(floor(length(Original.posterior_prediction)*0.975),:)';
%plot
median_plot2 = plot(observed_output,median_calcification,'.k','MarkerSize',10);
errorbar(observed_output,median_calcification,median_calcification-CI95_lower,CI95_upper-median_calcification,'.k','LineWidth',1)

%plot x=y line
plot([0 300],[0 300],'Color',[0.5 0.5 0.5])

% add legend
median_plot1.Annotation.LegendInformation.IconDisplayStyle = 'off';
median_plot2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('No transcellular pathway model','Original model','Location','northwest')

end