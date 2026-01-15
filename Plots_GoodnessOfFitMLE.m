function Plots_GoodnessOfFitMLE(posterior_prediction, observed_output, ...
                    MLE_prediction,n_highlight)
%This function plots the observed output against the modelled output as a
%measure of goodness of fit. It also includes the ability to plot MLEs. 

%This plot is currently set up for the calcification model. 

%sort posterior predictions
sorted_output = sort(posterior_prediction);

%find medians and credible intervals for all point predictions
median_output = sorted_output(size(posterior_prediction,1)*0.5,:)';
CI95_lower_bound = sorted_output(ceil(size(posterior_prediction,1)*0.025),:)'; %95% credible interval
CI95_upper_bound = sorted_output(floor(size(posterior_prediction,1)*0.975),:)'; %95% credible interval

%set up plot
figure
xlabel('Calcification Observations (nmol cm^{-2} h^{-1})')
ylabel('Model Calcification (nmol cm^{-2} h^{-1})')
hold on

%plots based on calcification output
median_plot = plot(observed_output,median_output,'.k','MarkerSize',15);

%error bar plots
x = observed_output;
y = median_output;
yneg = median_output-CI95_lower_bound;
ypos = CI95_upper_bound-median_output;
CI_plot = errorbar(x,y,yneg,ypos,'.k');

%plot MLE
all_MLE_plot = plot(observed_output,MLE_prediction,'*k');
highlight_MLE_plot = plot(observed_output,MLE_prediction(1:n_highlight,:),'*','Color',[0.8500 0.3250 0.0980]);

%plot x=y line
plot([0 300],[0 300],'Color',[0.5 0.5 0.5])

% add legend
legend([median_plot CI_plot all_MLE_plot(1) highlight_MLE_plot(1)],...
    {'Median prediction','95% central credible interval',...
    'Local MLE predictions','Selected local MLE predictions'},...
        'Location','northwest')


end