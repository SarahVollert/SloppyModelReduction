function Plots_GoodnessOfFit(posterior_prediction, observed_output)
%This function plots the observed output against the modelled output as a
%measure of goodness of fit. 
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
plot(observed_output,median_output,'.k','MarkerSize',15);

%error bar plots
x = observed_output;
y = median_output;
yneg = median_output-CI95_lower_bound;
ypos = CI95_upper_bound-median_output;
errorbar(x,y,yneg,ypos,'.k')

%plot x=y line
plot([0 300],[0 300],'Color',[0.5 0.5 0.5])

% add legend
legend('Median Calcification Prediction','95% Central Credible Interval','Location','northwest')

end