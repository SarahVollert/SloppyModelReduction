function Plots_PosteriorDensityMLE(prior_sample, posterior_sample, Priors, ...
            units, MLE_sample,n_highlight)
% This script plots the estimated posterior and prior densities of a
% sample, similarly to the function Plots_PosteriorDensity but including
% the MLEs. 

%locally extract prior bounds
prior_lowers = Priors.prior_lowers;
prior_uppers = Priors.prior_uppers;
[n_params,~] = size(prior_uppers);

%prepare to rescale 
posteriorSample_scale = posterior_sample;
priorSample_scale = prior_sample;
Prior_lowers_scale = prior_lowers;
Prior_uppers_scale = prior_uppers;
MLE_scale = MLE_sample;

%the rescale order of magnitude for each parameter
scalers = floor(log10(1./prior_uppers));

%rescaling
for i=1:n_params
    posteriorSample_scale(:,i) = posterior_sample(:,i).*10.^scalers(i);
    priorSample_scale(:,i) = prior_sample(:,i).*10.^scalers(i);
    Prior_lowers_scale(i) = prior_lowers(i).*10.^scalers(i);
    Prior_uppers_scale(i) = prior_uppers(i).*10.^scalers(i);
    MLE_scale(:,i) = MLE_sample(:,i).*10.^scalers(i);

    if not(scalers(i) ==0)
        units(i) = strcat(units(i),' $\times10^{',int2str(scalers(i)),'}$');
    end
end

%set up figure
figure;
figure_rows = ceil(n_params/3);
tiles = tiledlayout(figure_rows,3);
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

%for each parameter
for i=1:n_params
    
    %prep figure
    nexttile
    hold on
    
    %plot the prior distribution in grey
    priorPlot = ksdensity(priorSample_scale(:,i),'BoundaryCorrection','reflection','Support',[Prior_lowers_scale(i) Prior_uppers_scale(i)]);
    priorPlot(1) = priorPlot(2);priorPlot(end)=priorPlot(end-1);
    stepSize = (Prior_uppers_scale(i) - Prior_lowers_scale(i))/(length(priorPlot)-1);
    prior_plot = area(Prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),priorPlot,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    
    %plot the MLEs
    for j=1:size(MLE_scale,1)
        all_MLE_plot = xline(MLE_scale(j,i),'LineWidth',0.25);
        if j<=n_highlight
            highlight_MLE_plot = xline(MLE_scale(j,i),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
        end
    end

    %plot the posterior distribution
    posteriorPlot = ksdensity(posteriorSample_scale(:,i),'BoundaryCorrection','reflection','Support',[Prior_lowers_scale(i) Prior_uppers_scale(i)]);
    posteriorPlot(1) = posteriorPlot(2);posteriorPlot(end)=posteriorPlot(end-1);
    posterior_plot = plot(Prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),posteriorPlot,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
    
    %axis labels and limits
    xlabel(units(i),'Interpreter','latex','FontSize',10)
    set(gca,'ytick',[])
    xlim([Prior_lowers_scale(i) Prior_uppers_scale(i)]);
    
    %add legend
    if i==2
        legend([prior_plot posterior_plot all_MLE_plot(1) highlight_MLE_plot(1)],...
            {'Prior Density','Posterior Density','Local MLEs','Local MLEs for analysis'},...
            'Location','northoutside')
        
        if n_highlight ==0 %if we don't want to highlight any MLEs
            legend([prior_plot posterior_plot all_MLE_plot(1)],...
                {'Prior Density','Posterior Density','Local MLEs'},...
                'Location','northoutside')
        end

    end
end

%y axis label
ylabel(tiles,'Density (Arbitrary Units)','FontSize',12)

end