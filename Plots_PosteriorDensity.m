function Plots_PosteriorDensity(prior_sample, posterior_sample, Priors, units)
% This script plots the estimated posterior and prior densities of a sample.

[~, n_params] = size(posterior_sample);

%locally extract prior bounds
prior_lowers = Priors.prior_lowers;
prior_uppers = Priors.prior_uppers;

%prepare to rescale to the same order of magnitude
posteriorSample_scale = posterior_sample;
priorSample_scale = prior_sample;
Prior_lowers_scale = prior_lowers;
Prior_uppers_scale = prior_uppers;

%the rescale order of magnitude for each parameter
scalers = floor(log10(1./prior_uppers));

%rescaling
for i=1:n_params
    posteriorSample_scale(:,i) = posterior_sample(:,i).*10.^scalers(i);
    priorSample_scale(:,i) = prior_sample(:,i).*10.^scalers(i);
    Prior_lowers_scale(i) = prior_lowers(i).*10.^scalers(i);
    Prior_uppers_scale(i) = prior_uppers(i).*10.^scalers(i);

    if not(scalers(i) ==0)
        units(i) = strcat(units(i),' $\times10^{',int2str(scalers(i)),'}$');
    end
end

%set up figure
figure;
figure_rows = ceil(n_params/5);
tiles = tiledlayout(figure_rows,5);
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
    area(Prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),priorPlot,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);

    %plot the posterior distribution
    posteriorPlot = ksdensity(posteriorSample_scale(:,i),'BoundaryCorrection','reflection','Support',[Prior_lowers_scale(i) Prior_uppers_scale(i)]);
    posteriorPlot(1) = posteriorPlot(2);posteriorPlot(end)=posteriorPlot(end-1);
    plot(Prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),posteriorPlot,'Color','k','LineWidth',1)

    %axis labels and limits
    xlabel(units(i),'Interpreter','latex')
    set(gca,'ytick',[])
    xlim([Prior_lowers_scale(i) Prior_uppers_scale(i)]);
    
    %add legend
    if i==3
        legend('Prior Density','Posterior Density','Location','northoutside');
    end
end

%y axis label
ylabel(tiles,'Density (Arbitrary Units)')

end