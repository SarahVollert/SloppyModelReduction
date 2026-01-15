function Plots_ModelComparison_PosteriorDensity(Original, NoBatNoTrans)
% This script plots the estimated posterior and prior densities of a sample
% for both the original model and one without the BAT mechanism and without
% transcellular pathway. Note, that part of this function is static and is
% set up for these two models. 

%get info
[~, n_params_original] = size(Original.posterior_sample);
[n_particles_noBatNoTrans, ~] = size(NoBatNoTrans.posterior_sample);

%fill in the missing parameters with 0's
NoBatNoTrans_format.posterior_sample = [zeros(n_particles_noBatNoTrans,1) NoBatNoTrans.posterior_sample(:,1:12) ...
                                        zeros(n_particles_noBatNoTrans,7)   NoBatNoTrans.posterior_sample(:,13)];

%define units for each parameter in the full model
units = define_parameter_units_OriginalModel(n_params_original);

%prepare to rescale
posterior_sample_noBatNoTrans_scale = NoBatNoTrans_format.posterior_sample;
posterior_sample_original_scale = Original.posterior_sample;
prior_sample_scale = Original.prior_sample; 
prior_lowers_scale = Original.Priors.prior_lowers; 
Prior_uppers_scale = Original.Priors.prior_uppers; 

%the rescale order of magnitude for each parameter
scalers = floor(log10(1./Prior_uppers_scale));

%rescaling
for i=1:n_params_original
    posterior_sample_noBatNoTrans_scale(:,i) = posterior_sample_noBatNoTrans_scale(:,i).*10.^scalers(i);
    posterior_sample_original_scale(:,i) = posterior_sample_original_scale(:,i).*10.^scalers(i);
    prior_sample_scale(:,i) = prior_sample_scale(:,i).*10.^scalers(i);
    prior_lowers_scale(i) = prior_lowers_scale(i).*10.^scalers(i);
    Prior_uppers_scale(i) = Prior_uppers_scale(i).*10.^scalers(i);

    if not(scalers(i) ==0)
        units(i) = strcat(units(i),' $\times10^{',int2str(scalers(i)),'}$');
    end
end

%set up figure
figure;
figure_rows = ceil(n_params_original/3);
tiles = tiledlayout(figure_rows,3);
tiles.TileSpacing = 'compact';
tiles.Padding = 'compact';

%for each parameter in the full model
for i=1:n_params_original
    
    %prep figure
    nexttile
    hold on
    
    %plot the prior distribution in grey
    priorPlot = ksdensity(prior_sample_scale(:,i),'BoundaryCorrection','reflection','Support',[prior_lowers_scale(i) Prior_uppers_scale(i)]);  
    priorPlot(1) = priorPlot(2);priorPlot(end)=priorPlot(end-1);
    stepSize = (Prior_uppers_scale(i) - prior_lowers_scale(i))/(length(priorPlot)-1);
    area(prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),priorPlot,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    
    %Plot original model
    k1= ksdensity(posterior_sample_original_scale(:,i),'BoundaryCorrection','reflection','Support',[prior_lowers_scale(i) Prior_uppers_scale(i)]);
    k1(1) = k1(2);k1(end)=k1(end-1);
    plot(prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),k1,'Color','k','LineWidth',1)
    
    %Plot reduced model
    if sum(posterior_sample_noBatNoTrans_scale(:,i)>0)>0
        k4 = ksdensity(posterior_sample_noBatNoTrans_scale(:,i),'BoundaryCorrection','reflection','Support',[prior_lowers_scale(i) Prior_uppers_scale(i)]);
        k4(1) = k4(2);k4(end)=k4(end-1);
        plot(prior_lowers_scale(i):stepSize:Prior_uppers_scale(i),k4,'Color','#ba5eba','LineWidth',1)
    end
    
    %format axis labels and limits
    xlabel(units(i),'Interpreter','latex')
    xlim([prior_lowers_scale(i) Prior_uppers_scale(i)]);
    set(gca,'ytick',[])

    %insert the legend
    if i==2        
        legend('Prior Density','Posterior Density (original model)','Posterior Density (reduced model)','Location','northoutside')
    end
end

%y axis label
ylabel(tiles,'Density (Arbitrary Units)')

end