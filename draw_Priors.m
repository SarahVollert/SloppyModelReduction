function particles = draw_Priors(N,Priors)
% This function generates a set of N uniform prior, where the upper and
% lower bounds of this distribution are specified by the structure Priors. 

%generate N uniform priors as defined by Priors
particles = zeros(N,length(Priors.prior_lowers));
for i=1:length(Priors.prior_lowers)
    particles(:,i)= rand(N,1).*(Priors.prior_uppers(i) - Priors.prior_lowers(i)) + Priors.prior_lowers(i);
end

end