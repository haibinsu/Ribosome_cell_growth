function [ncog_fit,r_fit, cog_fit]= approach2_Capprox(Mgpep,Rtotal,Actcsol_store, Actncsol_store,Rsol_store,PRcsol_store, PRncsol_store)

Mgseries = 0:0.1:8;
% yourFolder = ['../OPTI-master/Figures/growth_HB/' 'allkmetsameplot_Mg'];
yourFolder = ['../' 'Approximation_Concentration']; %make at the same location as folder containing this script

if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end

%fit for ncog PRnc vs Mg 
[fitresult, gof] = fit(Mgpep, PRncsol_store,'exp1');
coeff = coeffvalues(fitresult);
u = coeff(1);
v = coeff(2);
%u = 0.2;
%v = 0.62;
u/Rtotal
ncog_fit = @(x) u*exp(v*x);

%fit for free ribosome vs Mg 
[fitresult, gof] = fit(Mgpep, Rsol_store,'exp1');
coeff = coeffvalues(fitresult);
p = coeff(1);
q = coeff(2);
p/Rtotal
r_fit = @(x) p*exp(q*x);

cog_fit = @(Rt,x) Rt - (ncog_fit(x) + r_fit(x));

figure
scatter(Mgpep, Rsol_store,'filled')
hold on
plot(Mgseries, r_fit(Mgseries))
xlabel('Free Mg^{2+} (mM)')
ylabel('[R] (\muM)')
legend('Model','Fitting')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'Rfree_Mg'),'png')

figure
scatter(Mgpep, PRcsol_store,'filled')
hold on
plot(Mgseries, cog_fit(Rtotal,Mgseries))
xlabel('Free Mg^{2+} (mM)')
ylabel('[PR_{c}] (\muM)')
legend('Model','Fitting')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'PRc_Mg'),'png')

figure
scatter(Mgpep, PRncsol_store,'filled')
hold on
plot(Mgseries, ncog_fit(Mgseries))
xlabel('Free Mg^{2+} (mM)')
ylabel('[PR_{nc}] (\muM)')
legend('Model','Fitting')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'PRnc_Mg'),'png')

close all
end
