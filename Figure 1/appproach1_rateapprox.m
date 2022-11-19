function [fit_kcatKMc,fit_kcatKMnc,fb,fa] = appproach1_rateapprox(Mgpep, kcatKMc, kcatKMnc, qnc, kpepnc, qc, kpepc, khyd)

Mgseries = 0:0.1:8;
% yourFolder = ['../OPTI-master/Figures/growth_HB/' 'allkmetsameplot_Mg'];
yourFolder = ['../' 'Approximation_Rate']; %make at the same location as folder containing this script

if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end


%fit for kcatKMc vs Mg2+ 
ft = fittype( 'a/(b+exp(-c*x))', 'independent', 'x', 'dependent', 'y' );
ft2 = fittype( 'a*x/(b+x)', 'independent', 'x', 'dependent', 'y' );

[fitresult, gof] = fit(Mgpep, kcatKMc,ft);
coeff = coeffvalues(fitresult);
a1 = coeff(1);
b1 = coeff(2);
c1 = coeff(3);

%fit for kcatKMnc vs Mg2+ 
% [fitresult, gof] = fit(Mgpep, kcatKMnc,'exp1');
% coeff = coeffvalues(fitresult);
% u = coeff(1);
% v = coeff(2);

% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,0],...
%                'Upper',[1,1,1],...
%                'StartPoint',[0.01 0.001 0.7]);
[fitresult_nc, gof] = fit(Mgpep, kcatKMnc,ft,'Start',[0.01 0.001 0.7]);
coeff_nc = coeffvalues(fitresult_nc);
a1_nc = coeff_nc(1);
b1_nc = coeff_nc(2);
c1_nc = coeff_nc(3);

[fitresult, gof] = fit(Mgpep, log(kcatKMnc),'poly1');
coeff = coeffvalues(fitresult);
u = coeff(1);
v = coeff(2);



% fit_kcatKMnc = @(x)  u*exp(v*x);

fit_kcatKMc = @(x)  a1./(b1+exp(-c1*x));
fit_kcatKMnc = @(x) exp(u*x+v);

figure
scatter(Mgpep, kcatKMc,46,'filled')
hold on
scatter(Mgpep, kcatKMnc,46,'filled')
plot(Mgseries, fit_kcatKMc(Mgseries),'Color',[0, 0.4470, 0.7410],'LineWidth',1.3)
plot(Mgseries, fit_kcatKMnc(Mgseries),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.3) 
set(gca,'YScale','log')
xlabel('Free Mg^{2+} (mM)')
ylabel('k_{cat}/K_M (\muM^{-1}s^{-1})')
legend('cognate AAA ','near cognate GAA','cognate fit', 'near cognate fit','Location','SouthEast')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'kcatKM'),'png')
% saveas(gca,fullfile(yourFolder,'kcatKM'),'eps')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'kcatKM.eps'),'epsc')

figure
scatter(Mgpep, kcatKMnc,46,'filled')
hold on
plot(Mgseries, fitresult_nc(Mgseries))
plot(Mgseries, fit_kcatKMnc(Mgseries),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.3) 

%compare approximation of kf/(q+kpep) for cog and ncog
%fit for near cognate 
[fitresult, gof] = fit(Mgpep, log(kcatKMnc./(kpepnc+qnc)),'poly1');
coeff = coeffvalues(fitresult);
a1effncs = coeff(1);
a2effncs = coeff(2);

fb = @(x)  exp(a1effncs*x+a2effncs); 
fa = @(x)  1/kpepc*a1./(b1+exp(-c1*x));

figure
scatter(Mgpep, qnc, 'filled')
hold on
plot(Mgseries, fit_kcatKMnc(Mgseries)./fb(Mgseries)-kpepnc)
xlabel('Free [Mg^{2+}] (mM)')
ylabel('k^{rej}_{nc} (aa/(R.s)')
xlim([1 8])
% saveas(gca,fullfile(yourFolder,'krejnc'),'eps')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'krejnc.eps'),'epsc')


figure
scatter(Mgpep, kcatKMc./(qc+kpepc),46,'filled');
hold on
scatter(Mgpep, kcatKMnc./(qnc+kpepnc),46,'filled');
plot(Mgseries, fa(Mgseries),'Color', [0 0.45 0.74]);
plot(Mgseries, fb(Mgseries),'Color',[0.8500, 0.3250, 0.0980]);
set(gca,'YScale','log')
xlabel('Free Mg^{2+} (mM)')
ylabel('Rate ratio k_{f}/(k_{rej}+k_{pep})(\muM^{-1})')
legend('cognate AAA ','near cognate GAA','cognate fit', 'near cognate fit','Location','SouthEast')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'kratio'),'png')

%compare ksynmax and KMeff 
fc = 0.02;
fnc = 0.15;
ksynmax_approx = fc*kcatKMc./(fc*kcatKMc/kpepc+fnc*kcatKMnc./qnc);
% ksynmax_exact = (fc*kcatKMc+fnc*kcatKMnc./(1+qnc/kpepnc))./(fc*kcatKMc/kpepc+fnc*kcatKMnc./qnc);

KMeff_approx = 1./(fc*kcatKMc/kpepc+fnc*kcatKMnc./qnc);
ksynmax_Mg = @(x) fc*fit_kcatKMc(x)./(fc*fit_kcatKMc(x)/kpepc+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));
ksynmax_Mg2 = @(x) (fc*fit_kcatKMc(x).*1.02.*(0.13+x)./(0.2+x))./(fc*fit_kcatKMc(x)/kpepc+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));

% ksynmax_Mg2 = @(x) (fc*fit_kcatKMc(x).*(2-exp(-0.006*x)))./(fc*fit_kcatKMc(x)/kpepc+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));
% ksynmax_Mg2 = @(x) (fc*fit_kcatKMc(x)+fc*exp(-0.002*x).*(1-exp(-1.6*x)))./(fc*fit_kcatKMc(x)/kpepc+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));

KMeff_Mg = @(x) 1./(fc*fit_kcatKMc(x)/kpepc+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));

figure
scatter(Mgpep, ksynmax_approx,'filled')
hold on
plot(Mgseries, ksynmax_Mg(Mgseries))
plot(Mgseries, ksynmax_Mg2(Mgseries),'k')
xlabel('Free Mg^{2+} (mM)')
ylabel('k_{syn}^{max} (aa/(R.s))')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'ksynmax_Mg'),'png')
saveas(gca,fullfile(yourFolder,'ksynmax_Mg.eps'),'epsc')

figure
scatter(Mgpep, KMeff_approx,'filled')
hold on
plot(Mgseries, KMeff_Mg(Mgseries))
plot([1 8], [2 2],'k')
xlabel('Free Mg^{2+} (mM)')
ylabel('K_{M}^{eff} (\muM)')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'KMeff_Mg'),'png')
saveas(gca,fullfile(yourFolder,'KMeff_Mg.eps'),'epsc')


%combine both ksynmax and KMeff together depend on Accuracy
%  try but it doesn't look good because it can't capture the first discrete
%  accuracy point 
% A_calculate = (fc*fit_kcatKMc(Mgseries))./(fnc*fit_kcatKMnc(Mgseries)).*kpepc./(qc+kpepc)./(kpepnc./(fit_kcatKMnc(Mgseries)./fb(Mgseries)-kpepnc+kpepnc));
% 
% figure
% yyaxis left
% % scatter(Mgpep, ksynmax_approx,'filled')
% scatter(Alist, ksynmax_approx,'filled')
% hold on
% % plot(Mgseries, ksynmax_Mg(Mgseries))
% % plot(Mgseries, ksynmax_Mg2(Mgseries),'k')
% plot(A_calculate, ksynmax_Mg(Mgseries))
% xlabel('Accuracy')
% ylabel('k_{syn}^{max} (aa/(R.s))')
% 
% yyaxis right
% scatter(Alist, KMeff_approx,'filled')
% hold on
% plot(A_calculate, KMeff_Mg(Mgseries))
% plot([1 8], [2 2],'k')
% 
% set(gca,'XScale','log')

figure
yyaxis left
scatter(Mgpep, ksynmax_approx,'filled')
hold on
plot(Mgseries, ksynmax_Mg(Mgseries),'-')
plot(Mgseries, ksynmax_Mg2(Mgseries),'k-')
 ylabel('k_{syn}^{max} (s^{-1})')
yyaxis right
scatter(Mgpep, KMeff_approx,'filled')
hold on
plot(Mgseries, KMeff_Mg(Mgseries),'-')
plot([1 8], [2 2],'k-')
xlim([1 8])
xlabel('Free Mg^{2+} (mM)')
ylabel('K_{M}^{eff} (\muM)')
box on
saveas(gca,fullfile(yourFolder,'ksynmax_KMeff_Mg'),'png')
saveas(gca,fullfile(yourFolder,'ksynmax_KMeff_Mg.eps'),'epsc')

%rescale KMeff 
figure
scatter(Mgpep, KMeff_approx/KMeff_approx(3),'filled')
hold on
plot(Mgseries, KMeff_Mg(Mgseries)/KMeff_approx(3))
xlabel('Free Mg^{2+} (mM)')
ylabel('K_{M}^{eff} (mM)')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'KMeff_Mg_rescale'),'png')
saveas(gca,fullfile(yourFolder,'KMeff_Mg_rescale.eps'),'epsc')

%change kpepc = 22 1/s and replot KMeff and ksynmax
KMeff_approxnew = 1./(fc*kcatKMc/22+fnc*kcatKMnc./qnc);
KMeff_Mgnew = @(x) 1./(fc*fit_kcatKMc(x)/22+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));

ksynmax_approxnew = fc*kcatKMc./(fc*kcatKMc/22+fnc*kcatKMnc./qnc);
ksynmax_Mgnew = @(x) fc*fit_kcatKMc(x)./(fc*fit_kcatKMc(x)/22+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));
ksynmax_Mg2new = @(x) (fc*fit_kcatKMc(x).*1.02.*(0.13+x)./(0.2+x))./(fc*fit_kcatKMc(x)/22+fnc*fit_kcatKMnc(x)./(fit_kcatKMnc(x)./fb(x)-kpepnc));

figure
scatter(Mgpep, KMeff_approxnew,'filled')
hold on
plot(Mgseries, KMeff_Mgnew(Mgseries))
plot([1 8], [30 30],'k')
xlabel('Free Mg^{2+} (mM)')
ylabel('K_{M}^{eff} (\muM)')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'KMeff_Mg_22aapersec'),'png')
saveas(gca,fullfile(yourFolder,'KMeff_Mg_22aapersec.eps'),'epsc')

figure
scatter(Mgpep, ksynmax_approxnew,'filled')
hold on
plot(Mgseries, ksynmax_Mgnew(Mgseries))
plot(Mgseries, ksynmax_Mg2new(Mgseries),'k')
xlabel('Free Mg^{2+} (mM)')
ylabel('k_{syn}^{max} (aa/(R.s))')
xlim([1 8])
saveas(gca,fullfile(yourFolder,'ksynmax_Mg_22aapersec'),'png')
saveas(gca,fullfile(yourFolder,'ksynmax_Mg_22aapersec.eps'),'epsc')

close all
end
