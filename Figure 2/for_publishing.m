%Fig 2a
%include coupled nutrient and accuracy effect
%see how ksyneff - accuracy depend on [T3]
%assumption
kpepnc = 0.3; %1/s

khyd = 500; %1/s GTP hydrolysis is fast 
kpepc = 22; %1/s
qc = 1;  %cognate PR rejection rate - assume constant with Mg2+ ,as long as << kpepc 
%assumption taken from 10.1016/j.molcel.2005.12.018
Rtotal = 10; %uM
%total ternary complex = 100 uM 
T3c = 2; %uM 
T3nc = 15; %uM 
%ksyneff expressed in Michaelis - Menten form 
fc = 0.02; %cognate number fraction
fnc = 0.15; %near cognate number fraction
T3 = 100; %uM  

kcatKMc_p1 = [40; 60; 88; 117; 147; 167; 180];  % cognate AAA uM^-1s^-1
kcatKMnc_p1 = [12; 19; 32; 66; 139; 327; 1750]; % near cognate GAA mM^-1s^-1 
kcatKMnc_p1 = kcatKMnc_p1/1000; %convert from mM^-1 to uM^-1
kcatKMpepnc_p1 = [2.018e-4; 3.9e-4; 8.03e-4; 2.7e-3; 9.86e-3; 3.67e-2; 2.5e-1];  %unit is uM^-1s^-1
%near cognate rejection rate constant from Eq. S11
qnc_p1 = (kcatKMnc_p1./kcatKMpepnc_p1-1)*kpepnc; 

%kfc fit based on free energy analysis in SI
kfc_fit = @(x) 23.1./(0.13 + exp(-1.11*x));

%kfnc fit based on free energy analysis in SI
kfnc_fit = @(x) 0.01*exp(0.7*x);

mglist = [0.7, 1.3, 1.8, 2.3, 3.4, 4.6, 7.5];

mgseries = 0:0.01:8;

%Fig. S9a
figure
scatter(mglist, kcatKMc_p1)
hold on
scatter([mglist(1), mglist(3)], [kcatKMc_p1(1), kcatKMc_p1(3)],'x','MarkerEdgeColor',[0, 0.4470, 0.7410])
plot(mgseries, kfc_fit(mgseries),'k')
yticks([20, 40, 60, 80, 100, 120, 140, 160, 180])
set(gca, 'YScale', 'log')
legend('Johansson et al. 2012, Zhang et al. 2016', 'Extrapolation','Theoretical line')
xlabel('[Mg^{2+}] (mM)') 
ylabel('k^c_f (\muM^{-1}s^{-1})')
ylim([20 190])
box on
set(gca,'FontName','Helvetica','FontSize',16)



%Fig. S9b
figure
scatter(mglist, kcatKMnc_p1)
hold on
scatter([mglist(1), mglist(3)], [kcatKMnc_p1(1), kcatKMnc_p1(3)],'x','MarkerEdgeColor',[0, 0.4470, 0.7410])
plot(mgseries, kfnc_fit(mgseries),'k')
set(gca, 'YScale', 'log')
% ylim([1e-3 1])
legend('Johansson et al. 2012, Zhang et al. 2016', 'Extrapolation', 'Theoretical line')
xlabel('[Mg^{2+}] (mM)')
ylabel('k^{nc}_f (\muM^{-1}s^{-1})')
box on
set(gca,'FontName','Helvetica','FontSize',16)

%Fig. S9c
figure
scatter(mglist, kcatKMpepnc_p1)
hold on
scatter([mglist(1), mglist(3)], [kcatKMpepnc_p1(1), kcatKMpepnc_p1(3)],'x','MarkerEdgeColor',[0, 0.4470, 0.7410])
set(gca, 'YScale', 'log')
legend('Johansson et al. 2012, Zhang et al. 2016', 'Extrapolation')
xlabel('[Mg^{2+}] (mM)')
ylabel('(k_{cat}/K_M)^{nc}_{pep} (\muM^{-1}s^{-1})')
box on
set(gca,'FontName','Helvetica','FontSize',16)

krejncfit = @(x) 15.2*exp(-0.3*x);

%Fig. S9d
figure
scatter(mglist, qnc_p1)
hold on
scatter([mglist(1), mglist(3)], [qnc_p1(1), qnc_p1(3)],'x','MarkerEdgeColor',[0, 0.4470, 0.7410])
plot(mgseries, krejncfit(mgseries), 'k')
legend('Johansson et al. 2012, Zhang et al. 2016', 'Extrapolation', 'Theoretical line')
xlabel('[Mg^{2+}] (mM)')
ylabel('k^{nc}_{rej} (s^{-1})')
box on
set(gca,'FontName','Helvetica','FontSize',16)

ksynmax = @(kfc,kpepc,kfnc,krejnc) (fc*kfc+fnc*kfnc./(1+krejnc./kpepnc))./(fc*kfc/kpepc+fnc*kfnc./krejnc);
KMeff = @(kfc,kpepc,kfnc,krejnc) 1./(fc*kfc/kpepc+fnc*kfnc./krejnc);
Accuracy = @(kfc,kpepc,kfnc,krejnc,krejc) fc*kfc./(fnc*kfnc).*(kpepc/(krejc+kpepc))./(kpepnc./(krejnc+kpepnc));
speed_f = @(kfc,kpepc,kfnc,krejnc,T3) ksynmax(kfc,kpepc,kfnc,krejnc)*T3./(KMeff(kfc,kpepc,kfnc,krejnc)+T3);

% ksynmax 
ksynmaxlist_p1 = ksynmax(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1);
% KMeff
KMefflist_p1 = KMeff(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1);
% Accuracy
Alist_p1 = Accuracy(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1,qc);

%Fig. 2b 
T3_100 = 100; %uM of ternary complex 
figure
yyaxis left
scatter(Alist_p1, ksynmaxlist_p1,'o','MarkerEdgeColor',[0, 0.4470, 0.7410])
hold on
scatter([Alist_p1(1) Alist_p1(3)], [ksynmaxlist_p1(1) ksynmaxlist_p1(3)],'x','MarkerEdgeColor',[0, 0.4470, 0.7410])
% scatter([Alist_p1(1) Alist_p1(3)], [ksynmaxlist_p1(1) ksynmaxlist_p1(3)],'o','MarkerEdgeColor',[0, 0.4470, 0.7410], 'LineWidth',1)
% s.AlphaData = [0.1 0.1];
% s.MarkerFaceAlpha = 'flat';
% s.MarkerEdgeColor = 'none';
plot(Alist_p1, speed_f(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1,T3_100), '^-', 'MarkerFaceColor',[0, 0.4470, 0.7410])
xlabel('Accuracy')
ylabel('Speed (1/s)')
legend('k^{cell}_{synmax}', 'k^{cell}_{syneff} at [T_3] = 100 \muM')
ylim([10 24])

yyaxis right
scatter(Alist_p1, KMefflist_p1, 40, 'sq','MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
scatter([Alist_p1(1) Alist_p1(3)], [KMefflist_p1(1) KMefflist_p1(3)], 40, 'x','MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
% scatter([Alist_p1(1) Alist_p1(2)], [KMefflist_p1(1) KMefflist_p1(2)],'sq','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'LineWidth',1)
% s2.AlphaData = [0.1 0.1];
% s2.MarkerFaceAlpha = 'flat';
% s2.MarkerEdgeColor = 'none';
ylabel('K_{M}^{eff} (\muM)')
set(gca,'XScale','log')
ylim([0 40])
yticks([0 : 5: 40])
xlim([10 1e5])
set(gca,'FontName','Helvetica','FontSize',16)


%% relate competition effect on growth rate

run('bionumers.m')

close all

% ksynmax = 20*3600; %aa/h
% Ka = 0.25e-3; %M
% Ke = 1e-2; %M %5e-2 1e-2
Ke = 2.5e-3;
Jmtm = 5.2e4/6 ; %convert ATP to aa because 6 ATP required for 1 peptide bond
%https://doi.org/10.1016/j.cels.2019.06.003

%use KMeff from part 1 for ternary complex instead of rescaling 
k2 = 100*3600; %1/h max charging speed 
% k2 = 80*3600;
R = 20; %uM 
% S = R/1.5; % --> charging synthase/ribosome ratio = 1/1.5 
S = R/(5); %S/R = 1/5 is more important 

c_factor = 0.0364; %constant factor between T3 and aa

%cell mass (aa)
M = @(y) NR*y(:,1) + NP*y(:,2) + y(:,3);

%volume (um^3)
V = @(y) M(y)/rho; %um^3 

%R mass fraction
Rmf = @(y) NR*y(:,1)./M(y);

%P mass fraction
Pmf = @(y) NP*y(:,2)./M(y);

%aa mass fraction
aamf = @(y) 1 - Rmf(y) - Pmf(y);

%% explicit simulation instead of putting in function

nSim = 100;
uRlist = 0.01:0.01:0.9;

%Michaelis constant for amino acid (not use)
% Kaselect = 2.5e-4*ones(length(ksynmaxlist_p1),1); %unit is M

%metabolic rate depends on nutrient 
kmetsample = [0.08 0.13 0.18 0.25 0.3]'*3600;   %1/h


%values for fitting mass vs growth rate to match experiment
Mx = 1.65e9;
My = 2e8;
M0f = @(x) My + Mx*x;   %x is growth rate (1/h)

%mass list with unit is amino acid (aa)
Mpicklist = 1e9*[1.18,    1.89,    2.3683,    3.0574,    3.4404];

%to indicate different nutrient quality
color_list = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

%the updated one that is consistent with part 1 rate constants vs Mg
%for modultion factor m 
blist = [18000, 15000, 5000, 1000, 300];
clist = [2.5e-4, 2.3e-4, 5e-4, 8e-4, 3e-3];
T3_A = @(x,b,c) 1.2./(2+exp(-c*(x-b)))+0.4;

%dummy - not important
qA = @(x) exp(-0.01*NR./x);  %x is accuracy 


% three times higher than max value in Table IV
% https://pubmed.ncbi.nlm.nih.gov/13271454/
% protein degradation rate constant
degrade = @(x) 0.03; %1/h 

odestore_kmet = cell(length(kmetsample),1);
odemax_kmet = cell(length(kmetsample),1);
kelongstore_kmet = cell(length(kmetsample),1);
korginalstore_kmet = cell(length(kmetsample),1);
Alist_kmet = cell(length(kmetsample),1);
lambda_ana_kmet = cell(length(kmetsample),1); 
koriginal_kmet = cell(length(kmetsample),1);
ksyn_kmet = cell(length(kmetsample),1);
flux_kmet = cell(length(kmetsample),1);
%store row index that gives max growth rate 
%among 5 nutrient conditions, each condition picks ksynmax-KMeff option
%that give fastest growth rate 
id_store = NaN*ones(length(kmetsample),1);  

%store the fastest growth rate and corresponding accuracy in each nutrient condition 
grmax_nutrient = NaN*ones(length(kmetsample),1);
ksynmax_nutrient = NaN*ones(length(kmetsample),1);
KM_nutrient = NaN*ones(length(kmetsample),1);
kelong_nutrient = NaN*ones(length(kmetsample),1);
koriginal_nutrient =  NaN*ones(length(kmetsample),1);
A_nutrient = NaN*ones(length(kmetsample),1);
flux_nutrient = NaN*ones(length(kmetsample),3);
ksyn_ana_nutrient = NaN*ones(length(kmetsample),1);
gr_ana_nutrient = NaN*ones(length(kmetsample),1);

%1st - 2nd - 3rd col: input flux - synthesis flux - degradation flux 
cost_nutrient = NaN*ones(length(kmetsample),1);
chosen = NaN*ones(length(kmetsample),11);

%step 1 - choose a nutrient condition by fixing kmet
%step 2 - choose accuracy, ksynmax, KMeff (i.e. fixing Mg2+ concentration)
%step 3 - vary ribosomal mass allocation to find fastest growth rate
%step 4 - store fastest growth rate across different Mg2+ concentration
%step 5 - select the fastest among the fastest 
for m = 1 : length(kmetsample)  % step 1

    kmet = kmetsample(m);
    Mth = Mpicklist(m);
    b = blist(m);
    c = clist(m);
    
    odestore = cell(length(ksynmaxlist_p1),1);
    odemax = NaN*ones(length(ksynmaxlist_p1),11);
    % 1st - 2nd - 3rd col: number of ribosomes R - metabolic proteins P - amino acid aa
    % 4th - 5th - 6th col: growth rate (1/h) - synthesis allocation for ribosomes uR - uP
    % 7th - 8th - 9th col: mass fraction of R - P - aa
    % 10th - 11th: number of ternary complex - cell mass
    kelongstore = NaN*ones(length(ksynmaxlist_p1),1);
    koriginal = NaN*ones(length(ksynmaxlist_p1),1);

    fluxstore = NaN*ones(length(ksynmaxlist_p1),1);
    %1st col: input flux, 2nd col: synthesis flux, 3rd col: degrade flux
    %each row corresponds to each ksynmax-KMeff-Accuracy
   
    for j = 1 : length(ksynmaxlist_p1)  %step 2 
        ksynmaxT3 = ksynmaxlist_p1(j)*3600; %1/h
        KM = KMefflist_p1(j)*1e-6; %M

        KMeff = KM;
%         Ka = Kaselect(j); %M 
        A = Alist_p1(j); %Accuracy from Fig. 2b 

        %cell mass (aa)
        Mf = @(y) NR*y(1)+NP*y(2)+y(3);
        %volume (um^3)
        Vf = @(y) Mf(y)/rho; %um^3 
        %R mass fraction
        Rmfc = @(y) NR*y(1)./Mf(y);

        %P mass fraction
        Pmfc = @(y) NP*y(2)./Mf(y);

        %aa mass fraction
        aamfc = @(y) 1 - Rmfc(y) - Pmfc(y);

        %full mapping [aa] to [T3]
        %T3convertf_T3 gives [T3] concentration from [aa] 
        T3convertf_T3 = @(y) c_factor*y(3);

        ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

        flist = {Mf,Vf,ktransf_T3,T3convertf_T3};

        %metabolic rate
        kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
        flist{5} = kmetf; 
        flist{6} = qA; 
        flist{7} = T3_A;
        flist{8} = degrade;
        
        store = NaN*ones(length(uRlist),10);
%         Mth = Mpick;
        
        parfor i =  1 : length(uRlist)  % step 3

            R0 = 2e4;
            P0 = 1e5;
            M0 = 3e8; %aa
            aa0 = M0 - (NR*R0 + P0*NP);
   
            y0 = [R0 P0 aa0];
            par = [ksynmaxT3 kmet uRlist(i) NR NP A b c];
            k = 1;
            t_tot = 0;
            y_tot = y0;
            tspan = [0 200];  
            gr = 1;
            Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));
            %run to reach exponential steady state 
            while k <= nSim
                [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),tspan,y0, Opt);
                if isempty(ye) ~= 1
                    t_tot = [t_tot ; [t;te] + t_tot(end)];
                    y_tot = [ y_tot ; [y;ye] ];    
                    y0 = ye/2;
                    k = k + 1;           
                else
                    gr = NaN;
                    break;
                end
            end
        %     semilogy(t_tot,y_tot)
            if isnan(gr) ~= 1 
            [t,y,te,ye,ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),[0:1e-3:2*te],y0,Opt);
        %     semilogy(t,y)
            %fitting to get growth rate 
            hfit = dy_dt_eval(t,par,flist,y,M);
            
            store(i,:) = [y0 hfit.p1 uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];

            end
         
        end
    %     plot(uRlist,store(:,4))
        odestore{j} = store;
        
        % step 4
        %pick uR parition that gives fastest growth rate 
        [val, idx] = max(store(:,4));
        odemax(j,:) = [store(idx,:) M(store(idx,:))];

        kelongstore(j) = ktransf_T3(odemax(j,:))*T3_A(A,b,c);
        koriginal(j) = ktransf_T3(odemax(j,:));
        
        fluxstore(j,1) = odemax(j,2)*kmetf(odemax(j,:));
        fluxstore(j,2) = odemax(j,1)*kelongstore(j);
        fluxstore(j,3) = degrade(A)*(NR*odemax(j,1)+NP*odemax(j,2));
    end 

    %analytical max growth rate for each Mg2+ (i.e. ksynmax, KMeff, ...)
    %concentration of amino acid 
    % Eq. S62
    aa_ana = @(ksynmax,KMeff,A,b,c) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*T3_A(A,b,c)*NP)).^(1/3);
    % Eq. S63
    gr_ana = @(ksynmax,KMeff,A,b,c)  aa_ana(ksynmax,KMeff,A,b,c)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*T3_A(A,b,c))+aa_ana(ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*T3_A(A,b,c))));
    % Eq. S64
    ksyn_function = @(ksynmax,KMeff,A,b,c)  3/2*NR*gr_ana(ksynmax,KMeff,A,b,c)./(1 - gr_ana(ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*T3_A(A,b,c))) + 3/2*NR.*gr_ana(ksynmax,KMeff,A,b,c)./ksynmax); 

    %list of solutions 
    lambda_ana = gr_ana(ksynmaxlist_p1*3600, KMefflist_p1, Alist_p1, b, c);
    ksyn_ana = ksyn_function(ksynmaxlist_p1*3600, KMefflist_p1, Alist_p1, b, c);

    odestore_kmet{m} = odestore; 
    odemax_kmet{m} = odemax; 
    kelongstore_kmet{m} = kelongstore;
    korginalstore_kmet{m} = koriginal;
    flux_kmet{m} = fluxstore;
    Alist_kmet{m} = Alist_p1;
    lambda_ana_kmet{m} = lambda_ana;
    ksyn_kmet{m} = ksyn_ana;
    
    %store fastest growth rate and corresponding quantities for each nutrient condition
    %due to variation in accuracy, ksynmax and KMeff

    %step 5
    %select the fastest growth rate among the fastesst
    [val, id] = max(odemax(:,4));  %max growth rate
    id_store(m) = id;
    grmax_nutrient(m) = val;
    ksynmax_nutrient(m) = ksynmaxlist_p1(id);
    KM_nutrient(m) = KMefflist_p1(id);

    kelong_nutrient(m) = kelongstore(id);
    koriginal_nutrient(m) = koriginal(id);
    A_nutrient(m) = Alist_p1(id);
    
    ksyn_ana_nutrient(m) = ksyn_ana(id);
    gr_ana_nutrient(m) = lambda_ana(id);
    
    flux_nutrient(m,:) = fluxstore(id,:);
    
    %the fastest among the fastest growth rate is stored here
    chosen(m,:) = odemax(id,:);
end


%check analytical growth rate vs numerical growth rate 
gr_ana_nutrient./chosen(:,4)
%check analytical synthesis speed vs numerical synthesis speed 
ksyn_ana_nutrient./kelong_nutrient

%c_factor to convert aa to Ternary complex
kelongf_all = @(kmax,KM,y,Alist,b,c) T3_A(Alist,b,c).*3600.*kmax.*c_factor.*y(:,3)./(c_factor*y(:,3)+KM.*1e-6.*NA.*V(y)*1e-15);
%metabolic rate
kmetf_all = @(y,kmet) kmet./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
%input flux
JM_f = @(y,kmet) y(:,2).*kmetf_all(y,kmet);
%output flux
JS_f = @(y,kmax,KM,Alist,b,c) y(:,1).*kelongf_all(kmax,KM,y,Alist,b,c);
%degradation flux
Jd_f = @(y,Alist) degrade(Alist).*(NR*y(:,1)+NP*y(:,2));

%check flux balances
ratio_result = cell(length(kmetsample),1);
for m = 1 : length(kmetsample)
    % m = 4;
    temp = NaN*ones(length(flux_kmet{m}(:,1)),5);
    % JM - degradation flux = lambda * M
    temp(:,1) = flux_kmet{m}(:,1)./(flux_kmet{m}(:,3)+odemax_kmet{m}(:,4).*odemax_kmet{m}(:,11));
    % R*kelong + lambda*aa = JM
    temp(:,2) = flux_kmet{m}(:,1)./(flux_kmet{m}(:,2)+odemax_kmet{m}(:,4).*odemax_kmet{m}(:,3));

    %input flux check: numerical vs analytical 
    temp(:,3) = flux_kmet{m}(:,1)./JM_f(odemax_kmet{m},kmetsample(m));
    %synthesis flux check
    temp(:,4) = flux_kmet{m}(:,2)./JS_f(odemax_kmet{m},ksynmaxlist_p1,KMefflist_p1,Alist_p1,blist(m),clist(m));
    %degradation flux check
    temp(:,5) = flux_kmet{m}(:,3)./Jd_f(odemax_kmet{m},Alist_p1);
    
    ratio_result{m} = temp;
end

%growth rate vs accuracy
%Fig 2c
figure
for i = 1 : length(kmetsample)
   plot(Alist_p1, odemax_kmet{i}(:,4), 'o-', 'Color', color_list(i,:), 'MarkerSize', 4) 
   hold on
   scatter(A_nutrient(i),chosen(i,4),50,color_list(i,:),'filled')
end
set(gca,'XScale','log')
yticks([0.4:0.4:2])
ylim([0.4 2])
xlim([10 1e5])
xlabel('Accuracy A')
ylabel('\lambda^* (1/h)')
set(gca,'FontName','Helvetica','FontSize',16)

%% Fig 2d
% compare number of ribosomes - growth rate with experiments 

%to get #R/cell vs growth rate and JM vs growth rate 
%and gr/JM vs gr 
%in Bremmer and Dennis  10.1128/ecosal.5.2.3
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h
RNApercell = [23 44 76 128 180 214]; %ug/10^9 cells
Rpercell = [8 15 26 44 61 73]*1e3; 
kselect = [13 18 21 22 22 22]*3600; %aa/(R.h)

%data from SI of  10.1038/nmicrobiol.2016.231 - Hwa data
load('dataextract')
%1st - 3rd col: #R/cell - #tRNA/cell - growth rate (1/h)
%4th col: peptide synthesis speed (aa/(R*s)) estimated from data_Dai.mat 
%5th col: R + P dry mass (pg) estimated from Liudata.mat 

%mass fraction growth law
figure
for j =  1 : length(kmetsample)
   datatemp = odestore_kmet{j}(id_store(j));
   Jmtemp = datatemp{1,1}(:,2).*kmetsample(j)./(1+(datatemp{1,1}(:,3)./(Ke*V(datatemp{1,1})*1e-15*NA)).^2);
   [val_check, id_check] = max(datatemp{1,1}(:,4)); %max growth rate
   if val_check == chosen(j,4)
        frange = find(abs(Jmtemp-Jmtemp(id_check))<=0.1*Jmtemp(id_check)); 
   end
   plot(datatemp{1,1}(frange,7),datatemp{1,1}(frange,4),'Color',color_list(j,:))
   hold on
end
%fitting numerical data to a straight line
mfit = fit(chosen(:,7),chosen(:,4),'poly1')
plot((0.1:0.1:0.5),mfit.p1*(0.1:0.1:0.5)+mfit.p2,'Color','k')
hold on
scatter(chosen(:,7),chosen(:,4),50,color_list,'filled')
scatter(glawdata(1:6,2),glawdata(1:6,1),'k','diamond','LineWidth',1.25) %Hwa
scatter(glawdata(7:11,2),glawdata(7:11,1),50,'k','<','LineWidth',1.25)  %Forchhammer & Lindahl
scatter(glawdata(12:end,2),glawdata(12:end,1),'k','o','LineWidth',1.25) %Bremer & Dennis
ylabel('\lambda^* (1/h)')
xlabel('\phi_R^*')
xlim([0.1 0.5])
ylim([0 2])
yticks([0:0.4:2])
xticks([0.1:0.1:0.5])
set(gca,'FontName','Helvetica','FontSize',16)
legend('','','','','','Theoretical line','Numerical results','Scott et al. 2010','Bremer & Dennis 2009','Forchhammer & Lindahl 1971','Location','SouthEast')

close 

%R vs gr to illustrate R is proportional to gr^2
%fit ribosomal mass fraction vs growth rate 
hfit = fit(chosen(:,4),chosen(:,7),'poly1')
% plot(hfit,chosen(:,4),chosen(:,7))
% hold on
% scatter(glawdata(1:6,1),glawdata(1:6,2),'k','diamond','LineWidth',1.25) %Hwa
% scatter(glawdata(7:11,1),glawdata(7:11,2),50,'k','<','LineWidth',1.25)  %Forchhammer & Lindahl
% scatter(glawdata(12:end,1),glawdata(12:end,2),'k','o','LineWidth',1.25) %Bremer & Dennis

% --> phiR = 0.13*gr + 0.15
% since M = 1.65e9*gr + 2e8
% --> R ~ 1/NR*(214500000*gr^2 + 273500000*gr + 30000000) 
%number of ribosomes vs growth rate
fR = @(x) 1/NR*M0f(x).*(hfit.p1*x + hfit.p2);



combined_data = NaN*ones(11,2);
%number of ribosome
combined_data(:,1) = [dataextract(:,1); Rpercell'];
%growth rate
combined_data(:,2) = [dataextract(:,3); gr'];
g = fittype('a*x^2', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a'});
gfit_expdata = fit(combined_data(:,2), combined_data(:,1), g);

R_fit = fit((0:0.1:2)', fR(0:0.1:2)', g);
%approximate R = a/ksynmax*gr^2
R_approx = @(x) R_fit.a*x.^2;
R_approx2 = @(x) 2.7166e4*x.^2; %obtained from a_ave in FIg. 2e section

figure
plot(fR(0:0.1:2),0:0.1:2,'k')
hold on
scatter(chosen(:,1), chosen(:,4),50,color_list,'filled')
scatter(dataextract(:,1), dataextract(:,3),'sq','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1)
scatter(Rpercell, gr, 50,'x','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1)
% plot(R_approx(0:0.1:2), 0:0.1:2)
% plot(R_approx2(0:0.1:2), 0:0.1:2, 'r')

ylabel('\lambda^* (1/h)')
xlabel('R^* (\times10^4 molecules/cell)')
yticks(0:0.4:2) 
ylim([-1.3 2.8])
xlim([0 13e4])
xticks([0:2e4:12e4]) 
legend('Theoretical results','Numerical results','Dai et al. 2016','Bremer & Dennis 2008','Location','NorthWest')
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)

p = get(gca, 'Position');
h = axes('Parent', gcf, 'Position', [p(1)+0.4 p(2)+0.05 0.45 0.33]);

for j =  1 : length(kmetsample)
   datatemp = odestore_kmet{j}(id_store(j));
   Jmtemp = datatemp{1,1}(:,2).*kmetsample(j)./(1+(datatemp{1,1}(:,3)./(Ke*V(datatemp{1,1})*1e-15*NA)).^2);
   [val_check, id_check] = max(datatemp{1,1}(:,4)); %max growth rate
   if val_check == chosen(j,4)
        frange = find(abs(Jmtemp-Jmtemp(id_check))<=0.1*Jmtemp(id_check)); 
   end
   plot(datatemp{1,1}(frange,7),datatemp{1,1}(frange,4),'Color',color_list(j,:))
   hold on
end
mfit = fit(chosen(:,7),chosen(:,4),'poly1')
plot((0.1:0.1:0.5),mfit.p1*(0.1:0.1:0.5)+mfit.p2,'Color','k')
hold on
scatter(glawdata(1:6,2),glawdata(1:6,1),25,'k','diamond','LineWidth',1.1) %Hwa
scatter(glawdata(7:11,2),glawdata(7:11,1),25,'k','<','LineWidth',1.1)  %Forchhammer & Lindahl
scatter(glawdata(12:end,2),glawdata(12:end,1),25,'k','o','LineWidth',1.1) %Bremer & Dennis
scatter(chosen(:,7),chosen(:,4),50,color_list,'filled')
ylabel('\lambda^* (1/h)')
xlabel('\phi_R^*')
xlim([0.1 0.6])
ylim([0 2])
yticks([0:0.4:2])
xticks([0.1:0.1:0.5])
set(gca,'FontName','Helvetica','FontSize',12)
legend('','','','','','Theoretical line','Scott et al. 2010','Bremer & Dennis 2009','Forchhammer & Lindahl 1971','Numerical results','Location','SouthEast')

%% Fig 2e

% compare the flux with metabolic fluxes from PNAS paper
% DeLong et al., PNAS, 2010
mW = readtable('data_scale.xlsx');
mW = table2array(mW);

fJMapprox2 = @(x) x.*M0f(x);
JM = chosen(:,2).*kmetsample./(1+(chosen(:,3)./(Ke*V(chosen)*1e-15*NA)).^2);
JS = chosen(:,1).*kelong_nutrient;

%fit JM vs M^2
JM_fit = fit(chosen(:,11), JM, g)
JM_approx = @(x) JM_fit.a*x.^2;
JM_approx2 = @(x) 1/(R_fit.a*ksynmax)*x.^2;

a_ave = (R_fit.a*ksynmax+1/JM_fit.a)/2;
JM_approx3 = @(x) 1/(a_ave)*x.^2;


f = figure
f.Position = [360,141,652,477]
plot(M0f(0.2:0.1:2.2),fJMapprox2(0.2:0.1:2.2),'k')
hold on
% plot(1e9*[0:0.01:4], JM_approx(1e9*[0:0.01:4]))
% plot(1e9*[0:0.01:4], JM_approx2(1e9*[0:0.01:4]))
% plot(1e9*[0:0.01:4], JM_approx3(1e9*[0:0.01:4]))

% plot(M0f(0.2:0.1:2.2),fJMapprox2(0.2:0.1:2.2),'k')
scatter(Mpercell*0.6,Rpercell.*kselect,50,'kx','LineWidth',1.25) %Bremer - Dennis 
scatter(dataextract(:,5)*0.6*1e-12/maa,dataextract(:,1).*dataextract(:,4)*3600,50,'k^','LineWidth',1.25) %Hwa-Liu
scatter(chosen(:,11), JM, 50, color_list, 'filled')
% scatter(chosen(:,11), JS)
xlabel('M^* (aa)')
ylabel('J_M^* (aa/h)')
xlim([0 9e9])
ylim([0 17e9])
yticks([0:2.5:15]*1e9)
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)
legend('Theory line','Bremer & Dennis 2008','Zheng et al. 2020, Dai et al. 2016','Numerical results', 'Location', 'SouthEast')

% axes('Position',[.1 .4 .5 .6])
p = get(gca, 'Position');
h = axes('Parent', gcf, 'Position', [p(1)+0.08 p(2)+0.45 0.28 0.25]);

plot(0.4:0.1:2,M0f(0.4:0.1:2),'LineWidth',1.2,'Color','k')
hold on
% scatter(chosen(:,4),chosen(:,11), 250, color_list, 'filled')
% scatter(gr,Mpercell*0.6,290,'kx','LineWidth',1.5)  
scatter(chosen(:,4),chosen(:,11), 50, color_list, 'filled')
scatter(gr,Mpercell*0.6,90,'kx','LineWidth',1.5)  
xlabel('\lambda^* (1/h)')
ylabel('M^* (aa)')
set(gca,'FontSize',12)
xticks(0:0.4:2)
yticks(0:1e9:4e9)
ylim([0 4e9])
xlim([0 2.2])
box on
set(gca,'FontName','Helvetica')

h = axes('Parent', gcf, 'Position', [p(1)+0.45 p(2)+0.45 0.25 0.25]);

scatter(mW(:,1)*dwconv*RPpercent/maa,mW(:,2)*nATP*3600/nconv,[],[169 169 169]/255,'sq','filled')
hold on
scatter(Mpicklist, flux_nutrient(:,1),50,color_list,'filled')
scatter(mW(32,1)*dwconv*RPpercent/maa,mW(32,2)*nATP*3600/nconv,50,[255 153 85]/255,'filled') %Ecoli
set(gca,'XScale','log','YScale','log')
% xlabel('M^* (aa)')
% ylabel('J_M^* (aa/h)')
xlabel('log(M^*)')
ylabel('log(J_M^*)')
xlim([1e7 1e11])
ylim([1e4 1e12])
yticks([1e4 1e6 1e8 1e10 1e12])
xticks([1e7 1e8 1e9 1e10 1e11])
% legend('','E Coli')
set(gca,'FontName','Helvetica','FontSize',12)
box on

% annotation('textarrow',[0.04 0.024],[0.85 0.85] ,'String','Optimal intake flux')
ticks = get(gca,'XTickLabel')
expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
replace = '$1'; 
exponents_X = regexprep(ticks,expression,replace);
ticks = get(gca,'YTickLabel')
exponents_Y = regexprep(ticks,expression,replace);
set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)

%% Figures in SI

%Fig. S6a
figure
load('Lidata')
%from https://doi.org/10.1038/s41564-020-0717-x
%1st - 2nd : growth rate (1/h) - mass per cell (OD600 ml per 1e9 cells)
id = find(Lidata(:,1) >= 0);
[mfit,gof] = fit(Lidata(id,1),Lidata(id,2),'poly1');

scatter(Lidata(:,1),Lidata(:,2),'filled')
hold on
plot(0.4:0.01:1.8,mfit(0.4:0.01:1.8))
xlabel('Growth rate (1/h)')
ylabel('Average cell mass (OD600 ml per 10^9 cells)')

load('growthdata')
% 10.15252/msb.20156178

%1st - 4th = growth rate, RNA mass/OD (ug), total protein mass/OD (ug), #cell/OD600 
%5th col: cell volume (um^3)
[~,id] = sort(growthdata(:,1)); %sort growth rate from slow to fast 
growthdata = growthdata(id,:);
% [~,id] = sort(growthCm(:,1)); %sort growth rate from slow to fast 
% growthCm = growthCm(id,:);

%compare protein density
%Fig. S6b
figure
plot(growthdata(:,1),growthdata(:,3)./(growthdata(:,4).*growthdata(:,5))/1e-8,'sq-')
ylim([8 18])
ylabel('Total protein density (10^{-8} \mug/\mum^3)')
xlabel('Growth rate (1/h)')
xlim([0.2 2])


%plot speed modulation factor m(A) vs accuracy at different nutrient quality 
%Fig. S7a
figure
for i = 1 : length(kmetsample)
    plot((10:1e5), T3_A((10:1e5),blist(i),clist(i)), 'Color', color_list(i,:))
    hold on    
end
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('m(A)')
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)

%gr/JM vs gr 
%Fig. S7b
figure
plot((0.2:0.1:2.5),1./M0f(0.2:0.1:2.5)*1e10)
hold on
scatter(chosen(:,4), chosen(:,4)./JM*1e10, 50, color_list, 'filled')
scatter(gr,gr./(Rpercell.*kselect)*1e10,50,'kx','LineWidth',1.25)
scatter(dataextract(:,3),dataextract(:,3)./(dataextract(:,1).*dataextract(:,4)*3600)*1e10,50,'ksq','LineWidth',1.25)
legend('Theory line','Numerical results','Bremer & Dennis 2008','Zheng et al. 2020, Dai et al. 2016')
xlim([0 2.1])
ylim([2e-10 16e-10]*1e10)
xticks(0:0.4:2)
set(gca,'FontName','Helvetica','FontSize',16)
ylabel('\lambda^*/J_M^* (\times 10^{-10} 1/aa)')
xlabel('\lambda^* (1/h)')

%Fig. S7c
%elongation speed vs growth rate from Hwa, PNAS, 2013 
speed_grexp = readtable('speed_grexp.xlsx');
speed_grexp = table2array(speed_grexp);
%1st col: growth rate (1/h)
%2nd col: translation speed (aa/s)

ft = fittype( '(a*x+c)/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [10 0.22 1];
opts.Upper = [22 1 5];
opts.Robust = 'Bisquare';
opts.StartPoint = [10 0.3 4];
[fitresult, gof] = fit(chosen(:,4),1/3600*kelong_nutrient, ft, opts );
ffit = @(x) 3600*(fitresult.a*x+fitresult.c)./(fitresult.b+x);
%based on 2 scaling relations --> both are approximation 
ftheoryapprox_fit = @(x) (NR/(hfit.p1)*x+NR*My/(Mx*hfit.p1))./((hfit.p1*My+hfit.p2*Mx)/(Mx*hfit.p1)+x);
% ksyn = gr*M/R
%fR = 0.1339*gr + 0.1462
%
ftheoryapprox = @(x) x.*M0f(x)./fR(x);


%based on growth rate optimization --> actual theory 
%x is gr 
KBlist = (3*NR./(2*ksynmax_nutrient*3600).*KM_nutrient/c_factor*1e-6./(3*NR./(2*ksynmax_nutrient*3600).*KM_nutrient/c_factor*1e-6 + (NR*KM_nutrient/c_factor.*1e-6.*kmetsample.*Ke^2./(2*ksynmax_nutrient*3600*NP)).^(1/3).*(NP./kmetsample+NR./(ksynmax_nutrient*3600))))./(3*NR./(2*(ksynmax_nutrient*3600)));
KB = mean(KBlist);
% KB = (3*NR/(2*ksynmax)*KMave_chosen/(3*NR/(2*ksynmax)*KMave_chosen + ()^(1/3)*(NP/kmetave+NR/ksynmax)))/(3*NR/(2*ksynmax));
ksynmax_set = 22*3600; %1/h
KC = 0.56; %1/h
ksyn_gr = @(x)  ksynmax_set*x./(KC + x); 

%Fig. S7c
figure
% figure
% plot(0:0.1:2,ffit(0:0.1:2)/3600)
% hold on
scatter(speed_grexp(:,1), speed_grexp(:,2))
hold on 
%quadratic approx line
plot(0.4:0.01:2,ftheoryapprox(0.4:0.01:2)/3600,'k')
%theoretical line
plot((0.4:0.1:2), ksyn_gr(0.4:0.1:2)/3600,'k--')
% hold on
% plot(0.4:0.01:2,ftheoryapprox_fit(0.4:0.01:2)/3600,'k--')
scatter(chosen(:,4),1/3600*kelong_nutrient, 50, color_list, 'filled')
% analytical growth rate vs synthesis speed 
% scatter(gr_ana_nutrient, ksyn_ana_nutrient/3600)  
% legend('Klump et al. 2013','Quadratic approx','Theory fit','Numerical results','Location','SouthEast')
legend('Klump et al. 2013','Eq. S48','Eq. S66','Numerical results','Location','SouthEast')
xlim([0 2])
xticks(0:0.4:2)
ylim([0 20])
yticks(0:5:20)
xlabel('\lambda^* (1/h)')
ylabel('k^{*}_{syn} (1/s)')
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)
box on


%Fig. S8a
%plot each kmet individually to visualize at a given accuracy 

%for SI figure to illustrate max 
m = 3; %index of nutrient quality --> kmet = 0.18 1/s
Mpicklist(m)
blist(m)
clist(m)
kmetsample(m)/3600
index_check = 3;  %index of the row of accuracy of interest
Alist_p1(index_check)
ksynmaxlist_p1(index_check)
KMefflist_p1(index_check)

%c_factor to convert aa to Ternary complex
kelongf_all = @(kmax,KM,y,Alist,b,c) T3_A(Alist,b,c).*3600.*kmax.*c_factor.*y(:,3)./(c_factor*y(:,3)+KM.*1e-6.*NA.*V(y)*1e-15);
%metabolic rate
kmetf_all = @(y,kmet) kmet./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
%input flux
JM_f = @(y,kmet) y(:,2).*kmetf_all(y,kmet);
%output flux
JS_f = @(y,kmax,KM,Alist,b,c) y(:,1).*kelongf_all(kmax,KM,y,Alist,b,c);
%degradation flux
Jd_f = @(y,Alist) degrade(Alist).*(NR*y(:,1)+NP*y(:,2));

JM_interest = JM_f(odestore_kmet{m}{index_check}, kmetsample(m));
JS_interest = JS_f(odestore_kmet{m}{index_check}, ksynmaxlist_p1(index_check), KMefflist_p1(index_check), Alist_p1(index_check), blist(index_check), clist(index_check));
Jd_f = Jd_f(odestore_kmet{m}{index_check}, Alist_p1(index_check));
%check flux
%P*kmet = lambda*M - degrade flux 
JM_interest./(odestore_kmet{m}{index_check}(:,4)*Mpicklist(m) + Jd_f);
%P*kmet = R*ksyn + lambda*aa
JM_interest./(odestore_kmet{m}{index_check}(:,4).*odestore_kmet{m}{index_check}(:,3) + JS_interest);

JM_max = kmetsample(m)*odestore_kmet{m}{index_check}(:,2);
JS_max = ksynmaxlist_p1(index_check)*3600*odestore_kmet{m}{index_check}(:,1);

%find max growth rate
odestore = odestore_kmet{m};
[val, id] = max(odestore{index_check}(:,4)); 

flux_pick = NaN*ones(3,6);
%1st - 2 - 3 - 4 - 5 col: R mass fraction - uR - JM - JS - JMmax - JSmax
rgap = [-5 0 5];
for i = 1 : 3
    if i == 3 
        if id - rgap(i) >= 1
            flux_pick(i,:) =  [odestore{index_check}(id-rgap(i),7) odestore{index_check}(id-rgap(i),5) JM_interest(id-rgap(i)) JS_interest(id-rgap(i)) JM_max(id-rgap(i)) JS_max(id-rgap(i))]; 
        else       
            flux_pick(i,:) =  [odestore{index_check}(id-rgap(i)+1,7) odestore{index_check}(id-rgap(i)+1,5) JM_interest(id-rgap(i)+1) JS_interest(id-rgap(i)+1) JM_max(id-rgap(i)+1) JS_max(id-rgap(i)+1,1)]; 
        end
    else
        flux_pick(i,:) =  [odestore{index_check}(id-rgap(i),7) odestore{index_check}(id-rgap(i),5) JM_interest(id-rgap(i)) JS_interest(id-rgap(i)) JM_max(id-rgap(i)) JS_max(id-rgap(i))]; 
    end
end

% frange = find(abs(JM_interest-JM_interest(id))<=0.3*JM_interest(id)); %find within certain range
frange = id - rgap(3)-3 : id - rgap(1) + 2;
%Fig. S8a
figure
b = bar(flux_pick(:,1),[flux_pick(:,2).*flux_pick(:,4) (1-flux_pick(:,2)).*flux_pick(:,4)],0.1,'stacked')
b(1).FaceColor = [105 105 105]/225;
b(2).FaceColor = [1 1 1];
hold on
scatter(flux_pick(:,1),flux_pick(:,3),'filled','MarkerFaceColor',[0, 0.4470, 0.7410])
plot(flux_pick(:,1),flux_pick(:,5:6),'sq--')
legend('R','P','J_M','J_M^{max}','J_S^{max}','Location','northwest')
xlabel('\phi_R')
ylabel('Flux (\times 10^9 aa/h)')

%Fig. S8b
%plot growth rate vs ribosomal mass fraction
figure
store = odestore{index_check};
plot(store(frange,7), store(frange,4))
xlabel('\phi_R')
ylabel('\lambda (1/h)')
xlim([0.26 0.36])
xticks([0.26:0.02:0.36])
yticks([1.26:0.02:1.34])

%Fig. S9 is at the top section of this script


%Fig. S10a and Fig. S10b
namelist = {'R','P'};
masslist = [NR NP];

for i = 1 : length(namelist)
    if i == 1
        yL = sprintf('%s^* (x 10^4 molecules/cell)',namelist{i}); 
    else 
        yL = sprintf('%s^* (x 10^6 molecules/cell)',namelist{i}); 
    end
    yR = ['\phi_' sprintf('{%s}^*',namelist{i})];       
    figure
    yyaxis left
    plot( chosen(:,4),   chosen(:,i),'.-','MarkerSize',20)
    if i == 1
        ylim([0 12e4])
    else
        ylim([2 7]*1e6)
    end
    ylabel(yL)
    xlabel('Growth rate (1/h)')
    yyaxis right
    plot( chosen(:,4),  chosen(:,i+6),'sq-')
    ylabel(yR)        
    set(gca,'FontName','Helvetica','FontSize',16)
end


%intake flux partition
%Fig. S10c
figure
b = bar(chosen(:,4), [chosen(:,4).*chosen(:,11) flux_nutrient(:,3) ], 'stacked')
% b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.9290    0.6940    0.1250];
hold on
scatter(chosen(:,4), flux_nutrient(:,1),[],'m','filled')
ylabel('Flux (\times 10^9 aa/h)')
xlabel('Growth rate (1/h)')
legend('Growth flux','Degradation flux','Input flux')
set(gca,'FontName','Helvetica','FontSize',16)
axes('Position',[.2 .5 .3 .3])   %[xstart ystart xend-xstart yend-ystart ]
box on
bar(chosen(:,4),flux_nutrient(:,3),'FaceColor',[0.9290    0.6940    0.1250])
xlim([0 2.5])
xticks(0:0.5:2.5)
yticks([0:2.5:10]*1e7)
ylabel('Flux (\times 10^7 aa/h)')
xlabel('Growth rate (1/h)')
set(gca,'FontName','Helvetica','FontSize',12)

%Fig. S11 is ran 
run('slow_growth_usethisone.m')


%plot a common factor for R = a*gr^2 and JM = 1/a*M^2
%a_ave obtained from Fig. 2d
%Fig. S12a
figure
plot(R_approx2(0:0.1:2), 0:0.1:2, 'k')
hold on
scatter(chosen(:,1), chosen(:,4),50,color_list,'filled')
scatter(dataextract(:,1), dataextract(:,3),'sq','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1)
scatter(Rpercell, gr, 50,'x','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1)

ylabel('\lambda^* (1/h)')
xlabel('R^* (\times10^4 molecules/cell)')
yticks(0:0.4:2) 
ylim([0 2])
xlim([0 13e4])
xticks([0:2e4:12e4]) 
legend('Theoretical line','Numerical results','Dai et al. 2016','Bremer & Dennis 2008','Location','NorthWest')
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)

%Fig. S12b
figure
plot(1e9*[0:0.01:4], JM_approx3(1e9*[0:0.01:4]),'k')
hold on
scatter(Mpercell*0.6,Rpercell.*kselect,50,'kx','LineWidth',1.25) %Bremer - Dennis 
scatter(dataextract(:,5)*0.6*1e-12/maa,dataextract(:,1).*dataextract(:,4)*3600,50,'k^','LineWidth',1.25) %Hwa-Liu
scatter(chosen(:,11), JM, 50, color_list, 'filled')
xlabel('M^* (\times 10^9 aa)')
ylabel('J_M^* (\times 10^9 aa/h)')
xlim([0 6e9])
ylim([0 8e9])
yticks([0:2:8]*1e9)
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)
legend('Theory line','Bremer & Dennis 2008','Zheng et al. 2020, Dai et al. 2016','Numerical results', 'Location', 'SouthEast')

