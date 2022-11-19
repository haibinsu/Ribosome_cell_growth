%for growth rate < 0.4 1/h 
%from Zheng et al  https://doi.org/10.1038/s41564-020-0717-x

m_small = [0.34, 0.45, 0.46, 0.49, 0.53, 0.57, 0.54]; % OD600ml/1e9 cell
gr_small = [0.06, 0.18, 0.19, 0.31, 0.34, 0.40, 0.40];

gr_m_fit = fit(m_small',gr_small','poly2');

%gr as function of mass 
f_gr_m = @(x) gr_m_fit.p1*x.^2 + gr_m_fit.p2*x + gr_m_fit.p3;
f_scaling = @(x) 2.85*x.^3.34;

masslist_small = 0.3:0.01:0.6;

%Fig. S11a
figure
scatter(m_small, gr_small, 'filled')
hold on
% plot(masslist_small,f_gr_m(masslist_small))
plot(masslist_small, f_scaling(masslist_small))
xlabel('Average cell mass (OD600 ml per 10^9 cells)')
ylabel('Growth rate (1/h)')
ylim([0.01 0.5])
yticks([0.05:0.05:0.5])
yticklabels({'','0.1','','0.2','','0.3','','0.4','','0.5'})
box on
set(gca,'FontName','Helvetica','FontSize',16)

%% assume a relation between fastest growth rate and initial mass
%use KMeff from part 1 for ternary complex instead of rescaling 

%starting point is growth_HB.m 

% remember to run bionumers.m
run('bionumers.m')
close all
% NR = 7336; %aa/ribosome
% NP = 300; %aa/metabolic protein 
% rho = 7.4e8; %aa/um^3 
% NA = 6.02e23;
% ksynmax = 22*3600; %aa/h
% lambdamax = 17.83*3600; %aa/h (synthesis speed in term of growth rate)
% % Ka = 2e-3; %M
% % Ke = 2e-1; %M
% RPpercent = 0.7;
% %dry weight = 1/3 wet weight
% dwconv = 0.33; 
% %new value of Ka Ke to try
% maa = 110*1.66e-24; % g = 110 Da %average mass of amino acid
%www.pnas.org/cgi/doi/10.1073/pnas.1421138111 

% ksynmax = 20*3600; %aa/h
Ka = 0.25e-3; %M
% Ke = 1e-2; %M %5e-2 1e-2
Ke = 2.5e-3;
Jmtm = 5.2e4/6 ; %convert ATP to aa because 6 ATP required for 1 peptide bond
%https://doi.org/10.1016/j.cels.2019.06.003

%use KMeff from part 1 for ternary complex instead of rescaling 
k2 = 100*3600; %1/h max charging speed 
% k2 = 80*3600;
R = 20; %uM 
% S = R/1.5; % --> charging synthase/ribosome ratio = 1/1.5 
S = R/5;
ksynmaxT3 = 22*3600; %1/h
% KM = 20e-6; %M
KM = 10e-6; %M 
KMeff = KM;

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

%mapping aa to T3
%exact mapping 
% T3convert_T3 = @(y)
% k2*S/R*KM*NA*V(y).*1e-15.*y(:,3)./(ksynmaxT3*Ka*NA*V(y)*1e-15+(ksynmaxT3-k2*S/R)*y(:,3));
%linear mapping
T3convert_T3 = @(y) k2*S/R*KM*y(:,3)./(ksynmaxT3*Ka);

%translation rate using T3
ktrans_T3 = @(y) ksynmaxT3*T3convert_T3(y)./(T3convert_T3(y)+KMeff*NA*V(y)*1e-15);

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

%mapping aa to T3
% T3convertf_T3 = @(y) k2*S/R*KM*NA*Vf(y)*1e-15*y(3)./(ksynmax*Ka*NA*Vf(y)*1e-15+(ksynmax-k2*S/R)*y(3));
T3convertf_T3 = @(y) k2*S/R*KM*y(3)./(ksynmax*Ka);
c_factor = (k2*S/R*KM)/(ksynmax*Ka)
%translation rate based on T3
ktransf_T3 = @(y) ksynmax*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

flist = {Mf,Vf,ktransf_T3,T3convertf_T3};

% kmetsample = [0.013 0.026 0.032 0.044 0.05]'*3600;

% kmetsample = [0.003 0.008 0.013 0.026 0.032]'*3600;

% kmetsample = [0.006 0.01 0.015 0.025 0.032]'*3600;

%function for gr - mass relation
%use to get Mpick
%x is growth rate
b = 1.09/4*1e-30;
% f_gr_mass = @(x) (x/b).^(1/3.34);
% f_gr_mass(odemax(:,4))

kmetsample = [0.016 0.02 0.025 0.029 0.045]'*3600;
Mpicklist = 1e9*[0.43, 0.51,  0.58,  0.69,    0.85];


nstore = NaN*ones(length(kmetsample),6);
%1 - 6th col: #R - #P - #aa - gr1 - gr2 - M0

nstore(:,6) = Mpicklist';

nSim = 70;
uRlist = 0.01:0.01:0.9;

odestore = cell(length(kmetsample),1);
odemax = NaN*ones(length(kmetsample),14);
%1st - 6th col: R - P - aa - gr - uR - uP 
%7th - 10th col: phiR - phiP - phiaa - ternary complex T3 
%11 - 14 col: JM - JS - JQ - total cell mass M0 
fluxtore = NaN*ones(length(kmetsample),3);


% uS = @(y) 0.38*Mf(y).^-3.34;

%function of cell mass
% uS = @(x) (0.38*x-4.17e7)./(b*x.^4.34);
uS = @(y) (0.38*Mf(y)-4.17e7)./(b*Mf(y).^4.34);
% maintenance flux = 0.38*x-4.17e7 
% synthesis flux = b*Mf(y).^4.34
% protein synthesis flux * uS_mass = maintenance flux \propto JS/lambda =

flist{6} = uS;
f_ratio = @(x) (0.38*x-4.17e7)./(b*x.^4.34);
 
% plot(mlist(non_idx_m),  f_ratio(mlist(non_idx_m)))

for j = 1 : length(kmetsample)
    kmet = kmetsample(j);
    %metabolic rate
    kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
    flist{5} = kmetf; 
%     kd = @(y) 0.1*Mf(y)/0.5e9;
%     flist{6} = kd;
    store = NaN*ones(length(uRlist),13);
    %1st - 10th col: R - P - aa - gr - uR - uP - phiR - phiP - phiaa - T3 
    %11 - 14: JM - JS - JQ - M 
    
    Mth = nstore(j,6); % note that 2*Mth = mass threshold for cell division
%     Mth = 2.2048e+09; %fixed mass threshold
    parfor i =  1 : length(uRlist)

        R0 = 2e4;
        P0 = 1e5;
        M0 = 3e8; %aa
        aa0 = M0 - (NR*R0 + P0*NP);
%         R0 = nstore(j,1);
%         P0 = nstore(j,2);
%         aa0 = nstore(j,3);
        y0 = [R0 P0 aa0];
        
%         store_uS = NaN*ones(length(uSlist),10);
        
%         for w = 1 : length(uSlist)
            
%             par = [ksynmax kmet uRlist(i) NR NP uSlist(w)];
        par = [ksynmax kmet uRlist(i) NR NP];

        k = 1;
        t_tot = 0;
        y_tot = y0;
        tspan = [0 2000];  
        gr = 1;
        Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));
        %run to reach exponential steady state 
        while k <= nSim
            [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_slow_partition(t,y,par,flist),tspan,y0, Opt);
            semilogy(t,y)
            semilogy(t,M(y))
            if isempty(ye) ~= 1
                t_tot = [t_tot ; [t;te] + t_tot(end)];
                y_tot = [ y_tot ; [y;ye] ];    
                y0 = ye/2;
                k = k + 1;
%                 dy = ode_growth_findmax(t(end),y(end,:),par);
%                 gr = sum(dy.*[NR;NP;1])/sum(y(end,:).*[NR NP 1]);
                %growth rate at current division sets future threshold mass 
%                 Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));            
            else
                gr = NaN;
                break;
            end
        end
    %     semilogy(t_tot,y_tot)
        if isnan(gr) ~= 1 
        [t,y,te,ye,ie] = ode15s(@(t,y)   ode_growth_slow_partition(t,y,par,flist),[0:1e-3:2*te],y0,Opt);
    %     semilogy(t,y)
%         semilogy(t,M(y))
        %check flux balance
%         ye(2)*kmetf(ye)./((1+uS(ye))*ye(1)*ktrans_T3(ye)+log(2)/t(end)*ye(3))
        gfit = fit(t,log(M(y)),'poly1');
        plot(gfit, t, log(M(y)))
        
               
        kmet_time = @(y) kmet./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
        uS_time = @(y)(0.38*M(y)-4.17e7)./(b*M(y).^4.34);
        
        JM_time = y(:,2).*kmet_time(y);
        JS_time = y(:,1).*ktrans_T3(y);
        JQ_time = uS_time(y).*JS_time;
        
%         plot(t,JM_time./(JS_time + JQ_time))
%         plot(t, JM_time)
%         hold on
%         plot(t,JS_time + JQ_time)
        
        JM_ave = trapz(t,JM_time)*log(2)/(t(end)-t(1));
        JS_ave = trapz(t,JS_time)*log(2)/(t(end)-t(1));
        JQ_ave = trapz(t,JQ_time)*log(2)/(t(end)-t(1));
        JM_ave/(JS_ave + JQ_ave);
        %evaluate dy/dt 
%             dy_dt = NaN*ones(length(t),3);
%             for p = 1:length(t)
%                 dy_dt(p,:) = ode_growth_slow(t(p),y(p,:),par,flist);
%             end
%             hfit = fit(M(y),M(dy_dt),'poly1');
%         plot(hfit, M(y),M(dy_dt))
%         store(i,:) = [y0 log(2)/t(end) uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];
%         store_uS(w,:) = [y0 gfit.p1 uSlist(w) 1 - uSlist(w) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];
%             m_index = 8;
%             store_uS(m_index,1)*ktrans_T3(store_uS(m_index,:))/(store_uS(m_index,4)*Mf(store_uS(m_index,:)))
        store(i,:) = [y0 gfit.p1 uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0) JM_ave JS_ave JQ_ave];
%         m_index = 1;
%         store(m_index,1)*ktrans_T3(store(m_index,:))/(store(m_index,4)*Mf(store(m_index,:)))

        end
    
end
    odestore{j} = store;
    [val, idx] = max(store(:,4));
%     JMmax = store(idx,2)*kmetf(store(idx,:));
%     JSmax = store(idx,1)*ktrans_T3(store(idx,:));
%     J_maint = 0.1*Mf(store(idx,:)) - 2.98e7;
%     J_maint = kd(store(idx,:))*(NR*store(idx,1) + NP*store(idx,2) + store(idx,3));
%     J_maint = uS(store(idx,:))*JSmax;

    JMmax = store(idx,11);
    JSmax = store(idx,12); 
    J_maint = store(idx,13);
    odemax(j,:) = [store(idx,:) M(store(idx,:))];
    fluxtore(j,:) = [JMmax JSmax J_maint];
end

masslist = 10.^(8:0.01:8.95);

flux_fit = fit(M(odemax), fluxtore(:,1), 'poly2');
JM_f = @(x) flux_fit.p1*x.^2 + flux_fit.p2*x  + flux_fit.p3;
JM_f2 = @(x) 3.25e-10*x.^2;


m_gr_f = @(x) b*x.^3.34;

%% fast growth
% remember to run bionumers.m
run('bionumers.m')
close all

Ka = 0.25e-3; %M
% Ke = 1e-2; %M %5e-2 1e-2
Ke = 2.5e-3;
Jmtm = 5.2e4/6 ; %convert ATP to aa because 6 ATP required for 1 peptide bond
%https://doi.org/10.1016/j.cels.2019.06.003

%use KMeff from part 1 for ternary complex instead of rescaling 
k2 = 100*3600; %1/h max charging speed 
% k2 = 80*3600;
R = 20; %uM 
% S = R/1.5; % --> charging synthase/ribosome ratio = 1/1.5 
S = R/5;
ksynmaxT3 = 22*3600; %1/h
% KM = 20e-6; %M
KM = 10e-6; %M 
KMeff = KM;

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

%mapping aa to T3
%exact mapping 
% T3convert_T3 = @(y)
% k2*S/R*KM*NA*V(y).*1e-15.*y(:,3)./(ksynmaxT3*Ka*NA*V(y)*1e-15+(ksynmaxT3-k2*S/R)*y(:,3));
%linear mapping
T3convert_T3 = @(y) k2*S/R*KM*y(:,3)./(ksynmaxT3*Ka);

%translation rate using T3
ktrans_T3 = @(y) ksynmaxT3*T3convert_T3(y)./(T3convert_T3(y)+KMeff*NA*V(y)*1e-15);

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

%mapping aa to T3
% T3convertf_T3 = @(y) k2*S/R*KM*NA*Vf(y)*1e-15*y(3)./(ksynmax*Ka*NA*Vf(y)*1e-15+(ksynmax-k2*S/R)*y(3));
T3convertf_T3 = @(y) k2*S/R*KM*y(3)./(ksynmax*Ka);
c_factor = (k2*S/R*KM)/(ksynmax*Ka)
%translation rate based on T3
ktransf_T3 = @(y) ksynmax*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

flist = {Mf,Vf,ktransf_T3,T3convertf_T3};


kmetsample_fast = [0.08 0.14 0.16 0.23 0.3]'*3600;
Mpicklist_fast = 1e9*[1.13,  1.36,    1.82,    2.26,    2.74];

nstore_fast = NaN*ones(length(kmetsample_fast),6);
%1 - 6th col: #R - #P - #aa - gr1 - gr2 - M0

nstore_fast(:,6) = Mpicklist_fast';

nSim = 50;
uRlist = 0.01:0.01:0.9;
odestore_fast = cell(length(kmetsample_fast),1);
odemax_fast = NaN*ones(length(kmetsample_fast),14);
%1st - 6th col: R - P - aa - gr - uR - uP 
%7th - 11th col: phiR - phiP - phiaa - ternary complex T3 - total cell mass M0 
fluxtore_fast = NaN*ones(length(kmetsample_fast),3);

% uS_fast = @(y) (0.38*Mf(y)-4.17e7)./(1/(1.6e9)*Mf(y).^2);
% uS_mass = @(x) (0.38*x-4.17e7)./(1/(1.6e9)*x.^2);
uS_fast = @(y) (0.38*Mf(y)-4.17e7)./(1/(1.65e9)*Mf(y).^2);
uS_mass = @(x) (0.38*x-4.17e7)./(1/(1.65e9)*x.^2);
% maintenance flux = 0.38*x-4.17e7 
% synthesis flux = 1/(1.65e9)*x.^2
% protein synthesis flux * uS_mass = maintenance flux \propto JS/lambda =

flist{6} = uS_fast; 
for j = 1 : length(kmetsample_fast)
    kmet_fast = kmetsample_fast(j);
    %metabolic rate
    kmetf_fast = @(y) kmet_fast./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
    flist{5} = kmetf_fast; 
%     kd = @(y) 0.1*Mf(y)/0.5e9;
%     flist{6} = kd;
    store_fast = NaN*ones(length(uRlist),13);
    %1st - 10th col: R - P - aa - gr - uR - uP - phiR - phiP - phiaa - T3 
    %11 - 14: JM - JS - JQ - M 
    
    Mth_fast = nstore_fast(j,6); % note that 2*Mth = mass threshold for cell division
%     Mth = 2.2048e+09; %fixed mass threshold
    parfor i =  1 : length(uRlist)

        R0 = 2e4;
        P0 = 1e5;
        M0 = 3e8; %aa
        aa0 = M0 - (NR*R0 + P0*NP);
        y0 = [R0 P0 aa0];
        

        par = [ksynmax kmet_fast uRlist(i) NR NP];

        k = 1;
        t_tot = 0;
        y_tot = y0;
        tspan = [0 2000];  
        gr = 1;
        Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth_fast,par));
        %run to reach exponential steady state 
        while k <= nSim
            [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_slow_partition(t,y,par,flist),tspan,y0, Opt);
            semilogy(t,y)
            semilogy(t,M(y))
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
        [t,y,te,ye,ie] = ode15s(@(t,y)   ode_growth_slow_partition(t,y,par,flist),[0:1e-3:2*te],y0,Opt);
 
        gfit = fit(t,log(M(y)),'poly1');
        plot(gfit, t, log(M(y)))
        
               
        kmet_time_fast = @(y) kmet_fast./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
        uS_time_fast = @(y)(0.38*M(y)-4.17e7)./(1/(1.65e9)*M(y).^2);
        
        JM_time_fast = y(:,2).*kmet_time_fast(y);
        JS_time_fast = y(:,1).*ktrans_T3(y);
        JQ_time_fast = uS_time_fast(y).*JS_time_fast;
        

        
        JM_ave_fast = trapz(t,JM_time_fast)*log(2)/(t(end)-t(1));
        JS_ave_fast = trapz(t,JS_time_fast)*log(2)/(t(end)-t(1));
        JQ_ave_fast = trapz(t,JQ_time_fast)*log(2)/(t(end)-t(1));
        JM_ave_fast/(JS_ave_fast + JQ_ave_fast);
        store_fast(i,:) = [y0 gfit.p1 uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0) JM_ave_fast JS_ave_fast JQ_ave_fast];

        end
    
end
    odestore_fast{j} = store_fast;
    [val, idx] = max(store_fast(:,4));


    JMmax_fast = store_fast(idx,11);
    JSmax_fast = store_fast(idx,12); 
    J_maint_fast = store_fast(idx,13);
    odemax_fast(j,:) = [store_fast(idx,:) M(store_fast(idx,:))];
    fluxtore_fast(j,:) = [JMmax_fast JSmax_fast J_maint_fast];
end

odemax_fast(:,14)./Mpicklist_fast'
M(odemax_fast)./odemax_fast(:,14)

% scatter(odemax_fast(:,4), odemax_fast(:,14))
mass_fit = fit(odemax_fast(:,4), odemax_fast(:,14),'poly1')
% plot(mass_fit, odemax_fast(:,4), odemax_fast(:,14))


masslist_fast = 10.^(8.95:0.01:9.8);

flux_fit_fast = fit(M(odemax_fast), fluxtore_fast(:,1), 'poly2')
JM_f_fast = @(x) flux_fit_fast.p1*x.^2 + flux_fit_fast.p2*x  + flux_fit_fast.p3;

JM_f2_fast = @(x) 5.94e-10*x.^2;


m_gr_f_fast = @(x) (x-mass_fit.p2)/mass_fit.p1;


%Fig. S11a
figure
scatter(m_small, gr_small, 'filled')
hold on
% plot(masslist_small,f_gr_m(masslist_small))
plot(masslist_small, f_scaling(masslist_small))
xlabel('Average cell mass (OD600 ml per 10^9 cells)')
ylabel('Growth rate (1/h)')
ylim([0.01 0.5])
yticks([0.05:0.05:0.5])
yticklabels({'','0.1','','0.2','','0.3','','0.4','','0.5'})
box on
set(gca,'FontName','Helvetica','FontSize',16)

%Fig. S11b
m0 = 1e10; %aa
figure
scatter(M(odemax), m0*odemax(:,4)./fluxtore(:,1), [], [0, 0.447, 0.741], 'filled')
hold on
plot(masslist, m0*m_gr_f(masslist)./JM_f2(masslist))
scatter(M(odemax_fast), m0*odemax_fast(:,4)./fluxtore_fast(:,1), [], [0.8500, 0.3250, 0.0980], 'filled')
% plot(masslist_fast, m_gr_f_fast(masslist_fast)./JM_f2_fast(masslist_fast))
plot(masslist_fast, m0*m_gr_f_fast(masslist_fast)./JM_f_fast(masslist_fast))
xlim([4e8 4e9])
ylim([2 9])
set(gca,'XScale', 'log')
% xticks(10.^[8.6021:0.1:9.6021])
xticks([0.4e9 0.5e9 0.6e9 0.7e9 0.8e9 0.9e9 1e9 2e9 3e9 4e9])
xticklabels({'0.4','','0.6','','0.8','','1','2','3','4'})
set(gca,'FontName','Helvetica','FontSize',16)
xlabel('M^\ast (\times 10^9 aa)')
ylabel('Z')
% ylabel('\lambda^\ast/J_{M}^\ast (\times 10^{-10} 1/aa)')
% ylabel('FoM (m_0\lambda^\ast/J_{M}^\ast)')
box on

