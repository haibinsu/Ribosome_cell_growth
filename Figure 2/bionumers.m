NA = 6.02e23; %Avogardo number
%mass of ribosome 
MR = 2.7e6*1.66e-24; %unit is g because 2.7 MDa , 1 Da = 1.66e-24 g
%https://pubmed.ncbi.nlm.nih.gov/16753153/
%https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=10&id=100118
% NR = 7336; % #aa/ribosome
%extended mass of ribosome table S1 of 10.1111/febs.13258
NR = 12307;  % #aa/ribosome 
NP = 300; % #aa/protein http://book.bionumbers.org/how-big-is-the-average-protein/
ksynmax = 22*3600; %synthesis speed #aa/(h.ribosome)
maa = 110*1.66e-24; % g = 110 Da %average mass of amino acid
%www.pnas.org/cgi/doi/10.1073/pnas.1421138111 
rho = 7.4e8; %aa/um^3  calculated from SI of Inflating bacterial cell
Ap = 14e-6; %um^2 surface area of transpoter SI of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2799844/
%http://ccdb.wishartlab.com/cgi-bin/STAT_NEW.cgi
vR = 3.4e6*1e-12; %um^3 volume of ribosome
Pr = 2.5; %nm %average radius of protein 
vP = 4/3*3.14*(Pr*1e-3)^3; %um^3 volume of protein 
vaa = 120e-12; %um^3 volume of amino acid http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html#refs

mtRNA = 25000; %Da https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=5&id=101177
mEFTu = 43000; %Da 
mGTP = 519; %Da
%mass of ternary complex 
mT3 = (mtRNA + mEFTu+mGTP+110)/110; %aa where maa = 110 Da
NT3 = 624; %aa 
%RNA + protein account for 70% of cell dry weight
RPpercent = 0.7;
% http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/#:~:text=Protein%20is%20evaluated%20at%20%E2%89%88,metabolites%2C%20ions%20etc.).
%https://link.springer.com/article/10.1007/BF02578899#:~:text=The%20amount%20of%20ATP%20required%20for%20the%20formation%20of%20microbial,under%20various%20conditions%20was%20calculated.&text=The%20results%20of%20the%20calculations,CoA%20is%20formed%20from%20glucose.

%dry weight = 1/3 wet weight
dwconv = 0.33; 
% http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/#:~:text=Protein%20is%20evaluated%20at%20%E2%89%88,metabolites%2C%20ions%20etc.).
% 10.1128/AEM.48.4.755-757.1984

NT = 400; %aa https://www.sciencedirect.com/science/article/pii/B0124437109005299
KaT = 3.9; % uM M-M constant for ternary complexes

nATP = 1e19; % 1 J = 1e19 ATP molecules
%http://book.bionumbers.org/how-much-energy-is-released-in-atp-hydrolysis/
%https://www.ncbi.nlm.nih.gov/books/NBK7919/

nconv = 6; %6 ATP needed for aa incorporation into peptide chain
%www.pnas.org/cgi/doi/10.1073/pnas.1421138111 

%10.15252/msb.20209478 
%dry mass
mdry = 509 ; %ug/(mL*OD600)
%number of cell
ncell = 1e9; %CFU/(mL*OD600)
%--> dry weight per cell
dwpercell = mdry/ncell ; %ug/CFU = ugDW/cell
%maintenance flux under starvation 
Jmtn_starve = 0.5*1e-15*6.02e23*15; % glycerol/day/cell which is converted from 0.5 fmol glycerol/(day*CFU)
%1 glycerol = 15 ATP 

%https://doi.org/10.1016/S0021-9258(18)94507-2
% wet weight
mwet = 0.89*1e-9; %mg/cell because 0.89 mg = 10^9 cell
% dry weight
dwpercell_2 = mwet*dwconv; %mg/cell

(dwpercell_2*1e-3)/(dwpercell*1e-6) %check ratio
%average dry weight per cell for E coli = 2.9e-13 g 
%SI 2 or 3 (tab soluable pool) of https://www.embopress.org/doi/full/10.1038/msb4100155
%meanwhile the dry weight details depending on growth rate is here
%http://book.bionumbers.org/how-big-is-an-e-coli-cell-and-what-is-its-mass/
%it varies from 150 fg (dbl time = 100 min) to 870 fg (dbl time = 24 min)

%maintenance flux under exp growth with aerobic glucose
%vary from 7.6 - 8.4 mmol ATP/gDW/h
%two references from 10.15252/msb.20209478 
%one of them is https://www.embopress.org/doi/full/10.1038/msb4100155
Jmtn_expgr = 8*1e-3*6.02e23/(1/(dwpercell*1e-6)) ; %ATP/h/cell
%vary within 8*1e-3*6.02e23/(1/(150e-15)) - 8*1e-3*6.02e23/(1/(870e-15))
% around 7e8 - 4.5e9 ATP/h/cell

%note that from https://www.embopress.org/doi/full/10.1038/msb4100155
%maintenance for non growth ~ 8.4 mmol ATP/gDW/h
%energy needed for growth ~ 60 mmol ATP/gDW/h = 22.87 for protein synthesis + 36.9 for other
% --> synthesis accounts for ~ 40% energy for growth
% --> maintenance = 0.14 energy for growth = 0.12 total energy
%note that energy required to maintain electrochemical gradient ~ 51% total energy
%https://openwetware.org/wiki/BioNumber_Of_The_Month

%1st growth law data - table S1 in SI of 10.1126/science.1192588
%growth rate (1/h) - Hwa
glawdataH(:,1) = [0.4;0.57;0.71;1;1.31;1.58];
%R-protein fraction mass fraction - Hwa 
glawdataH(:,2) = 0.76*[0.177;0.230;0.224;0.287;0.414;0.466];
%growth rate (1/h) - Forchhammer & Lindahl
glawdataF(:,1) = [0.38;0.6;1.04;1.46;1.73];
%R-protein fraction mass fraction - Forchhammer & Lindahl 
glawdataF(:,2) = 0.76*[0.189;0.224;0.295;0.421;0.469];
%growth rate (1/h) - Forchhammer & Lindahl
glawdataB(:,1) = [0.42;0.69;1.04;1.39;1.73];
%R-protein fraction mass fraction - Bremer & Dennis 2009
glawdataB(:,2) = 0.76*[0.2;0.255;0.331;0.391;0.471];

glawdata = vertcat(glawdataH,glawdataF,glawdataB);

scatter(glawdata(:,1),glawdata(:,2))        

load('data_Dai')  %10.1038/nmicrobiol.2016.231
%1st col: growth rate (1/h)
%2nd col: peptide synthesis speed (1/s)

load('Liudata')
%from https://doi.org/10.1038/s41564-020-0717-x
%1st col: growth rate (1/h)
%2nd col: (R+P)/cell (pg)
%3rd col: dry mass/cell (pg)

rhoLiu = 0.29; %pg/um^3
%find slope of dry mass vs growth rate 
id = find(isnan(Liudata(:,3)) ~= 1);
hfit = fit(Liudata(id,1),Liudata(id,2)*1e-12/maa,'poly1')

figure
scatter(Liudata(:,1),Liudata(:,2)*1e-12/maa,'filled')
xlabel('Growth rate (1/h)')
ylabel('R+P per cell (aa)')
xlim([0 2])
ylim([5e8 6e9])

figure
scatter(Liudata(:,1),Liudata(:,2),'filled')
xlabel('Growth rate (1/h)')
ylabel('R+P per cell (pg)')
xlim([0 2])

id = find(isnan(Liudata(:,3)) ~= 1);
gfit = fit(Liudata(id,3),Liudata(id,2),'poly1');
figure
scatter(Liudata(:,3),Liudata(:,2),'filled')
hold on
plot(0:0.1:1.5,0.73*(0:0.1:1.5)+0.02)
xlabel('Dry mass per cell (pg)')
ylabel('R+P per cell (pg)')

%cell dry mass in Bremmer and Dennis  10.1128/ecosal.5.2.3
Mpercell = [226 374 555 774 921 1023]'/1e9; %ugDW 
Mpercell = Mpercell*1e-6/maa;
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h

figure
scatter(gr,Mpercell*0.7,'sq','filled') %because R+P = 0.7 cell dry mass
hold on
scatter(Liudata(:,1),Liudata(:,2)*1e-12/maa,'filled')
plot(0:0.1:2,2e8 + 1.65e9*(0:0.1:2))  %the trend used in growth_HB.m


