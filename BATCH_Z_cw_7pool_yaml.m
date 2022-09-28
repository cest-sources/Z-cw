
%% analytic Z-spectra
%   Date: 2014/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  - moritz.zaiss@fau.de
%   CEST sources  Copyright (C) 2020  Moritz Zaiss
%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   **********************************
%
%   --SHORT  DOC--
% Figure 1 from the bookchapter Advanced Neuro MRI
% Amide pool in a water and MT environment reflecting white brain matter

%% This file requires the pulseq CEST librares in addition
 if exist('pulseq-cest', 'dir')
        disp('pulseq-cest-library already installed, skip...')
 else
    system('git clone -b major_update https://github.com/kherz/pulseq-cest'); 
    cd pulseq-cest;
    install_pulseqcest;
    cd ..
    pulseqCEST_simlib=[pwd '/pulseq-cest-library/sim-library/'];
 end

%% SETUP
clearvars P Pref Pstart
clc

% bmsimfile='WM_3T_Stanisz2005_5pool_bmsim.yaml';
bmsimfile='WM_3T_Wang2020_5pool_bmsim.yaml';

Psim = readSimulationParameters([pulseqCEST_simlib bmsimfile]);

% setup sequence parameters
    P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
    P.FREQ=300;             % static B0 field [MHz] ~7T ; ppm and µT are used for offsets and B1, therefore gamma=267.5153 is given in Hz.
    P.B1=0.75;                 % irradiation amplitude [µT]
    P.tp=2;                % pulse duration = saturation time [s]
    P.xZspec= [-6:0.05:6];   % chemical shift of the CEST pool in [ppm] 

    Pstart=P;
% PLOT Z_cw(P)
[Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw_yaml(P,Psim);
Pref=Psim;    
Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0; 
[Zref ,Rex_ref]=Z_cw_yaml(P,Pref);
    
figure(1), 
set(gcf,'Position',[587   262   560   720]);
t = annotation('textbox','String','a)','Position',[0 0.85 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','b)','Position',[0 0.63 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','c)','Position',[0 0.42 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','d)','Position',[0 0.2 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','e)','Position',[0.55 0.2 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
set(groot,'defaultLineLineWidth',1.2)
subplot(4,1,1), 
yyaxis left
% plot(P.xZspec,P.xZspec*0+R1obs,':', 'Displayname','R_{1,obs}') ;   hold on;
plot(P.xZspec,Rpw+Rpmt,'--',P.xZspec,Rpw+Rpmt+Rex,'-k') ;   hold on;
set(gca,'XDir','reverse'); ylabel('R_{1\rho}(\Delta\omega) [s^{-1}]'); set(gca,'yLim',[0 6]); 
yyaxis right
h_rex=plot(P.xZspec,Rex,'-') ;   hold on;
set(gca,'XDir','reverse'); ylabel('R_{ex}(\Delta\omega) [s^{-1}]'); set(gca,'yLim',[0 1]);
col_rex=get(h_rex,'Color');
hold on; legend({'R_{1\rho,wmt}','R_{1\rho}=R_{1\rho,wmt}+R_{ex}', 'R_{ex}'},'FontSize',8)

subplot(4,1,2),
plot(P.xZspec,Z,'-k',P.xZspec,Zref,'--') ;   hold on;
set(gca,'XDir','reverse'); ylabel('Z(\Delta\omega)'); set(gca,'yLim',[0 1]);
legend({'Z','Z_{ref}'},'FontSize',8)
subplot(4,1,3),
plot(P.xZspec,Zref-Z,'-','Color',col_rex) ; hold on;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('MTR_{LD}(\Delta\omega)'); set(gca,'yLim',[0 0.15]);
legend({'MTR_{LD}=Z_{ref}-Z'},'FontSize',8)

subplot(4,2,7), % vary B1 parameter
vary=[0.25  0.75 2 4 ]; % define value range for variation

for ii=1:numel(vary)
    P=Pstart;       % reset previous changes
    P.B1=vary(ii);  % define which parameter you want to vary
    
%     plot(P.xZspec,Z_cw(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
    Pref=Psim;    
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0; 
    plot(P.xZspec,Z_cw_yaml(P,Pref)-Z_cw_yaml(P,Psim),'-') ;   hold on;
  
end;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('MTR_{LD}(\Delta\omega)'); set(gca,'yLim',[0 0.15]);
set(gca,'xLim',[0 6]);
legend(strsplit(sprintf('B_1 = %.2f µT;',vary),';'),'FontSize',7)


subplot(8,2,14), % vary B1 parameter
%% alpha(B1,B0) verläufe  22feb
clear B1 Rexb1 Dilution MTRb1
B1=0.01:0.1:4;  
P=Pstart;
for ii=1:numel(B1)
P.B1=B1(ii);
P.xZspec=Psim.CESTPool(1).dw;
    Pref=Psim;    
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0; 
[Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw_yaml(P,Psim);
[Zref ,~, Rpwref, Rpmtref]=Z_cw_yaml(P,Pref);
Rexb1(ii)=Rex;
Dilution(ii)= (Zref)^2;
MTRb1(ii)=Zref-Z;
end
co = get(gca,'ColorOrder');
set(gca,'ColorOrder',circshift(co,3,1));
hold all
w1=Psim.Scanner.Gamma*B1;
alpha = w1.^2./(w1.^2+Psim.CESTPool(1).k*(Psim.CESTPool(1).k+Psim.CESTPool(1).R2));
plot(B1,alpha); hold on;
plot(B1,Dilution); hold on;
plot(B1,alpha.*Dilution); hold on;
plot(B1,MTRb1/(Psim.CESTPool(1).k*Psim.CESTPool(1).f)*R1obs); hold on;
legend({'\alpha','\sigma^\prime','\alpha\cdot\sigma^\prime','MTR\cdotR_{1}/k_sf_s'},'FontSize',7)

% Zrefsp = Rex* R1obs/(R1pwmt^2)          % with spillover
% Zref = Rex* R1obs/R1obs^2 =Rex/R1obs  % without spillover
% spillover term Zrefsp/Zref = R1obs^2/(R1pwmt^2) ~= Zref^2
% alpha from MTRLD: MTRLD = fb*kb*alpha * Zref^2 / R1obs
% alpha * Zref^2 = MTRLD / /fb*kb) *R1obs
set(gca,'XTicklabels',[]);
% xlabel('RF irradiation amplitude B_1 [µT]');
% ylabel('labeling efficiency');

subplot(8,2,16), % vary B1 parameter
%% alpha(B1,B0) verläufe  22feb
clear B1 Rexb1 Dilution MTRb1
B1=0.01:0.1:4;  
P=Pstart;
for ii=1:numel(B1)
P.B1=B1(ii);
P.xZspec=Psim.CESTPool(2).dw;
    Pref=Psim;    
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0; 
[Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw_yaml(P,Psim);
[Zref ,~, Rpwref, Rpmtref]=Z_cw_yaml(P,Pref);
Rexb1(ii)=Rex;
Dilution(ii)= (Zref)^2;
MTRb1(ii)=Zref-Z;
end
co = get(gca,'ColorOrder');
set(gca,'ColorOrder',circshift(co,3,1));
hold all
w1=Psim.Scanner.Gamma*B1;
alpha = w1.^2./(w1.^2+Psim.CESTPool(2).k*(Psim.CESTPool(2).k+Psim.CESTPool(2).R2));
plot(B1,alpha); hold on;
plot(B1,Dilution); hold on;
plot(B1,alpha.*Dilution); hold on;
plot(B1,MTRb1/(Psim.CESTPool(2).k*Psim.CESTPool(2).f)*R1obs); hold on;
% ylim([0,0.25]);
% legend({'\alpha','\sigma^\prime','\alpha\cdot\sigma^\prime','MTR\cdotR_{1obs}/k_sf_s'},'FontSize',7)

% Zrefsp = Rex* R1obs/(R1pwmt^2)          % with spillover
% Zref = Rex* R1obs/R1obs^2 =Rex/R1obs  % without spillover
% spillover term Zrefsp/Zref = R1obs^2/(R1pwmt^2) ~= Zref^2
% alpha from MTRLD: MTRLD = fb*kb*alpha * Zref^2 / R1obs
% alpa * Zref^2 = MTRLD / /fb*kb) *R1obs

xlabel('RF irradiation amplitude B_1 [µT]');
ylabel('labeling efficiency');
