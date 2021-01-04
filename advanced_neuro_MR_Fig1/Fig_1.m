
%% analytic Z-spectra
%   Date: 2014/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  - moritz.zaiss@uk-erlangen.de
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
% parameters taken from Table 1 in the same chapter

%% SETUP
clearvars P Pref Pstart
clc
% setup pool system parameters 7T  (https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25677)
    %water pool 'a'
    P.R1A=1/1.334;          % longitudinal relaxation rate 1/T1 of pool a  [s^-1] 
    P.R2A=1/(57e-3);            % transversal relaxation rate 1/T2 of pool a  [s^-1] 
    P.dwA=0;            % chemical shift of the water pool in [ppm] 

    % CEST pool 'b'
    P.fB=0.00156;     % proton fraction fB=M0B/M0A, WM value from Table 1  (e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;)
    P.kBA=165;            % exchange rate [s^-1]     WM value from Table 1 
    P.dwB=3.5;          % chemical shift of the CEST pool in [ppm] 
    P.R2B=73.1;           % transversal relaxation rate 1/T2 of pool b  [s^-1]  1/0.01368s, WM value from Table 1 
    P.R1B=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 

    % % semi-solid MT pool 'c'
    P.fC=0.11;               % proton fraction of the semi-solid pool (WM-like) WM value from Table 1 
    P.kCA=25;                 % exchange rate [s^-1] WM value from Table 1 
    P.dwC=0;                 % deltaW_B in [ppm} (chemical shift)
    P.R2C=1/(20.48*10^-6);      % transversal relaxation rate 1/T2 of pool c [s^-1] WM value from Table 1 
    P.R1C=1;                  % longitudinal relaxation rate 1/T1 of pool c [s^-1]
    P.MT_lineshape='SuperLorentzian';            % semi-solid Lineshape

% setup sequence parameters
    P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
    P.FREQ=300;             % static B0 field [MHz] ~7T ; ppm and µT are used for offsets and B1, therefore gamma=267.5153 is given in Hz.
    P.B1=1;                 % irradiation amplitude [µT]
    P.tp=5;                % pulse duration = saturation time [s]
    P.xZspec= [-6:0.05:6];   % chemical shift of the CEST pool in [ppm] 

    Pstart=P;
    
% PLOT Z_cw(P)
[Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw(P);
Pref=P;    
Pref.fB=0;  
[Zref ,Rex_ref]=Z_cw(Pref);
    
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
set(gca,'XDir','reverse'); ylabel('R_{1\rho}(\Delta\omega) [s^{-1}]'); set(gca,'yLim',[0 3]); 
yyaxis right
h_rex=plot(P.xZspec,Rex,'-') ;   hold on;
set(gca,'XDir','reverse'); ylabel('R_{ex}(\Delta\omega) [s^{-1}]'); set(gca,'yLim',[0 0.5]);
col_rex=get(h_rex,'Color');
hold on; legend({'R_{1\rho,wmt}','R_{1\rho}=R_{1\rho,wmt}+R_{ex}', 'R_{ex}'},'FontSize',8)

subplot(4,1,2),
plot(P.xZspec,Z,'-k',P.xZspec,Zref,'--') ;   hold on;
set(gca,'XDir','reverse'); ylabel('Z(\Delta\omega)'); set(gca,'yLim',[0 1]);
legend({'Z','Z_{ref}'},'FontSize',8)
subplot(4,1,3),
plot(P.xZspec,Zref-Z,'-','Color',col_rex) ; 
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('MTR_{LD}(\Delta\omega)'); set(gca,'yLim',[0 0.15]);
legend({'MTR_{LD}=Z_{ref}-Z'},'FontSize',8)

subplot(4,2,7), % vary B1 parameter
vary=[0.25  1 2 ]; % define value range for variation

for ii=1:numel(vary)
    P=Pstart;       % reset previous changes
    P.B1=vary(ii);  % define which parameter you want to vary
    
%     plot(P.xZspec,Z_cw(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
    
    Pref=P;    Pref.fB=0; %Pref.xZspec=-P.xZspec; %
    plot(P.xZspec,Z_cw(Pref)-Z_cw(P),'-') ;   hold on;
  
end;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('MTR_{LD}(\Delta\omega)'); set(gca,'yLim',[0 0.15]);
set(gca,'xLim',[0 6]);
legend(strsplit(sprintf('B_1 = %.1f µT;',vary),';'),'FontSize',7)


subplot(4,2,8), % vary B1 parameter
%% alpha(B1,B0) verläufe  22feb
clear B1 Rexb1 Dilution MTRb1
B1=0.01:0.1:4;  
P=Pstart;
for ii=1:numel(B1)
P.xZspec=P.dwB;
P.B1=B1(ii);
Pref=P; Pref.fB=0;
[Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw(P);
[Zref ,~, Rpwref, Rpmtref]=Z_cw(Pref);
Rexb1(ii)=Rex;
Dilution(ii)= (Zref)^2;
MTRb1(ii)=Zref-Z;
end
co = get(gca,'ColorOrder');
set(gca,'ColorOrder',circshift(co,3,1));
hold all
plot(B1,Rexb1./(P.kBA*P.fB)); hold on;
plot(B1,Dilution); hold on;
plot(B1,Rexb1./(P.kBA*P.fB).*Dilution); hold on;
plot(B1,MTRb1/(P.kBA*P.fB)*R1obs); hold on;
legend({'\alpha','\sigma^\prime','\alpha\cdot\sigma^\prime','MTR\cdotR_{1obs}/k_sf_s'},'FontSize',7)

% Zrefsp = Rex* R1obs/(R1pwmt^2)          % with spillover
% Zref = Rex* R1obs/R1obs^2 =Rex/R1obs  % without spillover
% spillover term Zrefsp/Zref = R1obs^2/(R1pwmt^2) ~= Zref^2
% alpha from MTRLD: MTRLD = fb*kb*alpha * Zref^2 / R1obs
% alpa * Zref^2 = MTRLD / /fb*kb) *R1obs

xlabel('RF irradiation amplitude B_1 [µT]');
ylabel('labeling efficiency');
