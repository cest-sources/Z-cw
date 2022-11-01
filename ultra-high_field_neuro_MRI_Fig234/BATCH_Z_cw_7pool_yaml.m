
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

figure(1),
set(gcf,'Position',[587   162   560   620]);
t = annotation('textbox','String','a)','Position',[0 0.85 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','b)','Position',[0 0.42 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';

set(groot,'defaultLineLineWidth',1.2)


bmsimfile='Stanisz_3T/WM_3T_Stanisz2005_5pool_bmsim.yaml';
% bmsimfile='Wang_3T/WM_3T_Wang2020_5pool_bmsim.yaml';
Psim = readSimulationParameters([pulseqCEST_simlib bmsimfile]);



% setup sequence parameters
P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4
P.FREQ=300;             % static B0 field [MHz] ~7T ; ppm and µT are used for offsets and B1, therefore gamma=267.5153 is given in Hz.
P.B1=0.75;                 % irradiation amplitude [µT]
P.tp=2;                % pulse duration = saturation time [s]
P.xZspec= [-6:0.05:6];   % chemical shift of the CEST pool in [ppm]

Pstart=P;

vary=[0.25 0.75 2 4 ]; % define value range for variation

for ii=1:numel(vary)
    
    % PLOT Z_cw(P)
    P=Pstart;       % reset previous changes
    P.B1=vary(ii);  % define which parameter you want to vary
    
    [Z ,Rex, Rpw, Rpmt,R1obs]=Z_cw_yaml(P,Psim);
    Pref=Psim;
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0;
    [Zref ,Rex_ref]=Z_cw_yaml(P,Pref);
    
    
    subplot(2,1,1),
    plot(P.xZspec,Z,'-k','Displayname',sprintf('Z_{lab}(B_1 = %.2f µT)',P.B1)); hold on;
    plot(P.xZspec,Zref,'--','Displayname',sprintf('Z_{ref}(B_1 = %.2f µT)',P.B1)) ;   
    set(gca,'XDir','reverse'); ylabel('Z(\Delta\omega)'); set(gca,'yLim',[0 1]);

    subplot(2,1,2),
    
    P=Pstart;       % reset previous changes
    P.B1=vary(ii);  % define which parameter you want to vary
    
    %     plot(P.xZspec,Z_cw(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
    Pref=Psim;
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0;
    plot(P.xZspec,Z_cw_yaml(P,Pref)-Z_cw_yaml(P,Psim),'-') ;   hold on;
    
end
subplot(2,1,1),
legend('Location','NorthEastOutside')
subplot(2,1,2),
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('CESTR(\Delta\omega)'); set(gca,'yLim',[0 0.07]);
set(gca,'xLim',[-6 6]);
legend(strsplit(sprintf('B_1 = %.2f µT;',vary),';'),'FontSize',8)
legend('Location','NorthEastOutside')
%% alpha(B1,B0) verläufe  22feb
figure(2),
set(gcf,'Position',[200   300   900   220]);
n=0
for kkk=[1 2 5]
    n=n+1;
    subplot(1,3,n), % vary B1 parameter
    
    clear B1 Rexb1 Dilution MTRb1
    B1=0.01:0.1:6;
    P=Pstart;
    for ii=1:numel(B1)
        P.B1=B1(ii);
        P.xZspec=Psim.CESTPool(kkk).dw;
        Pref=Psim;
        Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0;
        [Zref]=Z_cw_yaml(P,Pref);
        Dilution(ii)= Zref^2;
    end
    co = get(gca,'ColorOrder');
    set(gca,'ColorOrder',circshift(co,3,1));
    hold all
    w1=Psim.Scanner.Gamma*B1;
    alpha = w1.^2./(w1.^2+Psim.CESTPool(kkk).k*(Psim.CESTPool(kkk).k+Psim.CESTPool(kkk).R2));
    plot(B1,alpha); hold on;
    plot(B1,Dilution); hold on;
    plot(B1,alpha.*Dilution); hold on;
    legend({'\alpha','\sigma^\prime','\alpha\cdot\sigma^\prime'},'FontSize',10)
    xlabel('RF irradiation amplitude B_1 [µT]');
    if kkk==1
        ylabel({'labeling efficiency \alpha','spillover dilution \sigma'''});
    end
end