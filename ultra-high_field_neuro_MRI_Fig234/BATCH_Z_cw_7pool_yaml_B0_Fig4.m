
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
set(gcf,'Position',[100   162   700   180]);
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

vary=[66 123 300 400 ]; % define value range for variation

for ii=1:numel(vary)
    
    P=Pstart;       % reset previous changes
    P.FREQ=vary(ii);  % define which parameter you want to vary
    
    %     plot(P.xZspec,Z_cw(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
    Pref=Psim;
    Pref.CESTPool(1).f=0;  Pref.CESTPool(2).f=0;  Pref.CESTPool(3).f=0; Pref.CESTPool(4).f=0; Pref.CESTPool(5).f=0;
    plot(P.xZspec,Z_cw_yaml(P,Pref)-Z_cw_yaml(P,Psim),'-') ;   hold on;
    
end
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('CESTR(\Delta\omega)'); set(gca,'yLim',[0 0.07]);
set(gca,'xLim',[-6 6]);
legend(strsplit(sprintf('B_0 = %.1f T;',[1.5 3 7 9.4]),';'),'FontSize',8)

