%% analytic Z-spectra
%   Date: 2014/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  - moritz.zaiss@fau.de
%   CEST sources  Copyright (C) 2014  Moritz Zaiss
%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   **********************************
%
%   --SHORT  DOC--
%   The parameter struct P contains all system and sequence parameters.
%   Z_cw_2pool(P) calculates the Z-spectrum and returns it as a vector.
%
%   --references--:
%   2-pool-cw:      Zaiss and Bachert. NBM 2013;26(5):507�18. doi:10.1002/nbm.2887   and Zaiss et al JCP 2012;136:144106. doi:10.1063/1.3701178     
%   3-pool-cw:      Zaiss et al. NBM 2015 Feb;28(2):217-30. doi: 10.1002/nbm.3237.   
%   2-pool-pulsd SL: Roeloffs et al. NBM 2014; 28, 40�53, doi: 10.1002/nbm.3192.

%% SETUP
clearvars P Pref Pstart
clc
% setup pool system parameters
    %water pool 'a'
    P.R1A=1/3;          % longitudinal relaxation rate 1/T1 of pool a  [s^-1] 
    P.R2A=2;            % transversal relaxation rate 1/T2 of pool a  [s^-1] 
    P.dwA=0;            % chemical shift of the water pool in [ppm] 

    % CEST pool 'b'
    P.fB=0.0018018;     % proton fraction fB=M0B/M0A, e.g. 50mM creatine in water:  4*50/(2*55.5*1000)=0.0018018;
    P.kBA=50;            % exchange rate [s^-1]     % corresponds to creatine at ~ 22�C, pH=6.4 in PBS (Goerke et al.)
    P.dwB=1.9;          % chemical shift of the CEST pool in [ppm] 
    P.R2B=30;           % transversal relaxation rate 1/T2 of pool b  [s^-1] 
    P.R1B=1;            % longitudinal relaxation rate 1/T1 of pool b  [s^-1] 

    % % semi-solid MT pool 'c'
    P.fC=0.139;               % proton fraction of the semi-solid pool (WM-like)
    P.kCA=23;                 % exchange rate [s^-1]
    P.dwC=-2;                 % deltaW_B in [ppm} (chemical shift)
    P.R2C=1/(9.1*10^-6);      % transversal relaxation rate 1/T2 of pool c [s^-1]
    P.R1C=1;                  % longitudinal relaxation rate 1/T1 of pool c [s^-1]
    P.MT_lineshape='Lorentzian';            % semi-solid Lineshape

% setup sequence parameters
    P.Zi=1;                 % Z initial, in units of thermal M0, Hyperpol.: 10^4                  
    P.FREQ=300;             % static B0 field [MHz] ~7T ; ppm and �T are used for offsets and B1, therefore gamma=267.5153 is given in Hz.
    P.B1=1;                 % irradiation amplitude [�T]
    P.tp=80;                % pulse duration = saturation time [s]
    P.xZspec= [-30:0.1:30];   % chemical shift of the CEST pool in [ppm] 

    Pstart=P;
    
% PLOT Z_cw(P)
    
figure(1), plot(P.xZspec,Z_cw(P),'.-') ;   hold on;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)'); set(gca,'yLim',[0 1]);
text(0,0,evalc('P'),'FontSize',8)

% PLOT CEST:  Z_cw(Pref)-Z_cw(P) with reference struct Pref.
Pref=P;    
Pref.fB=0;                  % reference: Z of system without CEST pool
%Pref.xZspec=-P.xZspec;     % reference: Z of opposite frequency
plot(P.xZspec,Z_cw(Pref)-Z_cw(P),'r.-') ;   hold on; legend({'Z','Z_{ref}-Z'})

% PLOT relaxation-compensated CEST: AREX
% theta=atan(P.B1*267.5153./((P.xZspec-P.dwA)*2*pi*P.FREQ));
% plot(P.xZspec,cos(theta).^2.*(P.R1A+P.fC*P.R1C).*(1./Z_cw(P)-1./Z_cw(Pref)),'g.-') ;   hold on;
% legend({'Z','Z_{ref}-Z','AREX'})

%% vary a parameter
vary=[ 2 3 4 5 6]; % define value range for variation

for ii=1:numel(vary)
    P=Pstart;       % reset previous changes
    P.B1=vary(ii);  % define which parameter you want to vary
    
    figure(2), plot(P.xZspec,Z_cw(P),'.-','Color',cl(ii,numel(vary))) ;   hold on;
    
    Pref=P;    Pref.fB=0; %Pref.xZspec=-P.xZspec; %
    plot(P.xZspec,Z_cw(Pref)-Z_cw(P),'*-','Color',cl(ii,numel(vary))) ;   hold on;
  % PLOT relaxation-compensated CEST: AREX
  %  plot(P.xZspec,cos(theta).^2.*(P.R1A+P.fC*P.R1C).*(1./Z_cw(P)-1./Z_cw(Pref)),'d-','Color',cl(ii,numel(vary))) ;   hold on;
    
end;
set(gca,'XDir','reverse'); xlabel('\Delta\omega [ppm]'); ylabel('Z(\Delta\omega)'); set(gca,'yLim',[0 1]);
