% 2pool cw-solution for Z-spectra
%   Date: 2014/08/01 
%   Version for CEST-sources.de
%   Author: Moritz Zaiss  - m.zaiss@dkfz.de
%   Divison of Medical Physics in Radiology
%   German Cancer Research Center (DKFZ), Heidelberg, Germany , http://www.dkfz.de/en/index.html
%   CEST sources - Copyright (C) 2014  Moritz Zaiss

%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   CEST sources  Copyright (C) 2014  Moritz Zaiss
%   **********************************
%
% P= struct of system and sequence parameters
% xZspec = frequency offsets in ppm
% returns vector Z(xZspec)

function [Z ,Rex_n, Reff, Rex_MT, R1obs]=Z_cw_yaml(P,Psim)

w_ref=2*pi*P.FREQ;                  % reference frequency in rad/s
gamma=Psim.Scanner.Gamma;                     % for protons
w1 =P.B1*gamma;                     % Rabi frequency in rad/s
da=(P.xZspec)*w_ref;          % freq. offset in rad/s
theta=atan(w1./da);                 % theta in rad
checkP(P,Psim);                          % check if conditions are fullfilled

%% The individual relaxation rates in the rotating frame (R1p)

%Reff =  R1p of free water
Reff=Psim.WaterPool.R1*cos(theta).^2 +Psim.WaterPool.R2*sin(theta).^2 ;

%Rex, exchange dependent relaxation in the rotating frame
Rex_n=zeros(1,numel(P.xZspec));
if isfield(Psim,'CESTPool')
    for ii=1:numel(Psim.CESTPool)  
        Rex{ii}=Rex_Lorentz(da,w1,(P.xZspec-Psim.CESTPool(ii).dw)*w_ref,Psim.CESTPool(ii).f,Psim.CESTPool(ii).k,Psim.CESTPool(ii).R2);
        Rex_n=Rex_n + Rex{ii};
    end
end

%Rex_MT, exchange dependent relaxation in the rotating frame of the MT pool
Rex_MT=0; R1obs=Psim.WaterPool.R1;
if isfield(Psim,'MTPool')
    if Psim.MTPool.f>0
        [Rex_MT, R1obs]=Rex_MT_1(da,w1,Psim.WaterPool.R1,Psim.WaterPool.R2,(P.xZspec-Psim.MTPool.dw)*w_ref,Psim.MTPool.f,Psim.MTPool.k,Psim.MTPool.R1,Psim.MTPool.R2,Psim.MTPool.Lineshape);
    end
end
%% Superposition of R1p
% this line creates the multi-pool system
%       water +     MT   +   CEST
R1p  =  Reff +   Rex_MT +   Rex_n/(1+P.fC); % altered according to NBM 2014


Pz=cos(theta);    % for cw   ; Pz=1;    % for SL
Pzeff=cos(theta); % for cw   ; Pzeff=1; % for SL

%% The R1p-model for the steady-state Z-spectrum
Zss=cos(theta).*Pz.*R1obs./R1p;  % altered according to NBM 2014, R1a->R1obs

%% The R1p-model for the Z-spectrum after rf irradtaion of duration P.tp
Z= (Pz.*Pzeff.*P.Zi -  Zss).*exp(-(R1p.*P.tp)) +Zss;

end


%HyperCESTLimit %JCP paper;   assumes (R1B<<kBA, Reff<<R2B)
function Rex=Rex_Hyper_full(da,w1,db,fb,kb,r2b)
ka=kb*fb;
Rex=((ka.*kb.*w1.^2.*((-da+db).^2 + (r2b.*(da.^2 + (ka + kb).^2 + kb.*r2b + w1.^2))./kb))./...
    ((ka + kb).*(db.^2.*w1.^2 + ka.*r2b.*w1.^2) + ...
    (ka + kb).*((da.*db - ka.*r2b).^2 + (db.*ka + da.*(kb + r2b)).^2 + ...
    (ka + kb + r2b).^2.*w1.^2) + (ka + kb + r2b).*(da.^2.*w1.^2 + w1.^4)));
end

%LorentzLimit %NBM paper;    assumes (kAB<<kBA, R1B<<kBA, Reff<<R2B)
function Rex=Rex_Lorentz(da,w1,db,fb,kb,r2b)
ka=kb*fb;

REXMAX= ka.*w1^2./(da.^2+w1.^2).*((da-db).^2 +(da.^2+w1.^2).*r2b./kb + r2b.*(kb+r2b));
GAMMA=2*sqrt( (kb+r2b)./kb.*w1.^2 + (kb+r2b).^2);
Rex=REXMAX./((GAMMA./2).^2+db.^2);
end


function [Rex_MT, R1obs]=Rex_MT_1(da,w1,R1A,R2A,dc,fc,kca,R1C,R2C,MT_lineshape)

    theta   = atan(w1./da);
    Reff=R1A*cos(theta).^2 +R2A*sin(theta).^2;
    kac=kca*fc;
    
    
    r1a     = R1A-Reff;
    r2a     = R2A-Reff;
    r1c     = R1C-Reff;
    rfmt    = RF_MT(1/R2C,w1,dc,MT_lineshape);
    
    Rex_MT =    (((da.^2 + r2a.^2).*(kca.*r1a + (kac + r1a).*(r1c + rfmt)) + ...
        r2a.*(kca + r1c + rfmt).*w1.^2)./(da.^2.*(kac + kca + r1a + r1c + rfmt) + ...
        r2a.*(kca.*(2.*r1a + r2a) + r2a.*(r1c + rfmt) + ...
        kac.*(2.*r1c + r2a + 2.*rfmt) + ...
        r1a.*(2.*r1c + r2a + 2.*rfmt)) + (kca + r1c + r2a + rfmt).*w1.^2));


R1obs =    0.5*( kac + kca + R1A+ R1C - sqrt(( kac + kca + R1A + R1C )^2 - 4*( kca*R1A + kac*R1C + R1A*R1C )));
% R1obs=(P.R1A+P.fC*P.R1C)/(1+P.fC);  % approximation if R1<<k 
end



function rfmt=RF_MT(T2c,w1,dw,lineshape)

if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
    step=1/10000;
    sup_theta = 0:step:pi/2;
    
    cutoff=10;
    for j=1:length(dw)
        if abs(dw(j)) >= cutoff  %% the superlorentz has a pole: this avoids infinities
            % see Morrsion and Henkelman 1995.  units = s.  Seems weird to
            % me.  Need to multiply by w1^2 * pi to get saturation rate.
            du=.0001;
            u=0:du:1;
            integrand2=sqrt(2/pi)*T2c./abs(3*u.^2-1) .* exp(-2*(dw(j)*T2c./abs(3*u.^2-1)).^2);
            G(j)=w1.^2.*pi.*sum(integrand2)*du;
        else
            X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
            Y = RF_MT(T2c,w1,X,lineshape);
            G(j)=interp1(X,Y,dw(j),'spline');
        end
        rfmt=G;
        
    end

elseif strcmp(lineshape,'Gaussian')
    rfmt=w1.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); %Gaussian
else
        rfmt=w1.^2*T2c./(1+(dw.*T2c).^2); %Lorentzian
end;

end

function checkP(P,Psim)
if P.tp<1/Psim.WaterPool.R2
    warning('Assumption of single eigenvalue description might be violated')
end

if isfield(Psim,'CESTPool')
    if Psim.CESTPool(1).f>0.2
        warning('Assumption of small CEST pool might be violated')
    end


    if Psim.CESTPool(1).k<5*Psim.WaterPool.R1 || Psim.CESTPool(1).k<5*Psim.CESTPool(1).R1
        warning('Assumption k>>R1 might be violated')
    end
end

end