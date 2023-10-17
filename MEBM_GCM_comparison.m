% This code compares MEBM output with either AM2, CAM4, or CAM3 output
clear all; close all; clc;

dh_gms = 1.06; 

% define some constants
lv = 2.5*10^6; % latent heat of vaporization [J / kg] - this has a temperature dependence
ls = 2.834*10^6; % latent heat of sublimation
lf = ls-lv; % latent heat of fusion
rho_w = 1000; % density of water (kg / m**3)
jmx=101; 
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:); % sine latitude
phi = asin(x)*180/pi;                        % latitude

% weighting function
wt = normpdf(x,0,sind(15));
wt = 1-wt;

% choose GCM to compare
model_flg = 'CAM4';
if sum(strcmp(model_flg,{'AM2','CAM4','CAM3'}))<1
    warning('please enter a valid model flag')
    return;
end

dir = ['/Users/krist/OneDrive/Desktop/Rose Models/' model_flg];

% define some constants
lv = 2.5*10^6; % latent heat of vaporization [J / kg] - this has a temperature dependence
ls = 2.834*10^6; % latent heat of sublimation
lf = ls-lv; % latent heat of fusion
rho_w = 1000; % density of water (kg / m**3)
jmx=101; 
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:); % sine latitude
phi = asin(x)*180/pi;                        % latitude

qH0 = ((-299*2)/(90*cos(2*pi/9)))*sin((18/5)*(abs(deg2rad(phi))-(2*pi/9)));
qH = min(0,qH0);

qT0 = ((-16*2)/(3*sqrt(3)))*cos(3*deg2rad(phi));
qT = min(0,qT0);


% load files in directory

if strcmp(model_flg,'AM2')
    % You can look at the file info using "ncdisp":
    %ncdisp([dir '/2upH+CO2_1992_1996.atmos_month.zon.nc']);
    
    upH = '2upH_1994_2003.atmos_month.zon.nc';
    upT = '2upT_1994_2003.atmos_month.zon.nc';
    ctrl = 'ctrl_1994_2003.atmos_month.zon.nc';

    % load the necessary variables to calculate precipitation
    filenm = {[dir '/' upH],[dir '/' upT],[dir '/' ctrl]};

    for i = 1:3
        lat_AM2 = ncread(filenm{i},'lat');
        lon = ncread(filenm{i},'lon');
        time = ncread(filenm{i},'time');
        evap = ncread(filenm{i},'evap');
        precip = ncread(filenm{i},'precip');
        x = sin(deg2rad(lat_AM2));
        wt_AM2 = cos(deg2rad(lat_AM2)); 
        surf_t = ncread(filenm{i},'t_surf'); % surface temp
    

        netrad = ncread(filenm{i},'netrad_toa');

        % save variables
        if i==1
            evap_lat_upH= evap;
            precip_lat_upH= precip;
            netrad_upH= netrad;
            surf_t_upH= surf_t;
        elseif i==2
            evap_lat_upT= evap;
            precip_lat_upT= precip;
            netrad_upT= netrad;
            surf_t_upT= surf_t;
        else
            evap_lat_ctrl= evap;
            precip_lat_ctrl= precip;
            netrad_ctrl = netrad;
            surf_t_ctrl = surf_t;
        end
    end

elseif strcmp(model_flg,'CAM3')
    
    upH = 'QAqu_2upH.cam2.h0.zonclim.nc';
    upT = 'QAqu_2upT.cam2.h0.zonclim.nc';
    ctrl = 'QAqu_SNOWTEST.cam2.h0.zonclim.nc';

    % load the necessary variables to calculate precipitation
    filenm = {[dir '/' upH],[dir '/' upT],[dir '/' ctrl]};

    lat_CAM3 = ncread(filenm{1},'lat');
    wt_CAM3 = cos(deg2rad(lat_CAM3));

    for i = 1:3
        precip = ncread(filenm{i},'PRECT'); % total precipitation
        evap = ncread(filenm{i},'LHFLX'); % latent heat flux (evap)
        surf_t = ncread(filenm{i},'TS'); % surface temp
        swcf = ncread(filenm{i},'SWCF');

        precip = rho_w*lv*precip;

        net_sw = ncread(filenm{i},'FSNT'); % net solar radiation at top of model
        net_lw = ncread(filenm{i},'FLNT'); % net longwave at top of model

        netrad = net_sw-net_lw;

        % save variables
        if i==1
            evap_lat_upH = evap;
            precip_lat_upH = precip;
            netrad_upH = netrad;
            surf_t_upH = surf_t;
            swcf_upH = swcf;
            
        elseif i==2
            evap_lat_upT = evap;
            precip_lat_upT = precip;
            netrad_upT = netrad;
            surf_t_upT = surf_t;
            swcf_upT = swcf;
        else
            evap_lat_ctrl = evap;
            precip_lat_ctrl = precip;
            netrad_ctrl = netrad;
            surf_t_ctrl = surf_t;
            swcf_ctrl = swcf;
        end
    end

elseif strcmp(model_flg,'CAM4')

    upH = 'QAqu_qupH.cam.h0.zonclim.nc';
    upT = 'QAqu_qupT.cam.h0.zonclim.nc';
    ctrl = 'QAqu_ctrl.cam.h0.zonclim.nc';

    % load the necessary variables to calculate precipitation
    filenm = {[dir '/' upH],[dir '/' upT],[dir '/' ctrl]};
    lat_CAM4 = ncread(filenm{1},'lat');
    wt_CAM4 = cos(deg2rad(lat_CAM4));


    for i = 1:3
        precipC = ncread(filenm{i},'PRECC'); % convective precip
        precipL = ncread(filenm{i},'PRECL'); % large-scale precip
        precipSC = ncread(filenm{i},'PRECSC'); % convective snow
        precipSL = ncread(filenm{i},'PRECSL'); % large-scale snow
        surf_t = ncread(filenm{i},'TS'); % surface temp
        precip = precipC + precipL + precipSC + precipSL;
        evap = ncread(filenm{i},'LHFLX'); % latent heat flux (evap)
        swcf = ncread(filenm{i},'SWCF');

        precip = lv*rho_w*precip;

        net_sw = ncread(filenm{i},'FSNT'); % net solar radiation at top of model
        net_lw = ncread(filenm{i},'FLNT'); % net longwave at top of model

        netrad = net_sw-net_lw;
        
        % save variables
        if i==1
            evap_lat_upH = evap;
            precip_lat_upH = precip;
            netrad_upH = netrad;
            surf_t_upH = surf_t;
            swcf_upH = swcf;
        elseif i==2
            evap_lat_upT = evap;
            precip_lat_upT = precip;
            netrad_upT = netrad;
            surf_t_upT = surf_t;
            swcf_upT = swcf;
        else
            evap_lat_ctrl = evap;
            precip_lat_ctrl = precip;
            netrad_ctrl= netrad;
            surf_t_ctrl= surf_t;
            swcf_ctrl = swcf;
        end
    end

end

% calculate lambda, forcing
if strcmp(model_flg,'AM2')
    AM2 = load_rose14('AM2');
    lat = AM2.ctrl.lat;
    lam_upT = (-AM2.upT.netrad+AM2.ctrl.netrad)./(-AM2.upT.temp+AM2.ctrl.temp); % can separate into SW or LW also
    lam_upH = (-AM2.upH.netrad+AM2.ctrl.netrad)./(-AM2.upH.temp+AM2.ctrl.temp);
    lam_upT = -interp1(lat,lam_upT,phi,'linear');
    lam_upH = -interp1(lat,lam_upH,phi,'linear');
    anom_upT = netrad_ctrl-netrad_upT;
    anom_upH = netrad_ctrl-netrad_upH;
    surf_t_anom_upT = -surf_t_ctrl+surf_t_upT;
    surf_t_anom_upH = -surf_t_ctrl+surf_t_upH;
    p_anom_upT = -precip_lat_ctrl+precip_lat_upT;
    p_anom_upH = -precip_lat_ctrl+precip_lat_upH;
    e_anom_upT = -evap_lat_ctrl+evap_lat_upT;
    e_anom_upH = -evap_lat_ctrl+evap_lat_upH;
    e_m_p_upT = evap_lat_upT - precip_lat_upT;
    e_m_p_upH = evap_lat_upH - precip_lat_upH;
    e_m_p_anom_upT = e_m_p_upT - (evap_lat_ctrl - precip_lat_ctrl);
    e_m_p_anom_upH = e_m_p_upH- (evap_lat_ctrl - precip_lat_ctrl);
    forcing_upT = interp1(lat,anom_upT,phi,'linear');
    forcing_upH = interp1(lat,anom_upH,phi,'linear');
end
if strcmp(model_flg,'CAM3')
    CAM3 = load_rose14('CAM3');
    lat = CAM3.ctrl.lat;
    lam_upT = (-CAM3.upT.netrad+CAM3.ctrl.netrad)./(-CAM3.upT.temp+CAM3.ctrl.temp); % can separate into SW or LW also
    lam_upH = (-CAM3.upH.netrad+CAM3.ctrl.netrad)./(-CAM3.upH.temp+CAM3.ctrl.temp);
    lam_upT = -interp1(lat,lam_upT,phi,'linear');
    lam_upH = -interp1(lat,lam_upH,phi,'linear');
    anom_upT = netrad_ctrl-netrad_upT;
    anom_upH = netrad_ctrl-netrad_upH;
    surf_t_anom_upT = -surf_t_ctrl+surf_t_upT;
    surf_t_anom_upH = -surf_t_ctrl+surf_t_upH;
    p_anom_upT = -precip_lat_ctrl+precip_lat_upT;
    p_anom_upH = -precip_lat_ctrl+precip_lat_upH;
    e_anom_upT = -evap_lat_ctrl+evap_lat_upT;
    e_anom_upH = -evap_lat_ctrl+evap_lat_upH;
    e_m_p_upT = evap_lat_upT - precip_lat_upT;
    e_m_p_upH = evap_lat_upH - precip_lat_upH;
    e_m_p_anom_upT = e_m_p_upT - (evap_lat_ctrl - precip_lat_ctrl);
    e_m_p_anom_upH = e_m_p_upH- (evap_lat_ctrl - precip_lat_ctrl);
    forcing_upT = interp1(lat,anom_upT,phi,'linear');
    forcing_upH = interp1(lat,anom_upH,phi,'linear');
end
if strcmp(model_flg,'CAM4')
    CAM4 = load_rose14('CAM4');
    lat = CAM4.ctrl.lat;
    lam_upT = (-CAM4.upT.netrad+CAM4.ctrl.netrad)./(-CAM4.upT.temp+CAM4.ctrl.temp); % can separate into SW or LW also
    lam_upH = (-CAM4.upH.netrad+CAM4.ctrl.netrad)./(-CAM4.upH.temp+CAM4.ctrl.temp);
    lam_upT = -interp1(lat,lam_upT,phi,'linear');
    lam_upH = -interp1(lat,lam_upH,phi,'linear');
    anom_upT = netrad_ctrl-netrad_upT;
    anom_upH = netrad_ctrl-netrad_upH;
    surf_t_ctrl = surf_t_ctrl;
    surf_t_anom_upT = -surf_t_ctrl+surf_t_upT;
    surf_t_anom_upH = -surf_t_ctrl+surf_t_upH;
    p_anom_upT = -precip_lat_ctrl+precip_lat_upT;
    p_anom_upH = -precip_lat_ctrl+precip_lat_upH;
    e_anom_upT = -evap_lat_ctrl+evap_lat_upT;
    e_anom_upH = -evap_lat_ctrl+evap_lat_upH;
    e_m_p_upT = evap_lat_upT - precip_lat_upT;
    e_m_p_upH = evap_lat_upH - precip_lat_upH;
    e_m_p_anom_upT = e_m_p_upT - (evap_lat_ctrl - precip_lat_ctrl);
    e_m_p_anom_upH = e_m_p_upH- (evap_lat_ctrl - precip_lat_ctrl);
    forcing_upT = interp1(lat,anom_upT,phi,'linear');
    forcing_upH = interp1(lat,anom_upH,phi,'linear');
end


% Code for Moist/Dry Energy Balance Model
% Using old Crank Nicolson algorithm

% This code solves the energy balance model used in Roe et al. (Nat. Geosci., 2015)
% The model operates in perturbation mode (i.e., the anomaly with respect to a prescribed climatology):
% You can specify:-
% the pattern of climate forcing,Rf
% the pattern of climate feedbacks, -B
% the pattern of diffusivity, D
% whether you diffuse moist static energy, or just sensible heat

%time step in fraction of year
delt=1./50000; disp(['delt = ' num2str(delt)])
%delt=.1;
NMAX=60000; disp(['NMAX = ' num2str(NMAX)])

%set up x array (latitude).
jmx=101; 
%jmx = 101; disp(['jmx = ' num2str(jmx)])
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:);
phi = asin(x)*180/pi;

% climate parameters
A=203.3; % Size of longwave cooling constant [W/m2]
B=2.09*ones(size(x)); B=B(:); disp('B is constant for control')
% longwave cooling [W/(m2 K)]

% magnitude of diffusivity [units?]
Dmag = 0.2598; disp(['D = ' num2str(Dmag) ' W/(m2 K)'])% D = 0.2598 W/(m2 K) is the value used by TF10 
D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
 
% I think this C = rho * c * h_ml /(pi*1e7).
% this is consistent with a ~1m layer of dirt
% Note - heat capacity over LAND LAND LAND.
Cl = 0.2; disp(['Cl = ' num2str(Cl) ' J /(m2 K)'])

% Moisture parameters
relhum = 0.8;               % relative humidity
eps = 0.622;                % moisture coonstant
psfc = 9.8e4;               % (Pa)
e0 = 611.2;                 % vap. press (Pa)
a = 17.67; b = 243.5;       % sat vap constants !!T must be in temperature

L = 2.45e6;                 % latent heat of vaporization (J kg-1)
cp = 1004;                  % (J kg-1 K-1)
Re = 6.4e6;                 % [m] earth's radius
rho = 1e3;                  % [kg m-3] density

%-------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%For dry energy balance model, uncomment these values
% Dmag = 0.44; % magnitude of diffusivity [W/(m2 K)]
% D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
% relhum = 0;  % switch off humidity
% disp('Diffusing sensible heat only')
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%-------------------------------------------------
 
% Use setupfast to create a matrix that calculates D*d/dx[ (1-x^2)d/dx] of
% whatever it operates on.
 [~,Mdiv]=setupfastM(delx,jmx,D,0.,1.0,delt);
 

% this work ok. GHR mar 3,17
% load T_ctrl_FromGHR.mat;                        % control temp from Nick Siler, based on ERA
% T_ctrl = interp1(x_ghr,T_ghr,x,'pchip');  % interpolate to our grid
% T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl));            % symmetrize (optional)
T_ctrl = interp1(lat,surf_t_ctrl,phi,'linear')-273.15;
%T_ctrl = movmean(T_ctrl,10);
A1 = legendrefit(T_ctrl,12,'inv');
clear P; for i0 = 1:12+1; P(i0,:) = legendreP(i0-1,x); end
T_ctrl = A1'*P;
T_ctrl = T_ctrl';

q_ctrl = eps*relhum/psfc*e0*exp(a*T_ctrl./(b+T_ctrl)); q_ctrl=q_ctrl(:);% here T is in oC. q is g kg-1
theta_e_ctrl = 1/cp*(cp*(T_ctrl+273.15) + L*q_ctrl); % theta is mse divided by cp; note units of Kelvin are needed!!!

%%
%--------------------------------------------------------------
% Now do climate change experiment
disp('**Doing climate experiment**')

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% For flat forcing, uniform feedbacks, pick these values for R_frc and B
% R_frc = 4.0; % forcing in [W/m2]
% B = 2.09*ones(size(x)); B=B(:); disp('B is constant for expt.')

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% To use feedback patterns from an aquaplanet model (Feldl and Roe, GRL,
% 2013), uncomment these values
%%

B = lam_upT;
R_frc = qT;
% load NicoleNewH2O.mat
% lam_tot = 0.5*(lam_tot + fliplr(lam_tot));
% B = -interp1(sind(lat),lam_tot,x,'pchip')-0.5; % note nudge of 0.5 GHR jan21
% R_frc = interp1(sind(lat),delRf,x,'pchip');


%--------------------------------------------------------------
%set up inital T profile 
T = 0.5*(1-1*x.^2); T = T(:);
%load T_final
Tinit=T;
Tglob = mean(Tinit);

% Timestepping loop
for n=1:NMAX
   Tglob_prev = Tglob;
       
% Calculate src for this loop.
   src = R_frc/Cl; src=src(:);
   
% spec. hum, and theta_e
   q = eps*relhum/psfc*e0*exp(a*(T+T_ctrl)./(b+(T+T_ctrl))); q=q(:);% here T is in oC. q is g kg-1
   theta_e = 1/cp*(cp*((T+T_ctrl)+273.15) + L*q); % note units of Kelvin are needed!!!
   theta_e_pert = theta_e-theta_e_ctrl; %anomalous theta_e

% Calculate new T.
% Diffuse sensible heat energy (T)
% Diffusion of T anomaly
%   dT = delt/Cl*(R_frc - (B.*T) + Mdiv*T); if n == 1,disp(['Diffusing T anomaly']),end;

% Diffuse moist static energy (theta_e)
% diffusion of theta_e_anomaly
  dT = delt/Cl*(R_frc - (B.*T) + Mdiv*theta_e_pert); if n == 1, disp(['Diffusing theta_e anomaly']), end;

   T = T + dT; 
   
% Check to see if global mean temperature has converged
   Tglob=mean(T);
   Tchange = Tglob-Tglob_prev;
   if (abs(Tchange) < 1.0e-12), break; end
end

%divF_expt = -Mdiv*T;
q_expt = q;
q_pert = q_expt-q_ctrl;
divF_pert = -Mdiv*theta_e_pert;
h_pert = theta_e_pert*cp;
T_pert = T;
T_expt = T+T_ctrl;
src_pert = R_frc;

%  
%% Module to calculate terms in hyrological cycle
%define some things from the climatology first
idx = find(x>=0);
x15 = sind(15);xw = 0.15;
% wt = 0.5*(1+tanh((x(idx)-x15)/xw));     % weighting function for eddy fluxes
% wt = [flipud(wt); wt(2:end)]; %2-sided climate response
 
% some things about the control climate
Dhlf = 0.5*(D(1:end-1)+D(2:end));                           % [W m-2 K-1] calculate D on same grid as T,q, etc
h_ctrl = cp*(T_ctrl+273.15)+L*q_ctrl;                       % [J kg-1] control climate mse
F_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_ctrl,x);  % [W] control climate flux 
F_lh_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(L*q_ctrl,x);    % [W] control latent heat flux
F_hc_ctrl = (1-wt).*F_ctrl;                                 % [W] Hadley Cell Flux
% dh_gms = 1.5e4;                                             % [J kg-1] gross moist stability
heq_ctrl = h_ctrl(x==0);                                    % [J kg-1] moist static energy at the surface
V_ctrl = F_hc_ctrl./(heq_ctrl*dh_gms-h_ctrl);               % [kg s-1] Diagnosed mass transport in Hadley Cell
F_LH_ctrl = -L*V_ctrl.*q_ctrl + wt.*F_lh_ctrl;              % [W] latent heat (Hadley plus eddy)
divF_LH_ctrl = 1/(2*pi*Re^2)*gradient(F_LH_ctrl,x);         % [W m-2] divergence of latent heat flux
E_m_P_ctrl = divF_LH_ctrl;                   % [m yr-1]E-P

% code from Nick Siler to partition E-P;
alpha=L./461.5./(T_ctrl+273.15).^2;                         % Nick's alpha parameter
beta=cp/L./alpha./q_ctrl;                                   % beta parameter
RG=180*(1*(1-x.^2)-.4*exp(-(x./.15).^2));                   % [W m-2] idealized R-G
%RG=170*(1*(1-x.^2));%-.4*exp(-(x./.15).^2))
Ch=1.5e-3;                                                  % drag coefficient
LWfb=0;                                                     % LW feedback at surface, in W/m2/K
u=4+abs(sin(pi*x./1.5))*4;%-2.5*cos(3*asin(x));             % wind speed
rho_air=1.2;%psfc./287./(T_ctrl+273.15);                        % air density
E_ctrl=(RG.*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q_ctrl);
P_ctrl=E_ctrl-divF_LH_ctrl;


% some things about the new climate state
h_expt = h_ctrl+h_pert;                                 % [J kg-1] new total mse
F_expt = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_expt,x);      % [W] new total flux
F_lh_expt = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(L*q_expt,x);    % [W] control latent heat flux
F_hc_expt = (1-wt).*F_expt;                             % [W] Hadley Cell Flux
% dh_gms = ;                                           % [J kg-1] gross moist stability (unchanged)
heq_expt = h_expt(x==0);                                % [J kg-1] moist static energy at the surface
V_expt = F_hc_expt./(heq_expt*dh_gms-h_expt);           % [wkg s-1] Diagnosed mass transport in Hadley Cell
F_LH_expt = -L*V_expt.*q_expt + wt.*F_lh_expt;          % [W] latent heat (Hadley plus eddy)
divF_LH_expt = 1/(2*pi*Re^2)*gradient(F_LH_expt,x);     % [W m-2]
E_m_P_expt = divF_LH_expt;               % [m yr-1] E-P

alpha=L./461.5./(T_expt+273.15).^2;                         % Nick's alpha parameter
beta=cp/L./alpha./q_expt;                                   % beta parameter
E_expt=(RG.*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q_expt);
P_expt=E_expt-divF_LH_expt;
 

% P and E in the new climate state
%alpha=L./461.5./(T_expt+273.15).^2;
%E_pert=((RG+LWfb*T_pert).*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q)-E_ctrl;
%P_pert=E_pert-divFl;

% perturbation E-minus-P
E_m_P_pert =  E_m_P_expt-E_m_P_ctrl;                   % E-P in m/yr


figure
xt = [-60 -45 -30 -15 0 15 30 45 60];
xtl = ['-60'; '-45'; '-30'; '-15'; '  0'; ' 15'; ' 30'; ' 45'; ' 60'];
xt = sind(xt);
xlabel('Latitude')
ylabel('W m^{-2}')
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);
subplot(2,2,1)
yyaxis right
plot(x,-B,'linewidth',2); grid on
ylabel('W m^{-2} K^{-1}')
ylim([-11 6]);
yyaxis left
plot(x,R_frc,'linewidth',2); grid on
hold on; grid on;
% plot(x,G,'linewidth',2,'color','k','linestyle','-'); grid on
xlabel('Latitude')
ylabel('W m^{-2}')
ax = gca;
ylim([-11 6]);

ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);

title(model_flg,'qT R_{f}, \lambda')
set(gca,'fontsize',12)
%ylim([-15 15]);
legend('Forcing R_{f}','Feedback \lambda')


subplot(2,2,2)
plot(x,0.5*((T_expt-T_ctrl)+flipud(T_expt-T_ctrl)),'LineWidth',2); hold on;
plot(x,interp1(lat,0.5*(surf_t_anom_upT+flipud(surf_t_anom_upT)),phi,'linear'),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
ylim([-10 4])
title(model_flg,'qT T anomaly')
grid on;
ylabel('deg C')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

mebm_p_anom_qT = P_expt-P_ctrl;
mebm_e_anom_qT = E_expt-E_ctrl;
mebm_e_qT = E_expt;
mebm_p_qT = P_expt;

mebm_e_m_p_qT = E_m_P_expt;
mebm_e_m_p_anom_qT = E_m_P_expt-E_m_P_ctrl;






% Code for Moist/Dry Energy Balance Model
% Using old Crank Nicolson algorithm

% This code solves the energy balance model used in Roe et al. (Nat. Geosci., 2015)
% The model operates in perturbation mode (i.e., the anomaly with respect to a prescribed climatology):
% You can specify:-
% the pattern of climate forcing,Rf
% the pattern of climate feedbacks, -B
% the pattern of diffusivity, D
% whether you diffuse moist static energy, or just sensible heat

%time step in fraction of year
delt=1./50000; disp(['delt = ' num2str(delt)])
%delt=.1;
NMAX=60000; disp(['NMAX = ' num2str(NMAX)])

%set up x array (latitude).
jmx=101; 
%jmx = 101; disp(['jmx = ' num2str(jmx)])
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:);
phi = asin(x)*180/pi;

% climate parameters
A=203.3; % Size of longwave cooling constant [W/m2]
B=2.09*ones(size(x)); B=B(:); disp('B is constant for control')
% longwave cooling [W/(m2 K)]

% magnitude of diffusivity [units?]
Dmag = 0.2598; disp(['D = ' num2str(Dmag) ' W/(m2 K)'])% D = 0.2598 W/(m2 K) is the value used by TF10 
D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
 
% I think this C = rho * c * h_ml /(pi*1e7).
% this is consistent with a ~1m layer of dirt
% Note - heat capacity over LAND LAND LAND.
Cl = 0.2; disp(['Cl = ' num2str(Cl) ' J /(m2 K)'])

% Moisture parameters
relhum = 0.8;               % relative humidity
eps = 0.622;                % moisture coonstant
psfc = 9.8e4;               % (Pa)
e0 = 611.2;                 % vap. press (Pa)
a = 17.67; b = 243.5;       % sat vap constants !!T must be in temperature

L = 2.45e6;                 % latent heat of vaporization (J kg-1)
cp = 1004;                  % (J kg-1 K-1)
Re = 6.4e6;                 % [m] earth's radius
rho = 1e3;                  % [kg m-3] density

%-------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%For dry energy balance model, uncomment these values
% Dmag = 0.44; % magnitude of diffusivity [W/(m2 K)]
% D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity
% relhum = 0;  % switch off humidity
% disp('Diffusing sensible heat only')
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%-------------------------------------------------
 
% Use setupfast to create a matrix that calculates D*d/dx[ (1-x^2)d/dx] of
% whatever it operates on.
 [~,Mdiv]=setupfastM(delx,jmx,D,0.,1.0,delt);
 

% this work ok. GHR mar 3,17
load T_ctrl_FromGHR.mat;                        % control temp from Nick Siler, based on ERA
% T_ctrl = interp1(x_ghr,T_ghr,x,'pchip');  % interpolate to our grid
% T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl));            % symmetrize (optional)
% T_ctrl = interp1(lat,surf_t_ctrl,phi,'linear');

% 2-order polynomial; works for stylized climatology with ITCZ on the
% equator
% P2 = 0.5*(3*x.^2 - 1);
% T_ctrl = 13-31*P2;

q_ctrl = eps*relhum/psfc*e0*exp(a*T_ctrl./(b+T_ctrl)); q_ctrl=q_ctrl(:);% here T is in oC. q is g kg-1
theta_e_ctrl = 1/cp*(cp*(T_ctrl+273.15) + L*q_ctrl); % theta is mse divided by cp; note units of Kelvin are needed!!!

%%
%--------------------------------------------------------------
% Now do climate change experiment
disp('**Doing climate experiment**')

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% For flat forcing, uniform feedbacks, pick these values for R_frc and B
% R_frc = 4.0; % forcing in [W/m2]
% B = 2.09*ones(size(x)); B=B(:); disp('B is constant for expt.')

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% To use feedback patterns from an aquaplanet model (Feldl and Roe, GRL,
% 2013), uncomment these values
%%

B = lam_upH;
R_frc = qH;
% load NicoleNewH2O.mat
% lam_tot = 0.5*(lam_tot + fliplr(lam_tot));
% B = -interp1(sind(lat),lam_tot,x,'pchip')-0.5; % note nudge of 0.5 GHR jan21
% R_frc = interp1(sind(lat),delRf,x,'pchip');


%--------------------------------------------------------------
%set up inital T profile 
T = 0.5*(1-1*x.^2); T = T(:);
%load T_final
Tinit=T;
Tglob = mean(Tinit);

% Timestepping loop
for n=1:NMAX
   Tglob_prev = Tglob;
       
% Calculate src for this loop.
   src = R_frc/Cl; src=src(:);
   
% spec. hum, and theta_e
   q = eps*relhum/psfc*e0*exp(a*(T+T_ctrl)./(b+(T+T_ctrl))); q=q(:);% here T is in oC. q is g kg-1
   theta_e = 1/cp*(cp*((T+T_ctrl)+273.15) + L*q); % note units of Kelvin are needed!!!
   theta_e_pert = theta_e-theta_e_ctrl; %anomalous theta_e

% Calculate new T.
% Diffuse sensible heat energy (T)
% Diffusion of T anomaly
%   dT = delt/Cl*(R_frc - (B.*T) + Mdiv*T); if n == 1,disp(['Diffusing T anomaly']),end;

% Diffuse moist static energy (theta_e)
% diffusion of theta_e_anomaly
  dT = delt/Cl*(R_frc - (B.*T) + Mdiv*theta_e_pert); if n == 1, disp(['Diffusing theta_e anomaly']), end;

   T = T + dT; 
   
% Check to see if global mean temperature has converged
   Tglob=mean(T);
   Tchange = Tglob-Tglob_prev;
   if (abs(Tchange) < 1.0e-12), break; end
end

%divF_expt = -Mdiv*T;
q_expt = q;
q_pert = q_expt-q_ctrl;
divF_pert = -Mdiv*theta_e_pert;
h_pert = theta_e_pert*cp;
T_pert = T;
T_expt = T+T_ctrl;
src_pert = R_frc;

%  
%% Module to calculate terms in hyrological cycle
%define some things from the climatology first
idx = find(x>=0);
x15 = sind(15);xw = 0.15;
% wt = 0.5*(1+tanh((x(idx)-x15)/xw));     % weighting function for eddy fluxes
% wt = [flipud(wt); wt(2:end)]; %2-sided climate response
 
% some things about the control climate
Dhlf = 0.5*(D(1:end-1)+D(2:end));                           % [W m-2 K-1] calculate D on same grid as T,q, etc
h_ctrl = cp*(T_ctrl+273.15)+L*q_ctrl;                       % [J kg-1] control climate mse
F_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_ctrl,x);  % [W] control climate flux 
F_lh_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(L*q_ctrl,x);    % [W] control latent heat flux
F_hc_ctrl = (1-wt).*F_ctrl;                                 % [W] Hadley Cell Flux
% dh_gms = 1.06;                                             % [J kg-1] gross moist stability
heq_ctrl = h_ctrl(x==0);                                    % [J kg-1] moist static energy at the surface
V_ctrl = F_hc_ctrl./(heq_ctrl*dh_gms-h_ctrl);               % [kg s-1] Diagnosed mass transport in Hadley Cell
F_LH_ctrl = -L*V_ctrl.*q_ctrl + wt.*F_lh_ctrl;              % [W] latent heat (Hadley plus eddy)
divF_LH_ctrl = 1/(2*pi*Re^2)*gradient(F_LH_ctrl,x);         % [W m-2] divergence of latent heat flux
E_m_P_ctrl = divF_LH_ctrl;                   % [m yr-1]E-P

% code from Nick Siler to partition E-P;
alpha=L./461.5./(T_ctrl+273.15).^2;                         % Nick's alpha parameter
beta=cp/L./alpha./q_ctrl;                                   % beta parameter
RG=180*(1*(1-x.^2)-.4*exp(-(x./.15).^2));                   % [W m-2] idealized R-G
%RG=170*(1*(1-x.^2));%-.4*exp(-(x./.15).^2))
Ch=1.5e-3;                                                  % drag coefficient
LWfb=0;                                                     % LW feedback at surface, in W/m2/K
u=4+abs(sin(pi*x./1.5))*4;%-2.5*cos(3*asin(x));             % wind speed
rho_air=1.2;%psfc./287./(T_ctrl+273.15);                        % air density
E_ctrl=(RG.*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q_ctrl);
P_ctrl=E_ctrl-divF_LH_ctrl;


% some things about the new climate state
h_expt = h_ctrl+h_pert;                                 % [J kg-1] new total mse
F_expt = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_expt,x);      % [W] new total flux
F_lh_expt = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(L*q_expt,x);    % [W] control latent heat flux
F_hc_expt = (1-wt).*F_expt;                             % [W] Hadley Cell Flux
% dh_gms = ;                                           % [J kg-1] gross moist stability (unchanged)
heq_expt = h_expt(x==0);                                % [J kg-1] moist static energy at the surface
V_expt = F_hc_expt./(heq_expt*dh_gms-h_expt);           % [wkg s-1] Diagnosed mass transport in Hadley Cell
F_LH_expt = -L*V_expt.*q_expt + wt.*F_lh_expt;          % [W] latent heat (Hadley plus eddy)
divF_LH_expt = 1/(2*pi*Re^2)*gradient(F_LH_expt,x);     % [W m-2]
E_m_P_expt = divF_LH_expt;               % [m yr-1] E-P

alpha=L./461.5./(T_expt+273.15).^2;                         % Nick's alpha parameter
beta=cp/L./alpha./q_expt;                                   % beta parameter
E_expt=(RG.*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q_expt);
P_expt=E_expt-divF_LH_expt;


% P and E in the new climate state
%alpha=L./461.5./(T_expt+273.15).^2;
%E_pert=((RG+LWfb*T_pert).*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./L./q)-E_ctrl;
%P_pert=E_pert-divFl;

% perturbation E-minus-P
E_m_P_pert =  E_m_P_expt-E_m_P_ctrl;                   % E-P in m/yr








xt = [-60 -45 -30 -15 0 15 30 45 60];
xtl = ['-60'; '-45'; '-30'; '-15'; '  0'; ' 15'; ' 30'; ' 45'; ' 60'];
xt = sind(xt);
xlabel('Latitude')
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);
subplot(2,2,3)
yyaxis right
plot(x,-B,'linewidth',2); grid on
ylabel('W m^{-2} K^{-1}')
ylim([-11 6]);
yyaxis left
plot(x,R_frc,'linewidth',2); grid on
hold on; grid on;
% plot(x,G,'linewidth',2,'color','k','linestyle','-'); grid on
xlabel('Latitude')
ylabel('W m^{-2}')
ax = gca;
ylim([-11 6]);

ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);

title(model_flg,'qH R_{f}, \lambda')
set(gca,'fontsize',12)
%ylim([-15 15]);
legend('Forcing R_{f}','Feedback \lambda')



subplot(2,2,4)
plot(x,0.5*((T_expt-T_ctrl)+flipud(T_expt-T_ctrl)),'LineWidth',2); hold on;
plot(x,interp1(lat,0.5*((surf_t_anom_upH)+flipud(surf_t_anom_upH)),phi,'linear'),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
ylim([-10 4])
title(model_flg,'qH T anomaly')
grid on;
ylabel('deg C')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

figure
subplot(2,2,1)
plot(x,0.5*(mebm_e_m_p_qT+flipud(mebm_e_m_p_qT)),'LineWidth',2); hold on;
plot(x,interp1(lat,0.5*(e_m_p_upT+flipud(e_m_p_upT)),phi,'linear'),'LineWidth',2); hold on;grid on
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT E - P')
ylabel('W/m2')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,2)
plot(x,0.5*(E_m_P_expt+flipud(E_m_P_expt)),'LineWidth',2); hold on;
plot(x,interp1(lat,0.5*(e_m_p_upH+flipud(e_m_p_upH)),phi,'linear'),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH E - P')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,3)
plot(x,0.5*(mebm_e_m_p_anom_qT+flipud(mebm_e_m_p_anom_qT)),'LineWidth',2); hold on;
aa0 = interp1(lat,e_m_p_anom_upT,phi,'linear');
plot(x,0.5*(aa0+flipud(aa0)),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT E-P anomaly')
grid on;
ylabel('W/m2')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,4)
plot(x,0.5*((E_m_P_expt-E_m_P_ctrl)+flipud(E_m_P_expt-E_m_P_ctrl)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,e_m_p_anom_upH,phi,'linear'))+flipud(interp1(lat,e_m_p_anom_upH,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH E-P anomaly')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)


figure
subplot(2,2,1)
plot(x,0.5*((mebm_p_qT)+flipud(mebm_p_qT)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,precip_lat_upT,phi,'linear'))+flipud(interp1(lat,precip_lat_upT,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT Precip')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,2)
plot(x,0.5*((P_expt)+flipud(P_expt)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,precip_lat_upH,phi,'linear'))+flipud(interp1(lat,precip_lat_upH,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH Precip')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,3)
plot(x,0.5*((mebm_p_anom_qT)+flipud(mebm_p_anom_qT)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,p_anom_upT,phi,'linear'))+flipud(interp1(lat,p_anom_upT,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT Precip anomaly')
grid on;
ylabel('W/m2')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,4)
plot(x,0.5*((P_expt-P_ctrl)+flipud(P_expt-P_ctrl)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,p_anom_upH,phi,'linear'))+flipud(interp1(lat,p_anom_upH,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH Precip anomaly')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

figure
subplot(2,2,1)
plot(x,0.5*((mebm_e_qT)+flipud(mebm_e_qT)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,evap_lat_upT,phi,'linear'))+flipud(interp1(lat,evap_lat_upT,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT Evap')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,2)
plot(x,0.5*((E_expt)+flipud(E_expt)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,evap_lat_upH,phi,'linear'))+flipud(interp1(lat,evap_lat_upH,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH Evap')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,3)
plot(x,0.5*((mebm_e_anom_qT)+flipud(mebm_e_anom_qT)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,e_anom_upT,phi,'linear'))+flipud(interp1(lat,e_anom_upT,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qT Evap anomaly')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

subplot(2,2,4)
plot(x,0.5*((E_expt-E_ctrl)+flipud(E_expt-E_ctrl)),'LineWidth',2); hold on;
plot(x,0.5*((interp1(lat,e_anom_upH,phi,'linear'))+flipud(interp1(lat,e_anom_upH,phi,'linear'))),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title(model_flg,'qH Evap anomaly')
grid on;
ylabel('W/m2')
xlabel('latitude')
legend('MEBM',model_flg)
set(gca,'fontsize',12)

figure
plot(x,wt,'linewidth',2)
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
title('(gms =',string(dh_gms))
xlabel('latitude')
set(gca,'fontsize',12)
