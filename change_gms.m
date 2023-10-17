clear all; close all; clc;
% test effect of changing gms
gms = 0.6;

% Physical constants
psfc = 1.013e5;  % [Pa] surface pressure
g = 9.81;        % [m s^-2] er....
Re = 6.37e6;     % [m] Earth's radius 
cp = 1004;       % [K kg^-1 K^-1] specific heat constant pressure
Lv = 2.45e6;     % [J kg-1] latent heat of vaporization (J kg-1)

% Thermodynamics and moisture parameters
relhum = 0.8;   % relative humidity
eps = 0.622;    % moisture coonstant (Rd/Rv I think)
e0 = 611.2;     % vap. press (Pa)
a = 17.67; b = 243.5;   % sat vap constants !!T must be in temperature
Rv = 461.5;     % gas constant for water vapor
rho_w = 1e3;    % [kg m-3] water density

% Time array
ts = 0;     % [yr] Initial time (integer)
tf = 1000;   % [yr] Final time (integer)

%time step in fraction of year
% make denominator an integer multiple of 360
% if code blows up, try shortening the time step (this is explicit scheme)
delt=1./(360*6); disp(['delt = ' num2str(delt)]) % [yr^-1] trial and error 140 time steps a day.

% number of time steps
nts = (tf-ts)/delt + 1;

% time indicators
t = ts;             % [yr] time years (could be past climate for Milankovitch)
yr = 0;             % [yr] model years, gets updated immediately time loops starts
newyear_flg = 1;    % is it a new year? set to yes (=1) initially

% other flags
two_layer_flg = 0;      % two layer ocean model? (1=yes, 0=no)
noise_flg = 1;          % noise flag (can be specified in many ways)
moistEBM_flg = 1;       % flag for moist or dry diffusion (1 = moist)

% if dry, set Lv = 0;
if moistEBM_flg == 0; 
    Lv = 0; rel_hum = 0;                    % switch off latent-heat dependence
    disp('****Diffusing sensible heat only - dry EBM *****')
else
    disp('****Diffusing moist enthalpy - dry EBM *****')
end

%set up x array (latitude).
jmx=101; 
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';x = x(:); % sine latitude
phi = asin(x)*180/pi;                        % latitude

%% diffusivity; units conversion is always tricky
if moistEBM_flg==1
    D_HF10 = 1.06e6; % [m^2 s^-1] Hwang and Frierson (GRL, 2010)
else
    D_HF10 = 1.70e6;                % diffusivity for sensible heat.
end

% convert diffusivity units for this numerical code
Dmag = D_HF10*psfc*cp/(g*Re^2); % [W m^-2 K^-1] conversion to regular units here. For HF10, D~0.28; For dry EBM see Armouretal19, but should be ~0.44
disp(['D = ' num2str(Dmag) ' W/(m2 K)'])% D = 0.2598 W/(m2 K) is the value used by TF10 

D=Dmag*ones(jmx+1,1); D=D(:); % diffusivity array (at each cell edge)
D_mid = 0.5*(D(1:end-1)+D(2:end)); % diffusivity stipulated at cell midpoints.

%--------------------------------------------------------------------------------------
%% Feedbacks, forcing, and heat uptake - various choices
%--------------------------------------------------------------------------------------
% Load in CMIP5 feedbacks, forcing, and heat uptake here if desired.....
% CMIP5 ensemble-mean feedback and forcing values from 4xCO2 simulations (taken at year 100)
% Note these are different from Bonan et al., and Siler et al., so figrue
% that out at some point.
  load CMIP5_Rf_G_lambda.mat %feedback, forcing, and heat uptake for 11 models
  B = -interp1(CMIP5_lat,mean(CMIP5_lambda,2),phi,'linear');
  R_frc = interp1(CMIP5_lat,mean((CMIP5_Rf),2),phi,'linear'); %CO2 forcing Rf plus ocean heat uptake G
  G = interp1(CMIP5_lat,mean((CMIP5_G),2),phi,'linear'); %CO2 forcing Rf plus ocean heat uptake G

R_frc = mean(R_frc)*ones(size(R_frc)); %take mean value of forcing (optional)
B = mean(B)*ones(size(B)); %take mean value of feedbacks (optional)
G = 0*mean(G)*ones(size(G)); %take mean value of ht uptake (optional)

%--------------------------------------------------------------------------------------
%% Control climate - various choices
%--------------------------------------------------------------------------------------
% %% load in the climatological temperature from NCEP for the control climate
% load NCEP_ClimoObsTsfc.txt
% lat = NCEP_ClimoObsTsfc(:,1);
% T_ctrl = NCEP_ClimoObsTsfc(:,2);  % import observed T; units K
% T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl)); % average N & S hemispheres.
% T_ctrl = interp1(sind(lat),T_ctrl,x,'pchip')-273.15;

 %load T_ctrl_FromNick.mat;                        % control temp from Nick Siler, based on ERA
 %T_ctrl = interp1(x_nick,T_ctrl_nick,x,'pchip');  % interpolate to our grid
 %T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl));            % symmetrize (optional)
 %T_ctrl = T_ctrl + 2.*exp(-(x/0.15).^4);          % add a bump [n.b. still a little structure here]

% this work ok. GHR mar 3,17
load T_ctrl_FromGHR.mat;                        % control temp from GHR, 
T_ctrl = 1.1*interp1(x_ghr,T_ghr,x,'pchip');        % interpolate to our grid
T_ctrl = 0.5*(T_ctrl+flipud(T_ctrl));           % symmetrize (optional)

% control climate values of q and theta_e
q_ctrl = eps*relhum/psfc*e0*exp(a*T_ctrl./(b+T_ctrl)); q_ctrl=q_ctrl(:);% here T is in oC. q is g kg-1
theta_e_ctrl = 1/cp*(cp*(T_ctrl+273.15) + Lv*q_ctrl); % note units of Kelvin are needed!!! 

%--------------------------------------------------------------------------
%% heat capacity, mixed layer depth
%--------------------------------------------------------------------------
% two-layer coupling parameter (2 layer = 1, 1 layer = 0)
if two_layer_flg==1
    gamma = 0.67*ones(size(x)); % [W m^-2 K-1] coupling b/t surface and deep ocean, based on CMIP 5 fit from Armour (NCC 2017)
    disp('!!!!!!!! Using two layer ocean WATCH FOR DRIFT!!!!!!')
else
    gamma = 0*ones(size(x));    % no coupling b/t surface and deep ocean
end
disp(['two-layer coupling gamma = ' num2str(mean(gamma)) ' W m-2 K-1'] ); 

mix_depth = 35; %[m] mixed-layer depth
h_ml = mix_depth*ones(size(x)); h_ml = h_ml(:); % [m] depth of mixed layer (assuming water), can be specified with latitude.
rho_w = 1e3; % [kg m^-3]
cw = 4200;  % [J kg^-1 K^-1] heat capacity of water

% note: units of heat capacity allow units of years for time
C_L1 = rho_w * cw * h_ml /(pi*1e7); % heat capacity of layer 1 (upper layer)
disp(['h_mix_lyr = ' num2str(mix_depth) 'm; C_L1 = ' num2str(mean(C_L1)) ' J /(m2 K s yr^-1)'])

% deep layer for two-layer ocean model
h_d = 800; % [m] - deep-layer depth based roughly on CMIP5 fit from Armour(NCC, 2017)
h_ml = mix_depth*ones(size(x)); h_ml = h_ml(:); % [m] depth of deep-ocean layer, can be specified with latitude.
C_L2 = rho_w * cw * h_d /(pi*1e7); % heat capacity of layer 2 (upper layer)
disp(['h_deep = ' num2str(h_d) 'm; C_L2 = ' num2str(mean(C_L2)) ' J /(m2 K s yr^-1)'])

%set up inital T perturbation; set to zero, since this is perturbation
T_init = zeros(size(x)); T_init = T_init(:);
T = T_init;
Td = T; % give deep ocean temperature the initial temperature. Revise as necessary

% Use setupfast to create a matrix that calculates D*d/dx[ (1-x^2)d/dx] of
% whatever it operates on.
[~,Mdiv]=setupfastM(delx,jmx,D,0.,1.0,delt);

 %Global mean temperature
Tglob=mean(T);

%--------------------------------------------------------------------------
% Timestepping loop
%--------------------------------------------------------------------------
for n=1:nts-1;          % the -1 stops the loop before getting to yr_end +1
   Tglob_prev = Tglob;
   
% time indexing, output every month
    t(n) = ts + (n-1)*delt; % [yrs] time in years

% need to update days, months, and years
% is it a new year?
    tmp = floor(t(n))+1; 
    if ((tmp~=yr)||(n==1));
        yr = yr+1;
        disp(['happy new year! yr = ' num2str(yr)])
        % add interannual noise
        if noise_flg==1 % Need to choose noise forcing....
        %    noise = 15.0*randn(jmx,1);  % [W m-2] dfferent random noise at each latitude
        %    noise = 2.0*randn(1,1)*ones(jmx,1); % [W m-2] same random noise at all latitudes 
            noise = 0*15.0*randn(jmx,1) +2.0*randn(1,1)*ones(jmx,1); 
        else
            noise = zeros(jmx,1);
        end
        noise = noise(:);
        newyear_flg = 1;
    end
     
% spec. hum, and theta_e (n.b. If in dry-EBM mode, then Lv,q = 0, so theta_e = T)
% need to calculate theta_e total and then subtract the climatology to calculate the perturbation
   q = eps*relhum/psfc*e0*exp(a*(T+T_ctrl)./(b+(T+T_ctrl))); q=q(:);% here T is in oC. q is g kg-1
   theta_e = 1/cp*(cp*((T+T_ctrl)+273.15) + Lv*q); % note units of Kelvin are needed!!! 
   theta_e_pert = theta_e-theta_e_ctrl; % perturbation theta_e on which the diffusion acts
   
% Calculate new T.
% Diffuse moist static energy (theta_e_pert)
% come back to worry about minus sign for Mdiv term, but I tihnk it is ok....
   dT = delt./C_L1.* (R_frc + G - (B.*T) + Mdiv*theta_e_pert + gamma.*(Td-T)+noise);
   dTd = delt./C_L2.*gamma.*(T-Td);  % deep temperature update
   
% update temperature
    T = T + dT; 
    Td = Td + dTd;
   
% Check to see if global mean temperature has converged
%   Tglob=mean(T);
%   Tchange = Tglob-Tglob_prev;
%   if (abs(Tchange) < 1.0e-5), disp(['**Converged n = ' num2str(n)]); break; end

%---------------------------
%% if new year, output data 
%---------------------------
    if (newyear_flg == 1)
       time_out(yr) = yr; 
       Tout(:,yr) = T; 
       Tdout(:,yr) = Td;   % can be junked if using ony 1layer
       newyear_flg = 0;        % set new day flag to zero
       nout(:,yr) = noise; % [W m-2] noise output
       stor_out(:,yr) = C_L1.*dT/delt + C_L2.*dTd/delt; % [W m-2] storage output; 
%      note that this is the instantaneous storage rate not annual average 
    end


%% end loop over time
end
%%

% %% Hydrocycle and graph output
% % What is below just graphs the final state of the model
% % Monthly output from the model is stored in 
% % Tout (surface layer), and Tdout (deep-ocean layer)
% % Everything else (probably) can be diagnosed or added to the output
% %% climatological fields
% needed output fields
% Climatology, T,Td,E-P, Ftot, F_LH, P,E
% Perturbation for every latitude and year: Td, E-P, F_tot, F_LH, P, E

%------------------------------------------------------
% Output including hydrology
% Everything can be diagnosed from the temperature, 
% hydrology follows Siler et al., 2018
%------------------------------------------------------
% some things about the control climate
% !!! check weighting fcuntion is same as Nicks
%wt = 1-gaussmf(x,[sind(15) 0]); %define Gaussian Hadley Cell width
%wt = 1-normpdf(x,0,sind(15)); %define Gaussian Hadley Cell width %!!!!check this is the sams as Nick
% check this is the same as in Nick's paper - I did a Matlab guess, GHR 11Jan21!!!!
wt = normpdf(x,0,sind(15))/max(normpdf(x,0,sind(15))); %Gaussian-shaped weighting function, normalized to 1.
wt = 1-wt;
Dhlf = 0.5*(D(1:end-1)+D(2:end));                           % [W m-2 K-1] calculate D on same grid as T,q, etc
h_ctrl = cp*(T_ctrl+273.15)+Lv*q_ctrl;                       % [J kg-1] control climate enthalpy
F_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_ctrl,x);  % [W] control climate flux 
F_lh_ctrl = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(Lv*q_ctrl,x);    % [W] control latent heat flux
F_hc_ctrl = (1-wt).*F_ctrl;                                 % [W] Hadley Cell Flux
heq_ctrl = h_ctrl(x==0);                                    % [J kg-1] moist static energy at the surface
V_ctrl = F_hc_ctrl./(heq_ctrl*gms-h_ctrl);                 % [kg s-1] Diagnosed mass transport in Hadley Cell (Nick's way)

F_LH_ctrl = -Lv*V_ctrl.*q_ctrl + wt.*F_lh_ctrl;              % [W] latent heat (Hadley plus eddy)
F_LH_eddy_ctrl = wt.*F_lh_ctrl;                         % eddy latent heat fluxes including weighting function.
divF_LH_ctrl = 1/(2*pi*Re^2)*gradient(F_LH_ctrl,x);         % [W m-2] divergence of latent heat flux
E_m_P_ctrl = divF_LH_ctrl;                                  % [W m-2]E-P

% code from Nick Siler to partition E-P;
alpha=Lv./(Rv *(T_ctrl+273.15).^2);                      % Nick's alpha parameter
beta=cp/Lv./alpha./q_ctrl;                                   % beta parameter
RG=180*(1*(1-x.^2)-.4*exp(-(x./.15).^2));                   % [W m-2] idealized R-G
Ch=1.5e-3;               % drag coefficient                                   % drag coefficient
LWfb=0;                                                     % LW feedback at surface, in W/m2/K
u=4+abs(sin(pi*x./1.5))*4;%-2.5*cos(3*asin(x));             % wind speed
rho_air=1.2;%psfc./287./(T_ctrl+273.15);                        % air density
E_ctrl=(RG.*alpha+rho_air.*cp.*(1-relhum).*Ch.*u)./(alpha+cp./Lv./q_ctrl); % [W m-2]
P_ctrl=E_ctrl-divF_LH_ctrl;                                                % [W m-2]

%E_m_P_ctrl = E_m_P_ctrl/(Lv*rho_w)*pi*1e7;  % [m yr-1]. Convert from W m-2.
%P_ctrl=P_ctrl/(Lv*rho_w)*pi*1e7;            % [m yr-1]. Convert from W m-2.
%E_ctrl=E_ctrl/(Lv*rho_w)*pi*1e7;            % [m yr-1]  Convert from W m-2.

%------------------------------------------------
% hydrology of the output
% loop over every year of model output
% output 
% E_expt, P_expt, E_m_P_expt, T_out, F_expt, 
%------------------------------------------------
for i = 1:tf            % loop over years

    % % some things about the new climate state
    % note T_ctrl is in oC and  
    Ttry = Tout(:,i);
    q_expt = eps*relhum/psfc*e0*exp(a*(T_ctrl+Ttry)./(b+(T_ctrl+Ttry))); q_expt=q_expt(:);% here T is in oC. q is g kg-1
    h_expt(:,i) = cp*(Ttry+T_ctrl+273.15)+Lv*q_expt;                    % [J kg-1] control climate enthalpy
    F_expt(:,i) = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(h_expt(:,i),x);      % [W] new total flux
    divF_expt(:,i) = 1/(2*pi*Re^2)*gradient(F_expt(:,i),x);
    F_lh_expt = -2*pi*Re^2/cp*Dhlf.*(1-x.^2).*gradient(Lv*q_expt,x);    % [W] control latent heat flux
    F_hc_expt = (1-wt).*F_expt(:,i);                             % [W] Hadley Cell Flux
    heq_expt = h_expt(x==0,i);% [J kg-1] moist static energy at the surface
    V_expt = F_hc_expt./(heq_expt*gms-h_expt(:,i));           % [wkg s-1] Diagnosed mass transport in Hadley Cell

    % 
    F_LH_expt = -Lv*V_expt.*q_expt + wt.*F_lh_expt;          % [W] latent heat (Hadley plus eddy)
    F_LH_eddy_expt = wt.*F_lh_expt;                         % eddy latent heat fluxes including weighting function.
    divF_LH_expt = 1/(2*pi*Re^2)*gradient(F_LH_expt,x);     % [W m-2]
    E_m_P_expt(:,i) = divF_LH_expt;                              % [W m-2] E-P

    % note have chosen to keep alph and beta at climatology; could change if desired
    alpha=Lv./(Rv*(T_ctrl+273.15).^2);                         % Nick's alpha parameter 
    beta=cp./(Lv*alpha.*q_ctrl);                               % beta parameter

    % Equation 16 from Siler et al. (JClim,2018) 
    E_pert(:,i) = (E_ctrl.*Ttry.*beta.*(alpha-2./(T_ctrl+273.15))-G)./(1+beta); % perturbation from Nick's equation
    E_expt(:,i) = E_ctrl+E_pert(:,i);     % [W m-2] - evporation of the experiment

    P_expt(:,i) = E_expt(:,i) -divF_LH_expt; % [W m-2] - precipitation of the experiment

% end hydro loop over years
end

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
% ylim([-11 6]);
yyaxis left
plot(x,R_frc,'linewidth',2); grid on
hold on
plot(x,G,'linewidth',2,'color','k','linestyle','-'); grid on
xlabel('Latitude')
ylabel('W m^{-2}')
ax = gca;
% ylim([-11 6]);

ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl);

title('qT R_{f}, G, \lambda')
set(gca,'fontsize',12)
%ylim([-15 15]);
legend('Forcing R_{f}','G','Feedback \lambda')

subplot(2,2,2)
plot(x,Tout(:,end),'LineWidth',2); hold on;
%plot(x,interp1(lat,surf_t_anom_upT,phi,'linear'),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
% ylim([-10 4])
title('qT T anomaly')
grid on;
ylabel('deg C')
%legend('MEBM',model_flg)
set(gca,'fontsize',12)

mebm_p_anom_qT = P_expt(:,end)-P_expt(:,1);
mebm_e_anom_qT = E_expt(:,end)-E_expt(:,1);

mebm_e_m_p_qT = E_m_P_expt(:,end);
mebm_e_m_p_anom_qT = E_m_P_expt(:,end)-E_m_P_ctrl;


subplot(2,2,3)
%yyaxis("left")
plot(x,mebm_e_m_p_qT,'LineWidth',2); hold on;grid on;
%yyaxis("right")
%%plot(x,interp1(lat,e_m_p_upT,phi,'linear'),'LineWidth',2); hold on;grid on
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title('qT E - P')
ylabel('W/m2')
%legend('MEBM',model_flg)
set(gca,'fontsize',12)



subplot(2,2,4)
%yyaxis("left")
plot(x,mebm_e_m_p_anom_qT,'LineWidth',2); hold on;
%yyaxis("right")
%plot(x,interp1(lat,e_m_p_anom_upT,phi,'linear'),'LineWidth',2); hold on;
ax = gca;
set(gca,'XTick',xt); set(gca,'XTicklabel',xtl)
%ylim([-10 5])
title('qT E-P anomaly')
grid on;
ylabel('W/m2')
%legend('MEBM',model_flg)
set(gca,'fontsize',12)





