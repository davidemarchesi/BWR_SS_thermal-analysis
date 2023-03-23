%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Heat Transfer Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
BWR_parameters;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Pt-1 : Single Phase Heat Transfer Analysis %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Is made a discretization all along the length 'z' of the rod.
% At each length 'z' is made a single-phase liquid heat transfer
% calculation.
% The iteration stops at the point in which the water will
% start to nucleate, the so called ONB point. 
% Once that this point is reached, will have to be considered another heat 
% transfer description, the subcooled boiling, that will be made below.
%

% Declaration of the variables that are used
delta_Ts_sat = 0 ;
delta_T_ONB = 0 ;
Tf_ONB = 1 ;
Tf = 0 ;
Ti = T_in + 273.15 ; % Temperature in (K)
Cp_l = XSteam('Cp_pT',P_xsteam,T_in) ; % Specific Heat of liquid (kJ/kg K)

z_ONB = 0 ;
z1 = -1.8 ; % lower bound of the rod (m)
z2 = z1 + 0.001 ; % discretization in the order of the 'mm'
z_mean = (z1 + z2)./2 ;

k_fl = 0.567 ; % thermal conductivity of the fluid (W/m*K)
rho_g = XSteam('rhoV_p',P_xsteam) ; % saturated vapour density (kg/m^3)
rho_l = XSteam('rhoL_p',P_xsteam) ; % saturated liquid density (kg/m^3)
Tsat = XSteam('Tsat_p',71.4) ; % saturation temperature (deg C)
sigma = XSteam('st_T',Tsat) ; % surface tension
h_in = XSteam('h_pT',P_xsteam,T_in) ; % inlet enthalpy (kJ/kg)

% variables needed for the results and the iterations
it1 = 0 ;
Temp_bulk = [] ;
Temp_lo_wall = [] ;
alpha = [] ;
z_T = [] ;
delta_Ts_sat_vec = [] ;
x_eq_vec = [] ;
q_fl_vec_tot = [] ;

while  Tf <= Tf_ONB 

    % Determination of Tf through an energy balance on a bin 'delta_z'
    Tcpl = Ti - 273.15 ; % temperature used (in °C) to calculate the Cp of water
    Cp_l = XSteam('Cp_pT',P_xsteam,Tcpl) ;
    Q_z = ax_linear_heat_rod(z_mean).*abs(z2 - z1) ; % heat exchanged in the bin (z1,z2)
    Tf = ( Q_z./(m_ch.*Cp_l) ) + Ti ; % outlet bin temperature (K)
    
    % Determination of the wall T and of the heat transfer coefficient
    % using the heat transfer equation and the Dittus-Boelter correlation.
    Pr_lo = (Cp_l.*1000.*mu_f)./ k_fl ; % determination of the fluid Prandtl number
    Re_lo = (G_ch.*D_eq)./ mu_f ; % determination of the fluid Reynolds number

    if Re_lo < 30000 % the condition for the applicability of the Dittus Boelter would be Re >> 1e4
        error('\nThe Reynolds number do not respect the Dittus Boelter applicability ',...
                       ' due to the fact that is not >>10e4');
        break
    end

    Nu_lo = 0.023.*((Re_lo).^0.8).*((Pr_lo).^0.4) ; % determination of the Nusselt number
    alpha_lo = Nu_lo.* (k_fl./D_eq) ; % determination of the heat transfer coefficient (W/m^2 K)
    Twall = (((ax_linear_heat_rod(z_mean).*1000)./(pi.*D_r_out.*0.001))./alpha_lo) + Tf ; % wall temperature (K)
    qfl_i = (ax_linear_heat_rod(z_mean))./(pi.*D_r_out.*0.001) ; % average heat flux in the bin delta_z (kW/m^2)
    q_fl_vec_tot = [q_fl_vec_tot,qfl_i] ;
    delta_Ts_sat = Twall - (Tsat + 273.15) ;

    % After the determination of the Twall in the bin , has ot be checked
    % that we are continuing to stay above the nucleation limit.
    % In order to do this we have to check the ONB point.
    % Is used the Frost-Dzakovic model, for a more correct evaluation of
    % the Davis-Anderson formula for the ONB.
    delta_T_ONB_DH = sqrt((8.*sigma.*(Tsat+273.15).*((ax_linear_heat_rod(z_mean).*1000)./(pi.*D_r_out.*0.001)))./(k_fl.*h_lg.*1000.*rho_g)) ;
    delta_T_ONB = delta_T_ONB_DH.*Pr_lo ;
    Tf_ONB = (Tsat+273.15) + delta_T_ONB - (((ax_linear_heat_rod(z_mean).*1000)./(pi.*D_r_out.*0.001))./alpha_lo) ;
    
    % calculation of the vapour equilibrium quality in the bin
    Tf_degc = Tf -273.15 ;
    h_x = XSteam('h_pT',P_xsteam,Tf_degc);
    h_satl = XSteam('hL_p',P_xsteam);
    x_eq = (h_x-h_satl)./h_lg ;

    it1 = it1 + 1 ;
    
    z_T = [z_T, z2] ;

    x_eq_vec = [x_eq_vec , x_eq] ;

    Temp_bulk = [Temp_bulk , Tf] ;
    
    Temp_lo_wall = [Temp_lo_wall , Twall] ;

    alpha = [alpha , alpha_lo] ;

    delta_Ts_sat_vec = [delta_Ts_sat_vec , delta_Ts_sat] ;

    % change of variables to continue to iterate the cicle
    z1 = z2 ;
    z2 = z2 + 0.001 ;
    z_mean = (z1 + z2)./2 ;
    Ti = Tf ;

end

figure(2)
hold on
title('Temperature profile')
xlabel('Axial Position (m)')
ylabel('Temperature (K)')
plot(z_T,Temp_bulk,LineWidth=1.5,LineStyle="-")
hold on
plot(z_T,Temp_lo_wall,LineWidth=1.5,LineStyle="--")
legend('Bulk Temperature T_b','Wall Temperature T_w',Location='best')
hold off

figure(3)
hold on
title('Heat transfer coefficient')
xlabel('Axial Position (m)')
ylabel('\alpha (W/(m^2*K))')
plot(z_T,alpha,LineWidth=2)
hold off

figure(4)
hold on
title('Equlibrium quality x_{eq}')
xlabel('Axial Position (m)')
ylabel('x_{eq}')
plot(z_T,x_eq_vec,LineWidth=2)
hold off

z_ONB = z1-0.001;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Pt-2 : Subcooled Boiling Analysis %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Subcooled Boiling Analysis could be faced with different types of
% approaches.
% As above there's a discretization along the axial axes z and the
% equations are solved at each point up to the saturation temperature.
% A) Non-linear superposition method, with :
%       1 - The Jens-Lottes correlation ;
%       2 - The Thom correlation ;
%       3 - The Rohsenow correlation ;
%       4 - The Mostinski-Palens' corrected correlation
%
% B) Liu-Winterton subcooled correlation

% declaration of some parameters used in the below correlations

Csf = 0.0132 ; % Rohsenow correlation coefficient taken for mechanically polished stainless steel
Cp_lsat = XSteam('CpL_p',P_xsteam) ;
Pr_lsat = (Cp_lsat.*1000.*mu_f)./ k_fl ; % Liquid Prandtl number in Sat-cond.
p_crit = 22.064.*1000 ; % (kPa) critical pressure of water
p_crit_bar = 220.64 ;
p_rid = (P_nom./(p_crit./1000)) ; % reduced pressure
M = 18 ; % molar mass of water
Fp = 1.8.*(p_rid.^0.17) ; % Palen's correction factor for Mostinski formula

itJL_vec = [] ;
itT_vec = [] ;
itR_vec = [] ;
itM_vec = [] ;

alpha_JL = [] ;
alpha_T = [] ;
alpha_R = [] ;
alpha_LW = [] ;

Twall_JL = [] ;
Twall_T = [] ;
Twall_R = [] ;
Twall_LW = [] ;

qflux_vec = [];
toll = 10e-6 ;
nmax = 10e6 ;
a = Twall-1.5 ;
b = Twall+25.5 ;

while Tf < (Tsat + 273.15)

    % Calculation of the Tf with the energy balance
    Tcpl = Ti - 273.15 ;
    Cp_l = XSteam('Cp_pT',P_xsteam,Tcpl) ;
    Q_z = ax_linear_heat_rod(z_mean).*abs(z2 - z1) ; % (kW)
    q_flux_z = (ax_linear_heat_rod(z_mean).*1000)./(pi.*D_r_out.*0.001) ; % (W/m^2 K)
    qfl_i = (ax_linear_heat_rod(z_mean))./(pi.*D_r_out.*0.001) ;
    q_fl_vec_tot = [q_fl_vec_tot,qfl_i] ;

    Tf = ( Q_z./(m_ch.*Cp_l) ) + Ti ; % (K)
    Tf_degc = Tf -273.15 ;
    h_x = XSteam('h_pT',P_xsteam,Tf_degc);
    x_eq = (h_x-h_satl)./h_lg ;

    % Determination of the convective heat transfer coefficient
    Pr_lo = (Cp_l.*1000.*mu_f)./ k_fl ; % determination of the fluid Prandtl number
    Re_lo = (G_ch.*D_eq)./ mu_f ; % determination of the fluid Reynolds number

    if Re_lo < 30000 % the condition for the applicability of the Dittus Boelter would be Re >> 1e4
        error('\nThe Reynolds number do not respect the Dittus Boelter applicability ',...
                       ' due to the fact that is not >>10e4');
        break
    end

    Nu_lo = 0.023.*((Re_lo).^0.8).*((Pr_lo).^0.4) ; % determination of the Nusselt number
    alpha_lo = Nu_lo.* (k_fl./D_eq) ; % (W/m^2 K)
    
    % A) Non-linear superposition method
    %
    % 1 - Jens-Lottes
    K_JL = 1000000.*exp(4.*P_nom./6.2)./(25.^4) ;
    itJL = 0 ;
    
    qscb_JL =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + (((K_JL.*((Tw_tent-Tsat-273.15).^4)) - (K_JL.*(delta_T_ONB.^4))).^2)) ;
    f_JL =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + (((K_JL.*((Tw_tent-Tsat-273.15).^4)) - (K_JL.*(delta_T_ONB.^4))).^2))-q_flux_z ;

    [Tit_JL_vec,xdif1,fx1,itJL]=bisection_method(a,b,nmax,toll,f_JL) ;
    
    itJL_vec = [itJL_vec , itJL] ;
    
    Tit_JL = Tit_JL_vec(end) ;

    alpha_it_JL = qscb_JL(Tit_JL)./(Tit_JL-Tf) ;
    alpha_JL = [alpha_JL , alpha_it_JL] ;
    Twall_JL = [Twall_JL , Tit_JL] ;
    
    % 2 - Thom
    K_T = 1000000.*exp(2.*P_nom./8.7)./(22.7.^2) ;
    itT = 0 ;

    qscb_T =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + ((K_T.*((Tw_tent-Tsat-273.15).^2)) - (K_T.*(delta_T_ONB.^2))).^2) ;
    f_T =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + ((K_T.*((Tw_tent-Tsat-273.15).^2)) - (K_T.*(delta_T_ONB.^2))).^2)-q_flux_z ;

    [Tit_T_vec,xdif2,fx2,itT]=bisection_method(a,b,nmax,toll,f_T) ;
    
    itT_vec = [itT_vec , itT] ;
    
    Tit_T = Tit_T_vec(end) ;
    
    alpha_it_T = qscb_T(Tit_T)./(Tit_T-Tf) ;
    alpha_T = [alpha_T , alpha_it_T] ;
    Twall_T = [Twall_T , Tit_T] ;

    % 3 - Rohsenow
    K_R = mu_f.*h_lg.*1000.*((9.81.*(rho_l-rho_g)./sigma).^0.5).*((Cp_lsat./(h_lg.*Csf.*Pr_lsat)).^3) ;
    itR = 0 ;

    qscb_R =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + ((K_R.*((Tw_tent-Tsat-273.15).^3)) - (K_R.*(delta_T_ONB.^3))).^2) ;
    f_R =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + ((K_R.*((Tw_tent-Tsat-273.15).^3)) - (K_R.*(delta_T_ONB.^3))).^2)-q_flux_z ;
    
    [Tit_R_vec,xdif3,fx3,itR]=bisection_method(a,b,nmax,toll,f_R) ;
    
    itR_vec = [itR_vec , itR] ;
    
    Tit_R = Tit_R_vec(end) ;
    
    alpha_it_R = qscb_R(Tit_R)./(Tit_R-Tf) ;
    alpha_R = [alpha_R , alpha_it_R] ;
    Twall_R = [Twall_R , Tit_R] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%     % 4 - Mostinski (with Palen's correction)
% 
%     % ( Different possible form for the Palen correction factor :
%     % Fp = (2.1.*(p_rid.^0.27))+((((1-(p_rid.^2)).^(-1))+9).*(p_rid.^2)) 
%     % K_M = (0.00417.*(p_crit.^0.69).*Fp).^(1/0.3) ;
%     % )
% 
%     K_M = (0.1011.*(p_crit_bar.^0.69).*Fp).^(1/0.3);
%     itM = 0 ;
% 
%     qscb_M =@(Tw_tent) sqrt(((alpha_lo.^2).*((Tw_tent-Tf).^2)) + (((K_M.*((Tw_tent-Tsat-273.15).^3)) - (K_M.*(delta_T_ONB.^3))).^2)) ;
%     f_M =@(Tw_tent) (((alpha_lo.^2).*((Tw_tent-Tf).^2))+(((K_M.*((Tw_tent-Tsat-273.15).^3))-(K_M.*(delta_T_ONB.^3))).^2))-(q_flux_z.^2) ;
%     df_M = @(Tw_tent) ((2.*(alpha_lo.^2).*(Tw_tent-Tf))+(6.*(K_M.^2).*((Tw_tent-Tsat-273.15).^5))-(2.*(K_M.^2).*(delta_T_ONB.^3).*3.*(Tw_tent-Tsat-273.15))) ;
% 
% 
%     [Tit_M_vec,itM]=newton_method(595,nmax,toll,f_M,df_M) ;
%     
%     itM_vec = [itM_vec , itM] ;
%     
%     Tit_M = Tit_M_vec(end) ;
%     
%     alpha_it_M = qscb_M(Tit_M)./(Tit_M-Tf) ;
%     alpha_M = [alpha_M , alpha_it_M] ;
%     Twall_M = [Twall_M , Tit_M] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%     Chen - subcooled corrected
%
%     F = 1 ; % Chen factor accounting to vapuor turbulences
%     Re_lo = (G_ch.*D_eq)./ mu_f ;
%     Pr_lo = (Cp_l.*1000.*mu_f)./ k_fl ;
%     h_lo = (0.023.*((Re_lo).^0.8).*((Pr_lo).^0.4)).*(k_fl./D_eq) ;
%     S = (1+(2.53.*(10.^(-6)).*(Re_lo.^1.17))).^(-1) ;
%     Twall2 = Twall - 273.15 ;
%     P_atwall = XSteam('psat_t',Twall2) ;
%     dp = (P_atwall-P_nom).*(10.^6) ;
%     h_nb =@(Tw_tent) 0.00122.*((k_fl.^0.79).*((Cp_l.*1000).^0.45).*(rho_l.^0.49)./((sigma.^0.5).*(mu_f.^0.29).*((h_lg.*1000).^0.24).*(rho_g.^0.24))).*(dp.^0.75).*((Tw_tent-Tsat-273.15).^0.24).*S ;
%     q_CH =@(Tw_tent)  ((0.00122.*((k_fl.^0.79).*((Cp_l.*1000).^0.45).*(rho_l.^0.49)./((sigma.^0.5).*(mu_f.^0.29).*((h_lg.*1000).^0.24).*(rho_g.^0.24))).*(dp.^0.75).*((Tw_tent-Tsat-273.15).^0.24).*S).*(Tw_tent-Tsat-273.15))+(h_lo.*(Tw_tent-Tf)) ;
%     f_CH =@(Tw_tent)  (((0.00122.*((k_fl.^0.79).*((Cp_l.*1000).^0.45).*(rho_l.^0.49)./((sigma.^0.5).*(mu_f.^0.29).*((h_lg.*1000).^0.24).*(rho_g.^0.24))).*(dp.^0.75).*((Tw_tent-Tsat-273.15).^0.24).*S).*(Tw_tent-Tsat-273.15))+(h_lo.*(Tw_tent-Tf)))-q_flux_z ;
%     df_CH =@(Tw_tent) ((0.00122.*((k_fl.^0.79).*((Cp_l.*1000).^0.45).*(rho_l.^0.49)./((sigma.^0.5).*(mu_f.^0.29).*((h_lg.*1000).^0.24).*(rho_g.^0.24))).*(dp.^0.75).*S).*1.24.*((Tw_tent-Tsat-273.15).^0.24))+(h_lo) ;
% 
% 
%     [Tit_CH_vec,itCH]=bisection_method(560,580,nmax,toll,f_CH) ;
%     
%     itCH_vec = [itCH_vec , itCH] ;
%     
%     Tit_CH = Tit_CH_vec(end) ;
%     
%     alpha_it_CH = q_CH(Tit_CH)./(Tit_CH-Tf) ;
%     alpha_CH = [alpha_CH , alpha_it_CH] ;
%     Twall_CH = [Twall_CH , Tit_CH] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % B) Liu-Winterton subcooled correlation
    F = 1 ;
    S = (1+ (0.55.* (Re_lo.^(0.16)))).^-1 ;
    alpha_lo ;
    alpha_pool = 55.*(p_rid.^0.12).*(q_flux_z.^(2/3)).*((-log10(p_rid)).^(-0.55)).*(M.^(-0.5)) ;
    A_bp = alpha_lo./(S.*alpha_pool) ;
    A_qp = q_flux_z./(S.*alpha_pool.*(Tsat+273.15-Tf)) ;

    deltaT = ((Tsat+273.15-Tf)./(1+(A_bp.^2))).*(1+sqrt(1+(1+((A_bp.^2)).*(-1+(A_qp.^2))))) ;

    alpha_it_LW = q_flux_z./deltaT ;
    alpha_LW = [alpha_LW , alpha_it_LW] ;
    
    Tit_LW = Tf + (q_flux_z./alpha_it_LW) ;
    Twall_LW = [Twall_LW , Tit_LW] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Kandlikar subcooled correlation
%     alpha_lo ;
% 
%     f_K =@(Tw_tent) ((1058.*alpha_lo.*((G_ch.*h_lg.*1000).^(-0.7)).*(Tw_tent-Tsat-273.15)).^(1/0.3))-q_flux_z ;
%     
%     [Tit_K_vec,itK]=bisection_method(560,590,nmax,toll,f_K) ;
%     
%     itK_vec = [itK_vec , itK] ;
%     
%     Tit_K = Tit_K_vec(end) ;
%     
%     alpha_it_K = q_flux_z./(Tit_K-Tf) ;
%     alpha_K = [alpha_K , alpha_it_K] ;
%     Twall_K = [Twall_K , Tit_K] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    z_T = [z_T, z2] ;

    x_eq_vec = [x_eq_vec , x_eq] ;

    Temp_bulk = [Temp_bulk , Tf] ;

    % iterative cycle routine
    z1 = z2 ;
    z2 = z2 + 0.001 ;
    z_mean = (z1 + z2)./2 ;
    Ti = Tf ;

    qflux_vec = [qflux_vec , q_flux_z] ;
    
end

zsat = z1-0.001;

Twall_LW(end) = Twall_LW(end-1) ;
alpha_LW(end) = alpha_LW(end-1) ;

Temp_wall1 = [Temp_lo_wall , Twall_JL] ;
Temp_wall2 = [Temp_lo_wall , Twall_T] ;
Temp_wall3 = [Temp_lo_wall , Twall_R] ;
Temp_wall4 = [Temp_lo_wall , Twall_LW] ;

alpha1 = [alpha , alpha_JL] ;
alpha2 = [alpha , alpha_T] ;
alpha3 = [alpha , alpha_R] ;
alpha4 = [alpha , alpha_LW] ;

x_eq_vec = [x_eq_vec(1:end-1), 0] ;


figure(5)
hold on
title('Axial Temperature Profile')
xlabel('Axial Position (m)')
ylabel('Temperature (K)')
hold on
plot(z_T,Temp_bulk,LineWidth=2,LineStyle="-",Color='blue')
hold on
plot(z_T,Temp_wall1,LineWidth=1.25,LineStyle="--",Color=[0.9290 0.6940 0.1250])
hold on
plot(z_T,Temp_wall2,"LineWidth",1.25,"LineStyle","--","Color",[0.8500 0.3250 0.0980])
hold on
plot(z_T,Temp_wall3,"LineStyle","--","LineWidth",1.25,"Color",[0.4660 0.6740 0.1880])
hold on
plot(z_T,Temp_wall4,"LineStyle","--","LineWidth",1.25,"Color",[0.3010 0.7450 0.9330])
hold on
legend('Bulk Temperature','Wall Temperature with JL corr.','Wall Temperature with T corr.','Wall Temperature with R corr.','Wall Temperature with LW corr.')
hold off

figure(6)
hold on
title('Heat Transfer Coefficient Behaviour')
xlabel('Axial Position (m)')
ylabel('Alpha (W/m^2 K)')
hold on
plot(z_T,alpha1,LineWidth=1.5,LineStyle="-",Color=[0.9290 0.6940 0.1250])
hold on
plot(z_T,alpha2,"LineWidth",1.5,"LineStyle","-","Color",[0.8500 0.3250 0.0980])
hold on
plot(z_T,alpha3,"LineStyle","-","LineWidth",1.5,"Color",[0.4660 0.6740 0.1880])
hold on
plot(z_T,alpha4,"LineStyle","--","LineWidth",1.25,"Color",[0.3010 0.7450 0.9330])
hold on
legend('HTC with JL corr.','HTC with T corr.','HTC with R corr.','HTC with LW corr.')
hold off

figure(7)
hold on
title('Vapour Quality x_{eq}')
xlabel('Axial Position (m)')
ylabel('x_eq')
plot(z_T,x_eq_vec,LineWidth=2)
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Pt-3 : Saturated Boiling Analysis %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For the analysis are considered the Kandlikar and Aglar correlations

Cp_l = XSteam('CpL_p',P_xsteam) ;
Cp_v = XSteam('CpV_p',P_xsteam) ;
Tf = Tsat + 273.15 ;
delta_x_eq = 0 ;
x_eq = 0 ;

% Constants for the Aglar correlation
C1 = 3650 ;
C2 = 0.31 ;
m = 0.83 ;
n = 0.80 ;
z = 0.68 ;
p = 1.01 ;
r = 0.82 ;

it_satbol = 0 ;

alpha_kan = [] ;
alpha_ag = [] ;

Twall_kan = [];
Twall_ag = [] ;


while z2 <= 1.8

    % Energy balance in saturated condition
    Q_z = ax_linear_heat_rod(z_mean).*abs(z2 - z1) ;
    q_flux_z = (ax_linear_heat_rod(z_mean).*1000)./(pi.*D_r_out.*0.001) ;
    delta_x_eq = Q_z./(m_ch.*h_lg) ;
    x_eq = x_eq + delta_x_eq ;
    
    qfl_i = (ax_linear_heat_rod(z_mean))./(pi.*D_r_out.*0.001) ;
    q_fl_vec_tot = [q_fl_vec_tot,qfl_i] ;

    % 1 - Kandlikar correlation
    %
    % Determination of the factors in the Kandlikar correlation
    % - 'Co' , convective parameter
    Co = (((1-x_eq)./x_eq).^0.8).*((rho_g./rho_l).^0.5) ;
    % - 'Bo' , boiling parameter
    Bo = q_flux_z./(h_lg.*G_ch.*1000) ;
    % - 'Ffl' , fluid-surface parameter
    Ffl = 1 ; % parameter reccomended by kandlikar for stainless steel structures
    % - 'h_lo' , liquid only single phase HTC. Calculated through the
    % Dittus-Boelter correlation
    Re_lo = (G_ch.*D_eq)./ mu_f ;
    Pr_lo = (Cp_l.*1000.*mu_f)./ k_fl ;
    h_lo = (0.023.*((Re_lo).^0.8).*((Pr_lo).^0.4)).*(k_fl./D_eq) ;

    % determination of boiling and convective dominated heat transfer
    % coefficients
    h_tp_bd = (0.6683.*(Co.^(-0.2)).*h_lo.*((1-x_eq).^0.8))+((1058.*(Bo.^0.7)).*h_lo.*((1-x_eq).^0.8)) ;
    h_tp_cd = (1.136.*(Co.^(-0.9)).*h_lo.*((1-x_eq).^0.8))+(667.2.*(Bo.^0.7).*((1-x_eq).^0.8).*h_lo) ;
    h_tp_vec = [h_tp_bd,h_tp_cd] ;
    % evaluation of the two-phase heat transfer coefficient (the larger)
    h_tp = max(h_tp_vec) ;

    TwallK = (q_flux_z./h_tp)+(Tsat+273.15) ;

    alpha_kan = [alpha_kan , h_tp] ;
    Twall_kan = [Twall_kan , TwallK] ;

    % 2 - Aglar correlation
    %
    % Mc-Adams correlation for the evaluation of the bulk viscosity
    mu_b = ((x_eq./mu_g) + ((1-x_eq)./mu_f)).^(-1) ;

    % calculation of the friction factor for the two-phase and liquid alone
    Re_l_al = (G_ch.*D_eq.*(1-x_eq))./ mu_f ;
    fr_l_al = 0.014 + (0.125./(Re_l_al.^0.32)) ;
    Re_tp = (G_ch.*D_eq)./ mu_b ;
    fr_tp = 0.014 + (0.125./(Re_tp.^0.32)) ;

    % evaluation of the coefficient A for the aglar correlation
    A = (fr_tp./fr_l_al).*(1+((Cp_v./Cp_l).*(x_eq./(1-x_eq)))) ;
    % evaluation of the liquid alone heat transfer coefficient
    h_l_al = (0.023.*((Re_l_al).^0.8).*((Pr_lo).^0.4)).*(k_fl./D_eq) ;
    % determination of the final enhancement factor for the correlation
    E_ps = ((1 + (C1.*(Bo.^m)) + C2 .* ((x_eq./(1-x_eq)).^n).*((rho_l./rho_g).^z)).^p).*(A.^r) ;

    h_aglar = E_ps.*h_l_al ;

    TwallA = (q_flux_z./h_aglar)+(Tsat+273.15) ;

    alpha_ag = [alpha_ag , h_aglar] ;
    Twall_ag = [Twall_ag , TwallA] ;

    % iteration routine
    it_satbol = it_satbol + 1 ;

    x_eq_vec = [x_eq_vec, x_eq] ;
    z_T = [z_T, z2] ;
    Temp_bulk = [Temp_bulk , Tf] ;

    z1 = z2 ;
    z2 = z2 + 0.001 ;
    z_mean = (z1 + z2)./2 ;

end

% reminding the notation due to which each number is referred to a previous
% subcooled boiling casistic :
% 1 - Jens-Lottes
% 2 - Thom
% 3 - Rohsenow
% 4 - Liu-Winterton

alpha1_kan = [alpha1 , alpha_kan] ;
alpha2_kan = [alpha2 , alpha_kan] ;
alpha3_kan = [alpha3 , alpha_kan] ;
alpha4_kan = [alpha4 , alpha_kan] ;

alpha1_ag = [alpha1 , alpha_ag] ;
alpha2_ag = [alpha2 , alpha_ag] ;
alpha3_ag = [alpha3 , alpha_ag] ;
alpha4_ag = [alpha4 , alpha_ag] ;

Temp_wall1_kan = [Temp_wall1 , Twall_kan] ;
Temp_wall2_kan = [Temp_wall2 , Twall_kan] ;
Temp_wall3_kan = [Temp_wall3 , Twall_kan] ;
Temp_wall4_kan = [Temp_wall4 , Twall_kan] ;

Temp_wall1_ag = [Temp_wall1 , Twall_ag] ;
Temp_wall2_ag = [Temp_wall2 , Twall_ag] ;
Temp_wall3_ag = [Temp_wall3 , Twall_ag] ;
Temp_wall4_ag = [Temp_wall4 , Twall_ag] ;

figure(8)
hold on
title('Vapour Quality x_{eq}')
xlabel('Axial Position (m)')
ylabel('x_eq')
plot(z_T,x_eq_vec,LineWidth=2)
hold off

figure(9)
hold on
title('Kandlikar correlation different SCB HTC coupling')
xlabel('Axial Position (m)')
ylabel('Heat Transfer Coefficient \alpha (W/m^2 K)')
plot(z_T,alpha1_kan,"Color",[0.9290 0.6940 0.1250],LineWidth=1.5)
hold on
plot(z_T,alpha2_kan,"Color",[0.8500 0.3250 0.0980],LineWidth=1.5)
hold on
plot(z_T,alpha3_kan,"Color",[0.4660 0.6740 0.1880],LineWidth=1.5)
hold on
plot(z_T,alpha4_kan,"Color",[0.3010 0.7450 0.9330],LineWidth=1.5)
hold on
legend('Jens-Lottes','Thom','Rohsenow','Liu-Winterton')
hold off

figure(10)
hold on
title('Kandlikar correlation different SCB Wall T coupling')
xlabel('Axial Position (m)')
ylabel('Temperature (K)')
plot(z_T,Temp_wall1_kan,"Color",[0.9290 0.6940 0.1250],LineWidth=1.5)
hold on
plot(z_T,Temp_wall2_kan,"Color",[0.8500 0.3250 0.0980],LineWidth=1.5)
hold on
plot(z_T,Temp_wall3_kan,"Color",[0.4660 0.6740 0.1880],LineWidth=1.5)
hold on
plot(z_T,Temp_wall4_kan,"Color",[0.3010 0.7450 0.9330],LineWidth=1.5)
hold on
plot(z_T,Temp_bulk,"Color","k",LineWidth=1.5)
hold on
legend('Jens-Lottes','Thom','Rohsenow','Liu-Winterton','Bulk')
hold off

figure(11)
hold on
title('Aglar correlation different SCB HTC coupling')
xlabel('Axial Position (m)')
ylabel('Heat Transfer Coefficient \alpha (W/m^2 K)')
plot(z_T,alpha1_ag,"Color",[0.9290 0.6940 0.1250],LineWidth=1.5)
hold on
plot(z_T,alpha2_ag,"Color",[0.8500 0.3250 0.0980],LineWidth=1.5)
hold on
plot(z_T,alpha3_ag,"Color",[0.4660 0.6740 0.1880],LineWidth=1.5)
hold on
plot(z_T,alpha4_ag,"Color",[0.3010 0.7450 0.9330],LineWidth=1.5)
hold on
legend('Jens-Lottes','Thom','Rohsenow','Liu-Winterton')
hold off

figure(12)
hold on
title('Aglar correlation different SCB Wall T coupling')
xlabel('Axial Position (m)')
ylabel('Temperature (K)')
plot(z_T,Temp_wall1_ag,"Color",[0.9290 0.6940 0.1250],LineWidth=1.5)
hold on
plot(z_T,Temp_wall2_ag,"Color",[0.8500 0.3250 0.0980],LineWidth=1.5)
hold on
plot(z_T,Temp_wall3_ag,"Color",[0.4660 0.6740 0.1880],LineWidth=1.5)
hold on
plot(z_T,Temp_wall4_ag,"Color",[0.3010 0.7450 0.9330],LineWidth=1.5)
hold on
plot(z_T,Temp_bulk,"Color","k",LineWidth=1.5)
hold on
legend('Jens-Lottes','Thom','Rohsenow','Liu-Winterton','Bulk')
hold off

% plotting the two extremes that could be obtained using the two kandlikar
% or aglar correlations

figure(13)
hold on
title('Min and Max wall HTC behaviours that are predicted')
xlabel('Axial Position (m)')
ylabel('Heat Transfer Coefficient \alpha (W/m^2 K)')
plot(z_T,alpha1_ag,"Color","b",LineWidth=1.5)
hold on
plot(z_T,alpha4_kan,"Color","r",LineWidth=1.5)
hold on
legend('Jens-Lottes , Aglar','Liu-Winterton , Kandlikar')
hold off

figure(14)
hold on
title('Min and Max SCB Wall T coupling')
xlabel('Axial Position (m)')
ylabel('Temperature (K)')
plot(z_T,Temp_wall1_ag,"Color","b",LineWidth=1.5)
hold on
plot(z_T,Temp_wall4_kan,"Color","r",LineWidth=1.5)
hold on
plot(z_T,Temp_bulk,"Color","k",LineWidth=1.5)
hold on
legend('Jens-Lottes , Aglar','Liu-Winterton , Kandlikar','Bulk')
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Pt-4 : CHF Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For the CHF analysis are studied two different cases :
% 1 - CISE-4 correlation for dryout
% 2 - Groenveld Lookup tables for both DNB and Dryout

% 1 - CISE 4
% critical quality for dryout calculation
Lb = 1.8 - zsat ;
p_crit_MPa = p_crit./1000 ;
p_rid ;
G_ref = 3375.*((1-p_rid).^3) ;
a = 0;
ax_linear_heat_rod =@(zz) q0_l_peak.*cos(pi.*zz./L_rod);
Power_nom = trap_int_powerpurp(-1.8,1.8,3600,ax_linear_heat_rod) ;

if G_ch < G_ref
    a = (1+(1.481.*(10.^-4).*((1-p_rid).^-3).*G_ch)).^-1 ;
else
    a = (1-p_rid)./((G_ch./1000).^(1/3)) ;
end

b = 0.199.*(((p_rid.^-1)-1).^0.4) .*G_ch.*(D_eq.^1.4) ;

% critical quality
xcr_cise = a .* (Lb./(Lb+b)) ;
% linear heat peak to have the critical quality at the exit point: +L_rod
q0_lin_cr_cise = m_ch.*(h_satl-h_in+(h_lg.*xcr_cise)).*(pi./(2.*L_rod)) ;
% critical linear heat to have dryout
q_cise_prof =@(x)  q0_lin_cr_cise.*cos(pi.*(x./L_rod)) ;
% critical heat power of the channel
Q_cr_cise = trap_int_powerpurp(-1.8,1.8,3600,q_cise_prof) ;
% Critical Power Ratio in nominal condition
CPR_cise = Q_cr_cise./Power_nom ;

disp(xcr_cise)
disp(Q_cr_cise)
disp(CPR_cise)

% What done up to now can be extended to see with different operational
% parameters what happens in steady state conditions

% definition of the range to investigate (both varying the mass flow and 
% the linear heat of the rod)
q0_l_vec = [10:0.2:85] ; % prev step 3
m_ch_vec = [0.1:0.001:0.475] ; % prev step 0.015


LL = length(q0_l_vec) ;
KK = length(m_ch_vec) ;

CPR_matrix = ones(LL,KK) ;
x_cr_matrix = ones(LL,KK) ;

Ti = T_in+273.15 ;

% test of evaluation of the critical condition over which the water
% saturation temperature is not reached within the length of the rod
% f_z_i =@(z) ((Tsat+273.15-Ti).*(0.476.*Cp_l))-(10.*(L_rod./pi).*(sin(pi.*z./L_rod)+1)) ;
% zzz = linspace(-1.8,1.8);
% plot(zzz,f_z_i(zzz))

ii = 0;
jj = 0;

%BWR_parameters ;

for q_i = q0_l_vec

    ii = ii +1;
    jj = 0 ;
    for m_i = m_ch_vec
        
        jj = jj +1;

        % axial profile of the total heat (kW) given by the rod
        Q_z_i =@(z) q_i.*(L_rod./pi).*(sin(pi.*z./L_rod)+1) ;
        % function that must be zero in order to find the point at whose
        % height is reached the saturated condition. It is given by the
        % difference between the heat given by the rod and the heat that
        % the coolant absorbs ib steady state condition.
        f_z_i =@(z) ((Tsat+273.15-Ti).*(m_i.*Cp_l))-(q_i.*(L_rod./pi).*(sin(pi.*z./L_rod)+1)) ;
        
        % bisection method to find the length at which the saturated
        % condition is reached.
        [ris,xdif,fx,it_cise] = bisection_method(-1.8,1.8,nmax,toll,f_z_i) ;

        zsat = ris(end) ;

        Lb_i = 1.8 - zsat ;

        G_i = m_i./Area_flow ;

        ax_linear_heat_rod_i =@(zz) q_i.*cos(pi.*zz./L_rod);
        Power_nom_i = trap_int_powerpurp(-1.8,1.8,3600,ax_linear_heat_rod_i) ;

        if G_i < G_ref
            a = (1+(1.481.*(10.^-4).*((1-p_rid).^-3).*G_i)).^-1 ;
        else
            a = (1-p_rid)./((G_i./1000).^(1/3)) ;
        end
        
        b = 0.199.*(((p_rid.^-1)-1).^0.4) .*G_i.*(D_eq.^1.4) ;

        xcr_i = a .* (Lb_i./(Lb_i+b)) ;
        q0_lin_cr_cise_i = m_i.*(h_satl-h_in+(h_lg.*xcr_i)).*(pi./(2.*L_rod)) ;
        Q_cise_prof_i =@(x)  q0_lin_cr_cise_i.*cos(pi.*(x./L_rod)) ;
        Q_cr_cise_i = trap_int_powerpurp(-1.8,1.8,3600,Q_cise_prof_i) ;
        
        CPR_cise_i = Q_cr_cise_i./Power_nom_i ;

        CPR_matrix(ii,jj) = CPR_cise_i ;
        x_cr_matrix(ii,jj) = xcr_i ;
        
    end
end

tiledlayout(1,1)
[XX1,YY1] = meshgrid(q0_l_vec,m_ch_vec) ;

nexttile
s_CPR = surf(XX1,YY1,CPR_matrix,x_cr_matrix) ;
hold on
c = colorbar;
c.Label.String = 'Critical Quality ''x_{cr}'' ' ;
xlabel('Linear Heat Peak ''q_0'' (kW/m)')
ylabel('Channel Mass Flow ''m_{ch}'' (kg/s)')
zlabel('Critical Power Ratio ''CPR''')
s_CPR.EdgeColor = 'none';
colormap("turbo")

% nextile
% s_xcr = surf(XX1,YY1,x_cr_matrix,CPR_matrix) ;
% hold on
% c = colorbar;
% c.Label.String = 'Critical Power Ratio ''CPR'' ' ;
% xlabel('Linear Heat Peak ''q_0'' (kW/m)')
% ylabel('Channel Mass Flow ''m_{ch}'' (kg/s)')
% zlabel('Critical Quality ''x_{cr}''')
% s_CPR.EdgeColor = 'none';
% colormap("turbo")

% figure(16)
% [XX2,YY2] = meshgrid(q0_l_vec,m_ch_vec) ;
% x_cr_final = x_cr_matrix.* (XX2.*0 + YY2.*0 +1);
% surf(XX2, YY2, x_cr_final)
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 - Groenveld LUT
P_inlet = P_nom.*(10.^6) ;
x_inlet = (h_in - h_satl)./h_lg ;
A_heated = L_rod.*pi.*D_r_out./1000 ;
A_fl = Area_flow ;
h_vapzat = h_lg.*1000 ;
toggle_K = [1; 1; 0; 1; 0; 0; 0; 0];
D_htr = 11.2./1000 ;
pitch = Pitch_rodrod./1000;
K_g = 0;
toggle_vis = ['on'];
toggle_plotsave = 1 ;
plotfilepath = ['/Users/davidemarchesi/Desktop/backup multiphase/prove groenveld'];

m_ch_vec = [0.1:0.020:0.45] ;
m_ii = 0 ;

CHF_groen = [] ;
x_cr_gro = [] ;

for m_ii = m_ch_vec

    G_ch_i = m_ii./A_fl ;

    [CHF_LUT_i,CHF_LUT_Avg,L_LUT,x_LUT_local_i,K_Factors] = ...
    Groeneveld_LUT(G_ch_i,P_inlet,x_inlet,A_heated,A_fl,L_rod,h_vapzat,...
    toggle_K,D_eq,D_htr,pitch,K_g,rho_l,rho_g,toggle_vis,toggle_plotsave,plotfilepath) ;

    CHF_groen = [CHF_groen , CHF_LUT_i] ;
    x_cr_gro = [x_cr_gro , x_LUT_local_i] ;
end

figure()
hold on
xlabel('Quality ''x_{eq}''')
ylabel('Heat Flux ''q'''' (kW/m^2)')
plot(x_cr_gro,CHF_groen,LineWidth=2,Color="r",LineStyle="-") ;
hold on
plot(x_eq_vec,q_fl_vec_tot,LineWidth=2,Color="b",LineStyle="-.") ;
hold on
legend('Groenveld LUT Critical Heat Flux','Heat Flux in operational conditions')
hold off









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % colori pastello per i grafici :
% blu            [0 0.4470 0.7410]
% arancione      [0.8500 0.3250 0.0980]
% giallo         [0.9290 0.6940 0.1250]
% viola          [0.4940 0.1840 0.5560]
% verde          [0.4660 0.6740 0.1880]
% azzurro        [0.3010 0.7450 0.9330]
% rosso          [0.6350 0.0780 0.1840]