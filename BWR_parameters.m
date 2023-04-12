%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Nine Mile Point 2 G.E. BWR5 with GE11 Fuel %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Definition of the Geometrical and General Physical parameters

Q_th = 3323 ; % (MWth) total thermal power
P_nom = 7.14 ; % (MPa) nominal pressure
P_xsteam = 71.4 ; % (bar) nominal pressure in bar for XSteam usage
T_in_core = 278.3 ; % (°C) core inlet temperature
T_in = 278.3 ; % (°C) core inlet temperature
x_out_avg = 14.6 ; % (%) core average exit quality
m_core = 13671 ; % (kg/s) core coolant total flow rate
Pitch_rodrod = 14.37 ; % (mm) pitch rod to rod
D_r_out = 11.20 ; % (mm) fuel rod outside parameter
th_clad = 0.71 ; % (mm) cladding thickness
th_gap = 0.09 ; % (mm) fuel cladding gap in cold conditions
D_fp = 9.60 ; % (mm) fuel pellet diameter
L_fp = 10 ; % (mm) fuel pellet length
L_tot_rod = 4.10 ; % (m) total fuel rod length
L_rod = 3.6 ; % (m) heated length of the fuel rod (approx. from 3.588m)
q0_l_peak = 47.24 ; % (kW/m) core peak LHGR q'
q0_l_avg = 17.6 ; % (kW/m) core average LHGR q'
mu_f = 9.075e-5 ; % (kg/m*s) viscosity of the fluid
mu_g = 1.902e-5 ; % (kg/m*s) viscosity of the gas


% Geometrical calculations concerning one single rod (rivedi pag.716)

Per_rod_wet = pi .* D_r_out./1000 ; % (m) perimeter of a single rod
Area_rod = (pi .* (D_r_out).^2)/(4.*1000) ; % (m^2) area of a single rod
Area_flow = 1.42*(10^-4);% flow area (m^2)
m_ch = 0.2229 ; % mass flow (kg/s)
G_ch = m_ch/Area_flow ; % mass flux (kg/m^2 s)
D_eq = (4.*Area_flow)./Per_rod_wet ; % (m) Hydraulic equivalent diameter 


% All the TDN properties can be obtained using the XSteam functions or
% other tabulated values.
h_lg = 1499.63 + ((71.4-71).*(1493.35-1499.63)) ; % (kJ/kg) latent heat at 71.4 bar


% Calculation of the axial behaviour of the linear heat.
z_ax = -1.8 : 0.001 : 1.8 ;

% Note that the profile is extrapolated considering the length of the
% heated rod equal to the length of the neutronics
ax_linear_heat_rod =@(zz) q0_l_peak.*cos(pi.*zz./L_rod); % (kW/m)
integral_linear_heat = trap_int_powerpurp(-1.8,1.8,3600,ax_linear_heat_rod);
average_ax_linear_heat_rod =@(zz) integral_linear_heat./L_rod + (zz.*0) ;

% The values of the real axial point can be easily
% extrapolated as : Z_real = z + L_rod/2 .
% Thus are admitted only values of z s.t. : -1.8 < z < 1.8

figure(1)
hold on
title('Axial Linear Heat')
xlabel('Rod Length (m)')
ylabel('Linear Heat (kW/m)')
plot(z_ax,ax_linear_heat_rod(z_ax),"LineStyle","-","LineWidth",2,"Color",'blue')
hold on
plot(z_ax, average_ax_linear_heat_rod(z_ax),"LineStyle","--","LineWidth",1.5)
legend('Axial Linear Heat q''(z)', 'Mean Axial Linear Heat q''_{mean}',Location='northeast')
hold off
