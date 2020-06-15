%% Heat Transfer Model for Firefighters' Protective Clothing

clear all; close all; clc; %#ok<CLALL>

%% Code notation

% Subscript 1 refers to layer 1: Shell
% Subscript 2 refers to layer 2: Moisture Barrier
% Subscript 3 refers to layer 3: Thermal Layer

%% Table 1: Physical characteristics of fabric layers (20 C)

% Thickness (m)
th_1 = 0.00082; % +/- 0.00007
th_2 = 0.00055; % +/- 0.00005
th_3 = 0.0035;  % +/- 0.0004

% Density (kg/m^3)
rho_1 = 310;    % +/- 24
rho_2 = 800;    % +/- 60
rho_3 = 72;     % +/- 7

% Conductivity (W/mK)
k_1 = 0.047; 
k_2 = 0.012;
k_3 = 0.038;

% Specific Heat (J/kgK)
c_1 = 1300; 
c_2 = 2010;
c_3 = 700;

% Transmissivity 
tau_1 = 0.044;
tau_2 = 0.005; 
tau_3 = 0.0012;

% Reflectivity
r_1 = 0.090;
r_2 = 0.017;
r_3 = 0.002; 

%% Additional physical properties needed

% Surface Area (m^2)
A = 0.255^2;

% Seperation Between Layers (m)
sep = 0.001;

% Ambient Temperature
T_inf = 29.3 + 273.15;  % K

% Stefan Boltzmann Constant
sigma = 5.669*(10^-8);  % W/m^2K^4

% Emissivities
eps_1 = 1 - r_1 - tau_1;
eps_2 = 1 - r_2 - tau_2;
eps_3 = 1 - r_3 - tau_3;

% Flame Heat Transfer Rate
q = 2500;   % W/m^2 (radiative flux of typical flame)
Q = q*A;    % W (cooresponding heat transfer rate)

% Compute radiative coefficient (parallel wall case)
h_rad_front = @(T1) sigma*eps_1*(T_inf^2 + T1^2)*(T_inf + T1);
h_rad_12 = @(T1,T2) sigma*((eps_1*eps_2)/(eps_1+eps_2+(eps_1*eps_2)))*(T1^2+T2^2)*(T1+T2);
h_rad_23 = @(T2,T3) sigma*((eps_2*eps_3)/(eps_2+eps_3+(eps_2*eps_3)))*(T2^2+T3^2)*(T2+T3);
h_rad_back = @(T3)  sigma*eps_3*(T3^2 + T_inf^2)*(T3 + T_inf);

% Compute average conductive coefficient
k_a = @(T1,T2) (k_air(T1) + k_air(T2))/2;

%% Impedence Computation

% Conductive Resistances (in material)
R_k1 = th_1/(k_1*A);    % K/W
R_k2 = th_2/(k_2*A);    % K/W
R_k3 = th_3/(k_3*A);    % K/W

% Conductive Resistances (between materials)(Initial Guesses)
R_k12 = sep/(k_a(T_inf-273.15,T_inf-273.15)*A);  % K/W
R_k23 = sep/(k_a(T_inf-273.15,T_inf-273.15)*A);  % K/W

% Convective Resistances (Initial Guesses)
R_h_front = 1/(h(T_inf,T_inf+10)*A); % K/W
R_h_back = 1/(h(T_inf,T_inf+10)*A);  % K/W

% Radiative Resistances (Initial Guesses)
R_rad_front = 1/(h_rad_front(T_inf)*A);      % K/W
R_rad_12 = 1/(h_rad_12(T_inf+10,T_inf)*A);   % K/W
R_rad_23 = 1/(h_rad_23(T_inf+10,T_inf)*A);   % K/W
R_rad_back = 1/(h_rad_back(T_inf)*A);        % K/W

% Capacitances
s = tf('s');
C_1 = c_1*rho_1*th_1*A; % J/K
C_2 = c_2*rho_2*th_2*A; % J/K
C_3 = c_3*rho_3*th_3*A; % J/K
ZC_1 = 1/(s*C_1); ZC_2 = 1/(s*C_2); ZC_3 = 1/(s*C_3);

%% Equivilant Resistances

R_eq0 = (R_k1/2) + (R_rad_front*R_h_front/(R_rad_front+R_h_front)); % K/W
R_eq1 = (R_k1/2) + (R_k2/2) + (R_k12*R_rad_12/(R_k12+R_rad_12));    % K/W
R_eq2 = (R_k2/2) + (R_k3/2) + (R_k23*R_rad_23/(R_k23+R_rad_23));    % K/W
R_eq3 = (R_k3/2) + (R_rad_back*R_h_back/(R_rad_back+R_h_back));     % K/W

%% Nodal Analysis

% Admittance Matrix
B = [((1/R_eq0)+(1/R_eq1)+(1/ZC_1)),(-1/R_eq1),0;
    (-1/R_eq1),((1/R_eq1)+(1/R_eq2)+(1/ZC_2)),(-1/R_eq2);
    0,(-1/R_eq2),((1/R_eq2)+(1/R_eq3)+(1/ZC_3))];

% Forcing Vector
q_dot = [Q + T_inf*((1/R_eq0)+(1/ZC_1));T_inf/ZC_2;T_inf*((1/R_eq3)+(1/ZC_3))]/s;

% Find Temperature in Laplace Domain
T = B\q_dot;

% Convert to Time Domain
time = linspace(0,300,1000);  
t_K_init = impulse(T,time);     % K
t_C_init = t_K_init - 273.15;   % C

% Initial guesses response
t_1_init = t_C_init(:,1);   % C
t_2_init = t_C_init(:,2);   % C
t_3_init = t_C_init(:,3);   % C

%% Plot Response Using Initial Guesses

figure(1)
plot(time,t_1_init,time,t_2_init,time,t_3_init,'LineWidth',1.4)
title('Temperature of Layers Heating Up - Initial Guess','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (C)','interpreter','latex','FontSize',12)
legend_init = legend('Shell Layer (Closest to Flame)','Moisture Barrier (Middle Layer)',...
    'Thermal Layer (Furthest from Flame)','Location','northwest');
set(legend_init,'Interpreter','latex');
grid on
grid minor

%% Iterative Method for Finding Convection and Radiation

% Initialize the variables used while iterating
t_1_old = t_K_init(end,1);
t_2_old = t_K_init(end,2);
t_3_old = t_K_init(end,3);

% Initialize end condition
delta_t1 = 10;

% Check number of iterations
iterate1 = 0;

% While the temperature is still changing significantly
while delta_t1 > 0.01
    
    % Redefine resistances
    
    % Conductive Resistances (between materials)
    R_k12 = sep/(k_a(t_1_old-273.15,t_2_old-273.15)*A);  % K/W
    R_k23 = sep/(k_a(t_2_old-273.15,t_3_old-273.15)*A);  % K/W
    % Convective Resistances
    R_h_front = 1/(h(T_inf,t_1_old)*A); % K/W
    R_h_back = 1/(h(t_3_old,T_inf)*A);  % K/W
    % Radiative Resistances
    R_rad_front = 1/(h_rad_front(t_1_old)*A);      % K/W
    R_rad_12 = 1/(h_rad_12(t_1_old,t_2_old)*A);    % K/W
    R_rad_23 = 1/(h_rad_23(t_2_old,t_3_old)*A);    % K/W
    R_rad_back = 1/(h_rad_back(t_3_old)*A);        % K/W  
    % Equivilant Resistances
    R_eq0 = (R_k1/2) + (R_rad_front*R_h_front/(R_rad_front+R_h_front)); % K/W
    R_eq1 = (R_k1/2) + (R_k2/2) + (R_k12*R_rad_12/(R_k12+R_rad_12));    % K/W
    R_eq2 = (R_k2/2) + (R_k3/2) + (R_k23*R_rad_23/(R_k23+R_rad_23));    % K/W
    R_eq3 = (R_k3/2) + (R_rad_back*R_h_back/(R_rad_back+R_h_back));     % K/W
    
    % Redo nodal analysis
    
    % Admittance Matrix
    B = [((1/R_eq0)+(1/R_eq1)+(1/ZC_1)),(-1/R_eq1),0;
    (-1/R_eq1),((1/R_eq1)+(1/R_eq2)+(1/ZC_2)),(-1/R_eq2);
    0,(-1/R_eq2),((1/R_eq2)+(1/R_eq3)+(1/ZC_3))];
    % Forcing Vector
    q_dot = [Q + T_inf*((1/R_eq0)+(1/ZC_1));T_inf/ZC_2;T_inf*((1/R_eq3)+(1/ZC_3))]/s;
    % Find Temperature in Laplace Domain
    T = B\q_dot;
    % Convert to time domain
    t_K_it = impulse(T,time);     % K
    t_C_it = t_K_it - 273.15;     % C
    
    % Check end conditions
    
    % Update iterated temperatures
    t_1_new = t_K_it(end,1);
    t_2_new = t_K_it(end,2);
    t_3_new = t_K_it(end,3);
    
    % Update end condition
    delta_t1 = abs(t_1_new - t_1_old);
    
    % Update old temperatures
    t_1_old = t_1_new;
    t_2_old = t_2_new;
    t_3_old = t_3_new;
    
    % Update number of iterations
    iterate1 = iterate1 + 1;
end

% Iterated response
t_1_iterated = t_C_it(:,1);
t_2_iterated = t_C_it(:,2);
t_3_iterated = t_C_it(:,3);

%% Plot Response Using Iterated Values

figure(2)
plot(time,t_1_iterated,time,t_2_iterated,time,t_3_iterated,'LineWidth',1.4)
title('Temperature of Layers Heating Up - Iterated Values','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (C)','interpreter','latex','FontSize',12)
legend_final = legend('Shell Layer (Closest to Flame)','Moisture Barrier (Middle Layer)',...
    'Thermal Layer (Furthest from Flame)','Location','northwest');
set(legend_final,'Interpreter','latex');
grid on
grid minor

%% Removal of Heat Source - Initial Guess

% Admittance Matrix
B2 = [((1/R_eq0)+(1/R_eq1)+(1/ZC_1)),(-1/R_eq1),0;
(-1/R_eq1),((1/R_eq1)+(1/R_eq2)+(1/ZC_2)),(-1/R_eq2);
0,(-1/R_eq2),((1/R_eq2)+(1/R_eq3)+(1/ZC_3))];

% Forcing Vector
q_dot2 = [(T_inf/R_eq0)+(t_1_old/ZC_1);t_2_old/ZC_2;(T_inf/R_eq3)+(t_3_old/ZC_3)]/s;

% Find Temperature in Laplace Domain
T2 = B2\q_dot2;

% Convert to Time Domain 
t_K_init2 = impulse(T2,time);     % K
t_C_init2 = t_K_init2 - 273.15;   % C

% Initial guesses response
t_1_init2 = t_C_init2(:,1);   % C
t_2_init2 = t_C_init2(:,2);   % C
t_3_init2 = t_C_init2(:,3);   % C

%% Plot Response Using Initial Guesses

figure(3)
plot(time+300,t_1_init2,time+300,t_2_init2,time+300,t_3_init2,'LineWidth',1.4)
title('Temperature of Layers Cooling Down - Initial Guess','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (C)','interpreter','latex','FontSize',12)
legend_init2 = legend('Shell Layer (Closest to Flame)','Moisture Barrier (Middle Layer)',...
    'Thermal Layer (Furthest from Flame)','Location','northeast');
set(legend_init2,'Interpreter','latex');
grid on
grid minor

%% Iterative Method for Finding Convection and Radiation

% Initialize the variables used while iterating
t_1_old2 = t_K_init2(end,1);
t_2_old2 = t_K_init2(end,2);
t_3_old2 = t_K_init2(end,3);

% Initialize end condition
delta_t2 = 10;

% Check number of iterations
iterate2 = 0;

% While the temperature is still changing significantly
while delta_t2 > 0.01
    
    % Redefine resistances
    
    % Conductive Resistances (between materials)
    R_k12 = sep/(k_a(t_1_old2-273.15,t_2_old2-273.15)*A);  % K/W
    R_k23 = sep/(k_a(t_2_old2-273.15,t_3_old2-273.15)*A);  % K/W
    % Convective Resistances
    R_h_front = 1/(h(T_inf,t_1_old2)*A); % K/W
    R_h_back = 1/(h(t_3_old2,T_inf)*A);  % K/W
    % Radiative Resistances
    R_rad_front = 1/(h_rad_front(t_1_old2)*A);      % K/W
    R_rad_12 = 1/(h_rad_12(t_1_old2,t_2_old2)*A);    % K/W
    R_rad_23 = 1/(h_rad_23(t_2_old2,t_3_old2)*A);    % K/W
    R_rad_back = 1/(h_rad_back(t_3_old2)*A);        % K/W  
    % Equivilant Resistances
    R_eq0 = (R_k1/2) + (R_rad_front*R_h_front/(R_rad_front+R_h_front)); % K/W
    R_eq1 = (R_k1/2) + (R_k2/2) + (R_k12*R_rad_12/(R_k12+R_rad_12));    % K/W
    R_eq2 = (R_k2/2) + (R_k3/2) + (R_k23*R_rad_23/(R_k23+R_rad_23));    % K/W
    R_eq3 = (R_k3/2) + (R_rad_back*R_h_back/(R_rad_back+R_h_back));     % K/W
    
    % Redo nodal analysis
    
    % Admittance Matrix
    B2 = [((1/R_eq0)+(1/R_eq1)+(1/ZC_1)),(-1/R_eq1),0;
    (-1/R_eq1),((1/R_eq1)+(1/R_eq2)+(1/ZC_2)),(-1/R_eq2);
    0,(-1/R_eq2),((1/R_eq2)+(1/R_eq3)+(1/ZC_3))];
    % Forcing Vector
    q_dot2 = [(T_inf/R_eq0)+(t_1_old/ZC_1);t_2_old/ZC_2;(T_inf/R_eq3)+(t_3_old/ZC_3)]/s;
    % Find Temperature in Laplace Domain
    T2 = B2\q_dot2;
    % Convert to time domain
    t_K_it2 = impulse(T2,time);     % K
    t_C_it2 = t_K_it2 - 273.15;     % C
    
    % Check end conditions
    
    % Update iterated temperatures
    t_1_new2 = t_K_it2(end,1);
    t_2_new2 = t_K_it2(end,2);
    t_3_new2 = t_K_it2(end,3);
    
    % Update end condition
    delta_t2 = abs(t_1_new2 - t_1_old2);
    
    % Update old temperatures
    t_1_old2 = t_1_new2;
    t_2_old2 = t_2_new2;
    t_3_old2 = t_3_new2;
    
    % Update number of iterations
    iterate2 = iterate2 + 1;
end

% Iterated response
t_1_iterated2 = t_C_it2(:,1);
t_2_iterated2 = t_C_it2(:,2);
t_3_iterated2 = t_C_it2(:,3);

%% Plot Response Using Iterated Values

figure(4)
plot(time+300,t_1_iterated2,time+300,t_2_iterated2,time+300,t_3_iterated2,'LineWidth',1.4)
title('Temperature of Layers Cooling Down - Iterated Values','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (C)','interpreter','latex','FontSize',12)
legend_final2 = legend('Shell Layer (Closest to Flame)','Moisture Barrier (Middle Layer)',...
    'Thermal Layer (Furthest from Flame)','Location','northeast');
set(legend_final2,'Interpreter','latex');
grid on
grid minor

%% Full Turnout Coat Simulation

% Concatinate two parts of simulation
time_full = [time, time+300];
t1_full = [t_1_iterated;t_1_iterated2];
t2_full = [t_2_iterated;t_2_iterated2];
t3_full = [t_3_iterated;t_3_init2];

% Plot the full response
figure(5)
plot(time_full,t1_full,time_full,t2_full,time_full,t3_full,'LineWidth',2)
title('Full Turnout Coat Simulation','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (C)','interpreter','latex','FontSize',12)
legend_full = legend('Shell Layer (Closest to Flame)','Moisture Barrier (Middle Layer)',...
    'Thermal Layer (Furthest from Flame)','Location','northeast');
set(legend_full,'Interpreter','latex');
grid on
grid minor
xlim([0,600]);
ylim([0,180]);

