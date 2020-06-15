%% Midterm Resubmission - ME342 - March 28, 2020

clear all; close all; clc;

%% Define Parameters

% Dimensions
t = 0.001; % m
w = 0.06; % m
l = 0.04; % m
g = 0.00275; % m
b = 0.006; % m
F = 16; % number of fins

% Thermal Properties
h = 50; % W/m^2K
% Material chosen by commenting out some variables
% Copper
k = 400; % W/mK
rho = 8900; % kg/m^3
c = 390; % J/kgK
    % % Aluminium
    % k = 200; % W/mK
    % rho = 2700; % kg/m^3
    % c = 910; % J/kgK
    % % Steel
    % k = 20; % W/mK
    % rho = 7900; % kg/m^3
    % c = 490; % J/kgK

% System Conditions
Q = 250/F; % W
T_inf = 25+273.15; % K

%% Impedances

% Note: These are calculated for half of one fin.

% Resistances
R_k = (2*l)/(k*w*t); % K/W
R_h = 1/(h*l*w); % K/W
R_h_base = 2/(h*w*g); % K/W
R_h_end = 2/(h*w*t); % K/W

% Capacitances
s = tf('s');
C = (c*l*w*t*rho)/2; % J/K
C_base = (c*g*t*b*rho)/4; % J/K
Z_C = 1/(s*C);
Z_C_base = 1/(s*C_base);

%% Find Initial Conditions 

% Define terms of admittance matrix in steady-state
mho_0_11 = (1/R_h_base) + (6/R_k);
mho_0_22 = (9/R_k) + (3/R_h);
mho_0_33 = (6/R_k) + (3/R_h);
mho_0_44 = (3/R_k) + (3/R_h) + (1/((R_k/6) + R_h_end));
mho_0_12 = -6/R_k;
mho_0_23 = -3/R_k; mho_0_34 = -3/R_k;

% Define admittance matrix
mho_0 = [mho_0_11, mho_0_12, 0, 0; 
         mho_0_12, mho_0_22, mho_0_23, 0;
         0, mho_0_23, mho_0_33, mho_0_34;
         0, 0, mho_0_34, mho_0_44];

% Define terms of forcing vector
q_0_1 = (Q/2) + (T_inf/R_h_base);
q_0_2 = 3*T_inf/R_h;
q_0_3 = 3*T_inf/R_h;
q_0_4 = (3*T_inf/R_h) + (T_inf/((R_k/6) + R_h_end));

% Define forcing vector
q_0 = [q_0_1; q_0_2; q_0_3; q_0_4];

% Find initial temperatures
T_0 = mho_0\q_0;

% Define variable names for intial temperatures
T_00 = T_0(1);
T_01 = T_0(2);
T_02 = T_0(3);
T_03 = T_0(4);

%% Find Transient Response

% Define terms of admittance matrix
mho_11 = (6/R_k) + (1/R_h_base) + (1/Z_C_base);
mho_22 = (9/R_k) + (3/R_h) + (1/(3*Z_C));
mho_33 = (6/R_k) + (3/R_h) + (1/(3*Z_C));
mho_44 = (3/R_k) + (3/R_h) + (1/((R_k/6) + R_h_end)) + (1/(3*Z_C));
mho_12 = -6/R_k; mho_23 = -3/R_k; mho_34 = -3/R_k;

% Define admittance matrix
mho = [mho_11, mho_12, 0, 0;
       mho_12, mho_22, mho_23, 0;
       0, mho_23, mho_33, mho_34;
       0, 0, mho_34, mho_44];

% Define terms of forcing vector   
q_1 = (T_00/Z_C_base) + (T_inf/R_h_base);
q_2 = (T_01/(3*Z_C)) + (3*T_inf/R_h);
q_3 = (T_02/(3*Z_C)) + (3*T_inf/R_h);
q_4 = (T_03/(3*Z_C)) + (3*T_inf/R_h) + T_inf/((R_k/6)+R_h_end);

% Define forcing vector
q = [q_1; q_2; q_3; q_4]/s;

% Find temperatures at nodes
T = mho\q;

% Convert to time domain and switch to Celcius
time = linspace(0,25,10000);
t_0 = impulse(T(1),time) - 273.15;
t_1 = impulse(T(2),time) - 273.15;
t_2 = impulse(T(3),time) - 273.15;
t_3 = impulse(T(4),time) - 273.15;

%% Plot Temperature of Nodes

% Color
four_colors = [0, 229, 59; 0, 229, 201; 0, 152, 229; 31, 0, 229]/255;
set(groot,'defaultAxesColorOrder',four_colors)

figure(1)
plot(time,t_0,time,t_1,time,t_2,time,t_3,'LineWidth',1)
title('Temperature of Steel Nodes','interpreter','latex','FontSize',14)
xlabel('Time (s)','interpreter','latex','FontSize',12)
ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
key = legend('Base Plate','Node 1','Node 2','Node 3');
set(key,'Interpreter','latex');
grid on
grid minor
xlim([0,time(end)]);
ylim([24,51]);



