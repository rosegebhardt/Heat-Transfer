%% Rose Gebhardt and Harris Paspuleti -- February 12, 2020 -- Homework 3 Question 1

clear all; close all; clc;

%% Parameters Used in Both Problems

s = tf('s');
time_1 = linspace(0,6,1000);
time_2 = linspace(0,2,1000);

% Define material properties
rho_st = 8000; % kg/m^3
rho_cu = 8960; % kg/m^3
c_st = 502.416; % J/kgK
c_cu = 376.812; % J/kgK
k_st = 17; % W/mK
k_cu = 372; % W/mK

% Citations
    % https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm
    % https://hypertextbook.com/facts/2004/KarenSutherland.shtml
    % https://amesweb.info/Materials/Density_of_Copper.aspx

% Define dimensions
A = 0.3*0.5; % m^2
x_st = 0.002; % m
x_cu = 0.003; % m
x = 0.002; % m
V_st = A*x_st; % m^3
V_cu = A*x_cu; % m^3

% Define constant temperatures
T_H = 400 + 273.15; % K
T_L = 100 + 273.15; % K 
T_0 = 400 + 273.15; % K
    
%% Example 2.1

% Define variables
R_st1 = x_st/(k_st*A); % K/W
R_cu1 = x_cu/(k_cu*A); % K/W
C_st1 = rho_st*V_st*c_st; % J/K
C_cu1 = rho_cu*V_cu*c_cu; % J/K

% Define impedeances
Z_st1 = 1/(s*C_st1);
Z_cu1 = 1/(s*C_cu1);

% Define terms used in mesh analysis
z_1 = R_st1/2 + Z_st1;
z_2 = -1*Z_st1;
z_3 = R_st1/2 + R_cu1/2 + Z_st1 + Z_cu1;
z_4 = -1*Z_cu1;

% Do mesh analysis
Z = [z_1, z_2, 0, 0; z_2, z_3, z_4  0; 0, z_4, z_3, z_2; 0, 0, z_2, z_1];
f = [(T_H - T_0)/s; 0; 0; (T_0 - T_L)/s];

% Find heat transfer rate
q_dot_1 = Z\f;

% Calculate temperature at nodes
T_11 = (T_H/s) - (R_st1/2)*q_dot_1(1);
T_21 = (T_11) - (R_st1 + R_cu1)*q_dot_1(2)/2;
T_31 = (T_21) - (R_st1 + R_cu1)*q_dot_1(3)/2;

% Convert to time domain
t_11 = impulse(T_11,time_1); % K
t_21 = impulse(T_21,time_1); % K
t_31 = impulse(T_31,time_1); % K
q_11 = impulse(q_dot_1(1),time_1)/1000; % kW
q_21 = impulse(q_dot_1(2),time_1)/1000; % kW
q_31 = impulse(q_dot_1(3),time_1)/1000; % kW
q_41 = impulse(q_dot_1(4),time_1)/1000; % kW

%% Example 2.1 Plots

set(gcf, 'Position', get(0, 'Screensize'));

% Temperature at Nodes
figure(1)
subplot(1,3,1)
plot(time_1,t_11,'LineWidth',1)
title('Temperature of Steel Block ($T_1$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([425,700]);

subplot(1,3,2)
plot(time_1,t_21,'LineWidth',1)
title('Temperature of Copper Block ($T_2$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([425,700]);

subplot(1,3,3)
plot(time_1,t_31,'LineWidth',1)
title('Temperature of Steel Block ($T_3$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([425,700]);

set(gcf, 'Position', get(0, 'Screensize'));

% Heat Transfer Rate at Boundaries
figure(2)
subplot(1,4,1)
plot(time_1,q_11,'LineWidth',1)
title('Heat Transfer Rate ($\dot{q_1}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([0,400]);

subplot(1,4,2)
plot(time_1,q_21,'LineWidth',1)
title('Heat Transfer Rate ($\dot{q_2}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([0,400]);

subplot(1,4,3)
plot(time_1,q_31,'LineWidth',1)
title('Heat Transfer Rate ($\dot{q_3}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([0,400]);

subplot(1,4,4)
plot(time_1,q_41,'LineWidth',1)
title('Heat Transfer Rate ($\dot{q_4}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,6]);
ylim([0,400]);

%% Example 2.2

% Define variables
R_st2 = x/(k_st*A); % K/W
R_cu2 = x/(k_cu*A); % K/W
C_st2 = rho_st*V_st*c_st; % J/K
C_cu2 = rho_cu*V_cu*c_cu; % J/K

% Define impedeances
Z_st2 = 1/(s*C_st2);
Z_cu2 = 1/(s*C_cu2);

% Define numerator and denominator used in nodal analysis
d_1 = (4/R_st2) + (1/Z_st2);
d_2 = (4/R_cu2) + (1/Z_cu2);
n_1 = ((2*T_H)/(s*R_st2)) + ((2*T_L)/(s*R_st2)) + T_0*C_st2;
n_2 = ((2*T_H)/(s*R_cu2)) + ((2*T_L)/(s*R_cu2)) + T_0*C_cu2;

% Calculate temperature at nodes
T_12 = n_1/d_1;
T_22 = n_2/d_2;
T_32 = n_1/d_1;

% Calculate heat transfer rates
Q_12_in = ((T_H/s) - T_12)/(R_st2/2);
Q_22_in = ((T_H/s) - T_22)/(R_cu2/2);
Q_32_in = ((T_H/s) - T_32)/(R_st2/2);
Q_12_out = ((T_12 - T_L/s))/(R_st2/2);
Q_22_out = ((T_22 - T_L/s))/(R_cu2/2);
Q_32_out = ((T_32 - T_L/s))/(R_st2/2);

% Convert to time domain
t_12 = impulse(T_12,time_2); % K
t_22 = impulse(T_22,time_2); % K
t_32 = impulse(T_32,time_2); % K
q_12_in = impulse(Q_12_in,time_2)/1000; % kW
q_22_in = impulse(Q_22_in,time_2)/1000; % kW
q_32_in = impulse(Q_32_in,time_2)/1000; % kW
q_12_out = impulse(Q_12_out,time_2)/1000; % kW
q_22_out = impulse(Q_22_out,time_2)/1000; % kW
q_32_out = impulse(Q_32_out,time_2)/1000; % kW

%% Example 2.2 Plots

set(gcf, 'Position', get(0, 'Screensize'));

% Temperature at Nodes
figure(3)
subplot(1,3,1)
plot(time_2,t_12,'LineWidth',1)
title('Temperature of Steel Block ($T_1$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,2]);
ylim([520,680]);

subplot(1,3,2)
plot(time_2,t_22,'LineWidth',1)
title('Temperature of Copper Block ($T_2$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,2]);
ylim([520,680]);

subplot(1,3,3)
plot(time_2,t_32,'LineWidth',1)
title('Temperature of Steel Block ($T_3$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
grid on
grid minor
xlim([0,2]);
ylim([520,680]);

set(gcf, 'Position', get(0, 'Screensize'));

% Heat Transfer Rates at Boundaries
figure(4)
subplot(3,2,1)
plot(time_2,q_12_in,'LineWidth',1)
title('Heat Transfer Rate into Steel ($\dot{q_{1,in}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,800]);

subplot(3,2,2)
plot(time_2,q_12_out,'LineWidth',1)
title('Heat Transfer Rate out of Steel ($\dot{q_{1,out}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,800]);

subplot(3,2,3)
plot(time_2,q_22_in,'LineWidth',1)
title('Heat Transfer Rate into Copper ($\dot{q_{2,in}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,20000]);

subplot(3,2,4)
plot(time_2,q_22_out,'LineWidth',1)
title('Heat Transfer Rate out of Copper ($\dot{q_{2,out}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,20000]);

subplot(3,2,5)
plot(time_2,q_32_in,'LineWidth',1)
title('Heat Transfer Rate into Steel ($\dot{q_{3,in}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,800]);

subplot(3,2,6)
plot(time_2,q_32_out,'LineWidth',1)
title('Heat Transfer Rate out of Steel ($\dot{q_{3,out}}$)','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer Rate (kW)','interpreter','latex')
grid on
grid minor
xlim([0,1.5]);
ylim([0,800]);
