%% Homework 06 - ME342 - March 4, 2020

clear all; close all; clc;

%% Define parameters

s = tf('s');

% Boundary conditions
T_0 = 37 + 273.15; % K
T_inf = 25 + 273.15; % K

% Geometry 
delta_x = 0.002; % m
A = 0.01; % m^2

% Thermal properties
k = 14.4; % W/mK
c = 502.416; % J/kgK
rho = 8000; % kg/m^3
C = rho*delta_x*A*c; % J/K

% Stainless steel 304
% Conductivity: https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
% Specific Heat Capacitance: https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm
% Density: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mq304a

% Define impedeances
z_r = delta_x/(k*A); 
z_c = 1/(s*C);

%% Calculations

% Define forcing temperature in s domain
T_osc = ((650+273.15)/s) + ((300+273.15)*s)/(s^2 + 256*pi^2);

% Number of nodes
n = 8;

% Define impedeances that repeat
Z_R = z_r/(2*n);
Z_C = n*z_c;

% Define impedeance matrix
Z(1,1) = Z_R + Z_C; 
Z(1,2) = -1*Z_C;
for index = 2:n
    Z(index,index-1) = -1*Z_C;
    Z(index,index) = 2*Z_R + 2*Z_C;
    Z(index,index+1) = -1*Z_C;
end
Z(n+1,n) = -1*Z_C; 
Z(n+1,n+1) = Z_R + Z_C;

% Define forcing vector
F(1,1) = T_osc- (T_0/s);
F(n+1,1) = (T_0-T_inf)/s;

% Find flow rates
Q_dot = Z\F;

% Define temperature of nodes 
T(1) = (T_osc) - (Q_dot(1)*Z_R);
for index = 2:n+1
    T(index) = T(index-1) - (2*Q_dot(index)*Z_R);
end

% Create matrix over time and distance
time_steps = 10000;
time = linspace(0,1,time_steps);   
t = zeros(time_steps,n+2);
t(:,1) = impulse(T_osc,time);
for index = 1:(n+1)
    t(:,index+1) = impulse(T(index),time);
end

% Define position vector
pos_frac = zeros(n+2,1);
pos_frac(1) = 0;
pos_frac(n+2) = 1;
for index = 1:n
    pos_frac(index+1) = (2*index - 1)/(2*n);
end
pos = pos_frac*delta_x;

%%  Find Amplitude 

% Initialize amplitude vector
A = zeros(length(pos),1);

% Find amplitude at steady-state
for index = 1:length(pos)
    A(index) = (max(t(5000:6000,index)) - min(t(5000:6000,index)))/2;
end

figure(1)
plot(pos*100,A,'blue','LineWidth',1)
title('Amplitude as a Function of Depth','interpreter','latex','FontSize',14)
xlabel('Depth (mm)','interpreter','latex','FontSize',12)
ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,0.2]);
ylim([0,600]);
