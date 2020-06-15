%% Rose Gebhardt and Harris Paspuleti -- February 5, 2020 -- Homework 01

clear all; close all; clc;

%% Define parameters used

% Define system
m = 200; % kg
c = 0.45; % kJ/kgK
R = 1; % K/kW
C = m*c; % kJ/K
T_inf = 25 + 273.15; % K
T_0 = 300 + 273.15; % K

% Vectors of varying time, resistance, and capacitance
t = linspace(0,600,1000);

% Color
chromatic = [66, 135, 245; 69, 66, 245; 138, 66, 245; 176, 66, 245; 221, 66, 245;...
    245, 66, 227; 245, 66, 144; 245, 66, 93; 245, 75, 66; 245, 117, 66]/255;
set(groot,'defaultAxesColorOrder',chromatic)

%% 1. Plot System in Time Domain

% Define temperature as a function of time
T = T_inf - (T_inf-T_0)*exp((-1*t)/(R*C));

% Define heat transfer as a function of time
Q_dot = (1/R)*(T - T_inf);

% Plot temperature
figure(1)
plot(t,T,'LineWidth',1)
title('Temperature of Block','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
ylim([250,600]);

% Plot heat transfer
figure(2)
plot(t,Q_dot,'LineWidth',1)
title('Heat Transfer of System','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer (kW)','interpreter','latex')
ylim([-50,300]);

%% 2. Vary resistance values

% Get temperature and heat transfer behavior for several resistances
T_R = zeros(length(t),10);
Q_R = zeros(length(t),10);
for index = 1:10
    T_R(:,index) = T_inf - (T_inf-T_0)*exp((-1*t)/(index*0.2*C));
    Q_R(:,index) = (1/(index*0.2))*(T_R(:,index) - T_inf);
end

% Plot temperature for several resistances
figure(3)
plot(t,T_R,'LineWidth',1)
title('Temperature of Block','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
legendTR = legend('R = 0.2 K/kW','R = 0.4 K/kW','R = 0.6 K/kW','R = 0.8 K/kW','R = 1.0 K/kW'...
    ,'R = 1.2 K/kW','R = 1.4 K/kW','R = 1.6 K/kW','R = 1.8 K/kW','R = 2.0 K/kW');
set(legendTR,'Interpreter','latex');
ylim([250,600]);

% Plot heat transfer for several resistances
figure(4)
plot(t,Q_R,'LineWidth',1)
title('Heat Transfer of System','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer (kW)','interpreter','latex')
legendQR = legend('R = 0.2 K/kW','R = 0.4 K/kW','R = 0.6 K/kW','R = 0.8 K/kW','R = 1.0 K/kW'...
    ,'R = 1.2 K/kW','R = 1.4 K/kW','R = 1.6 K/kW','R = 1.8 K/kW','R = 2.0 K/kW');
set(legendQR,'Interpreter','latex');
ylim([-50,300]);

%% 3. Vary capacitance values

% Get temperature and heat transfer behavior for several capacitances
T_C = zeros(length(t),10);
Q_C = zeros(length(t),10);
for index = 1:10
    T_C(:,index) = T_inf - (T_inf-T_0)*exp((-1*t)/(R*(index*20)));
    Q_C(:,index) = (1/R)*(T_C(:,index) - T_inf);
end

% Plot temperature for several resistances
figure(5)
plot(t,T_C,'LineWidth',1)
title('Temperature of Block','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Temperature (K)','interpreter','latex')
legendTC = legend('c = 0.1 kJ/kgK','c = 0.2 kJ/kgK','c = 0.3 kJ/kgK','c = 0.4 kJ/kgK','c = 0.5 kJ/kgK'...
    ,'c = 0.6 kJ/kgK','c = 0.7 kJ/kgK','c = 0.8 kJ/kgK','c = 0.9 kJ/kgK','c = 1.0 kJ/kgK');
set(legendTC,'Interpreter','latex');
ylim([250,600]);

% Plot heat transfer for several resistances
figure(6)
plot(t,Q_C,'LineWidth',1)
title('Heat Transfer of System','interpreter','latex')
xlabel('Time (sec)','interpreter','latex')
ylabel('Heat Transfer (kW)','interpreter','latex')
legendQC = legend('c = 0.1 kJ/kgK','c = 0.2 kJ/kgK','c = 0.3 kJ/kgK','c = 0.4 kJ/kgK','c = 0.5 kJ/kgK'...
    ,'c = 0.6 kJ/kgK','c = 0.7 kJ/kgK','c = 0.8 kJ/kgK','c = 0.9 kJ/kgK','c = 1.0 kJ/kgK');
set(legendQC,'Interpreter','latex');
ylim([-50,300]);
