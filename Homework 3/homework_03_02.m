%% Rose Gebhardt and Harris Paspuleti -- February 12, 2020 -- Homework 3 Question 2

clear all; close all; clc;

%% Define parameters used

% Define system
s = tf('s');
m = 200; % kg
c = 0.45; % kJ/kgK
R = 2; % K/kW
C = m*c; % kJ/K
T_room = 25 + 273.15; % K
T_0 = 300 + 273.15; % K
omega = pi/60; % rad/s

% Vectors of varying time, resistance, and capacitance
time = linspace(0,1000,1000);

% Color
chromatic = [66, 135, 245; 69, 66, 245; 138, 66, 245; 176, 66, 245; 221, 66, 245;...
    245, 66, 227; 245, 66, 144; 245, 66, 93; 245, 75, 66; 245, 117, 66]/255;
set(groot,'defaultAxesColorOrder',chromatic)

%% 1. Plot System in Frequency Domain

% Define ambient temperature in frequency domain
T_inf = (25*omega)/(s^2 + omega^2) + (T_room/s);

% Define heat transfer rate in frequency domain
Q_dot = (T_inf - (T_0/s))/(R + (1/(s*C)));

% Define temperature of block in frequency domain
T = T_inf - R*Q_dot;

% Switch to time domain
t_inf = impulse(T_inf,time); % Used to compare phase
q_dot = impulse(Q_dot,time);
t = impulse(T,time);

% Plot temperature
figure(1)
plot(time,t,'LineWidth',1)
title('Temperature of Block','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (K)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,1000]);
ylim([250,600]);

% Plot heat transfer
figure(2)
plot(time,q_dot,'LineWidth',1)
title('Heat Transfer of System','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Heat Transfer (kW)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,1000]);
ylim([-140,20]);

%% 2. Vary resistance values

% Get temperature and heat transfer behavior for several resistances
t_R = zeros(length(time),10);
q_R = zeros(length(time),10);
for index = 1:10
    Q_it = (T_inf - (T_0/s))/((index*0.2) + (1/(s*C)));
    T_it = T_inf - index*0.2*Q_it;
    q_R(:,index) = impulse(Q_it,time);
    t_R(:,index) = impulse(T_it,time);
end

% Plot temperature for several resistances
figure(3)
plot(time,t_R,'LineWidth',1)
title('Temperature of Block','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature (K)','interpreter','latex','FontSize',12)
legendTR = legend('R = 0.2 K/kW','R = 0.4 K/kW','R = 0.6 K/kW','R = 0.8 K/kW','R = 1.0 K/kW'...
    ,'R = 1.2 K/kW','R = 1.4 K/kW','R = 1.6 K/kW','R = 1.8 K/kW','R = 2.0 K/kW');
set(legendTR,'Interpreter','latex');
grid on
grid minor
xlim([0,1000]);
ylim([250,600]);

% Plot heat transfer for several resistances
figure(4)
plot(time,q_R,'LineWidth',1)
title('Heat Transfer of System','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Heat Transfer (kW)','interpreter','latex','FontSize',12)
legendQR = legend('R = 0.2 K/kW','R = 0.4 K/kW','R = 0.6 K/kW','R = 0.8 K/kW','R = 1.0 K/kW'...
    ,'R = 1.2 K/kW','R = 1.4 K/kW','R = 1.6 K/kW','R = 1.8 K/kW','R = 2.0 K/kW','Location','southeast');
set(legendQR,'Interpreter','latex');
grid on
grid minor
xlim([0,1000]);
ylim([-400,150]);
