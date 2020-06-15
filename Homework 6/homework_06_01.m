%% Homework 06 - ME342 - March 4, 2020

clear all; close all; clc;

%% Define parameters

s = tf('s');

% Boundary conditions
T_0 = 37 + 273.15; % K
T_inf = 25 + 273.15; % K

% Geometry 
delta_x = 0.5; % m
A = 0.1; % m^2

% Thermal properties
k = 14.4; % W/mK
c = 502.416; % J/kgK
rho = 8000; % kg/m^3
C = rho*delta_x*A*c; % J/K

% Stainless steel 304
% Conductivity: https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
% Specific Heat Capacitance: https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm
% Density: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=mq304a

% Define impedances
z_r = delta_x/(k*A); 
z_c = 1/(s*C);

% Define thermal diffusivity
alpha = k/(rho*c);

%% Part A Calculations

% Define terms used in impedance matrix
z3_11 = (z_r/6) + (3*z_c);
z3_12 = -3*z_c;
z3_22 = (z_r/3) + (6*z_c);

% Mesh analysis equations
z3 = [z3_11, z3_12, 0, 0; 
      z3_12, z3_22, z3_12, 0; 
      0, z3_12, z3_22, z3_12; 
      0, 0, z3_12, z3_11];
f3 = [(T_inf-T_0)/s; 0; 0; (T_0-T_inf)/s]; 
q_dot3 = z3\f3;

% Find temperature at nodes
T3_1 = (T_inf/s) - q_dot3(1)*(z_r/6);
T3_2 = T3_1 - q_dot3(2)*(z_r/3);
T3_3 = T3_2 - q_dot3(3)*(z_r/3);

% Convert to time domain
time = linspace(0,30000,10000);
y3_1 = impulse(T3_1,time);
y3_2 = impulse(T3_2,time);
y3_3 = impulse(T3_3,time);
t3_1 = y3_1 - 273.15;
t3_2 = y3_2 - 273.15;
t3_3 = y3_3 - 273.15;

%% Part A Plots

figure(1)
subplot(1,3,1)
plot(time,t3_1,'blue','LineWidth',1)
title('Temperature at Node 1','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,30000]);
ylim([24,38]);

subplot(1,3,2)
plot(time,t3_2,'blue','LineWidth',1)
title('Temperature at Node 2','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,30000]);
ylim([24,38]);

subplot(1,3,3)
plot(time,t3_3,'blue','LineWidth',1)
title('Temperature at Node 3','interpreter','latex','FontSize',14)
xlabel('Time (sec)','interpreter','latex','FontSize',12)
ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
grid on
grid minor
xlim([0,30000]);
ylim([24,38]);

%% Part B Calculations

% Define range of position and time to consider
t_range = linspace(0,30000,10000);
x_range = linspace(0,0.5,100);

% Create a matrix with the Cartesian product of time and position
[T_range,X_range] = meshgrid(t_range,x_range);
TX = [T_range(:),X_range(:)];

% Initialize temperature vector
U = zeros((length(t_range)*length(x_range)),1);

% Get "infinite" sum
for n = 1:100
    c_n = (2*(T_0-T_inf)*(1-((-1)^n)))/(n*pi);
    e_n = exp(((-1*alpha*n^2*pi^2)/(delta_x^2))*TX(:,1));
    s_n = sin(n*pi*TX(:,2)/delta_x);
    U = U + c_n.*e_n.*s_n;
end

% Convert to Celcius
U = U + T_inf - 273.15;

% Reshape vector into matrix
u = reshape(U,length(x_range),length(t_range));

%% Part C Plots

% Position of the three nodes
x_range3 = [0, 1/12, 1/4, 5/12, 1/2];

% Plot temperature using both solutions at various times
for index = [1,1001,2001,5001]
    apx3 = [25, t3_1(index), t3_2(index), t3_3(index), 25];
    figure(index + 1)
    plot(x_range3,apx3,'blue','LineWidth',1);
    hold on
    plot(x_range, u(:,index),'red','LineWidth',1)
    hold off
    title('Temperature of Exact and Approximate Solutions','interpreter','latex','FontSize',14)
    xlabel('Position (m)','interpreter','latex','FontSize',12)
    ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
    legend_animation = legend('Three Node Approximation','Exact Solution');
    set(legend_animation,'Interpreter','latex');
    grid on
    grid minor
    xlim([0,delta_x]);
    ylim([24,38]);
    pause(1);
end

%% Part D Calculations

% Number of nodes
n = 9;

% Define impedances that repeat
Z_R = z_r/(2*n);
Z_C = n*z_c;

% Define impedance matrix
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
F(1,1) = (T_inf-T_0)/s;
F(n+1,1) = (T_0-T_inf)/s;

% Find flow rates
Q_dot = Z\F;

% Define temperature of nodes 
T(1) = (T_inf/s) - (Q_dot(1)*Z_R);
for index = 2:n
    T(index) = T(index-1) - (2*Q_dot(index)*Z_R);
end

% Create matrix over time and distance
time_steps = 1000000;
time = linspace(0,30000,time_steps);   
t = zeros(time_steps,n+2);
t(:,1) = T_inf;
t(:,n+2) = T_inf;
for index = 1:n
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

%% Animation

% Plot temperature over distance 
for index = 1:100:10000
    figure(4)
    plot(pos,t(index*100,:)-273.15,'blue',x_range, u(:,index),'red','LineWidth',1);
    title('Temperature of Exact and Approximate Solutions','interpreter','latex','FontSize',14)
    xlabel('Position (m)','interpreter','latex','FontSize',12)
    ylabel('Temperature ($^{\mathrm o}$C)','interpreter','latex','FontSize',12)
    legend_animation = legend('Approximate Solutions','Exact Solution');
    set(legend_animation,'Interpreter','latex');
    grid on
    grid minor
    xlim([0,delta_x]);
    ylim([24,38]);
    video(((index-1)/100)+1) = getframe(gcf);
    drawnow
end

% Create video from plots
writerObj = VideoWriter('plots.avi');
writerObj.FrameRate = 10;

open(writerObj);
for index=1:length(video)
    % convert the image to a frame
    frame = video(index);    
    writeVideo(writerObj, frame);
end
close(writerObj);