%% Wind Chill Effect - Results Table

clear all; close all; clc; 

%% Tabulate Cases

% Fluid Properties of Air
air_k = 0.0257; % W/mK
air_rho = 1.205; % kg/m^3
air_cp = 1009; % J/kgK
air_mu = 1.82*(10^-5); % kg/ms
air_beta = 3.41*(10^-3); % 1/K
air = [air_k,air_rho,air_cp,air_mu,air_beta];

% Fluid Properties of Water
water_k = 0.6; % W/mK
water_rho = 1000; % kg/m^3
water_cp = 4200; % J/kgK
water_mu = 10^-3; % kg/ms
water_beta = 2.07*(10^-4); % 1/K
water = [water_k,water_rho,water_cp,water_mu,water_beta];

% Diameter For Each Case
D = [0.02;0.02;0.02;0.02;0.08;0.08;0.02;0.02;0.08;0.08]; % m

% Velocity For Each Case
V = [0;0;0.447;4.47;0;4.47;0;1;0;1]; % m/s

% Skin Temperature
T_skin = (32+273.15)*ones(10,1); % K

% Ambient Temperature (only needed for natural convection)
T_inf = [0;-15;NaN;NaN;0;NaN;0;NaN;0;NaN] + 273.15; % K

% Forced/Natural Covection
f_n = [1;1;0;0;1;0;1;0;1;0];

% Table of All Values Needed
table = zeros(10);
table(:,1) = D; table(:,2) = V; table(:,3) = T_skin; table(:,4) = T_inf;
table(:,5:9) = [air;air;air;air;air;air;water;water;water;water];
table(:,10) = f_n;

%% Find Convection Coefficients Using Function

% Initialize Values
h = zeros(10,1);

% Apply Function to Each Case
for index = 1:10
    h(index) = find_convection(table(index,:));
end
