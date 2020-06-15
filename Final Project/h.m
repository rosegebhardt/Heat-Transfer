function h = h(T1,T2)
   
% This function is for computing the convection coefficient of a flat plate
% assuming the convective fluid is air at 1 atm and 25 degrees C. Ideal
% gas assumption is used.

% Fluid Properties
% Source for values can be found in the link below:
% https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
rho = 1.184;            % kg/m^3
c_p = 1007;             % J/kgK
k = 0.02551;            % W/mK
mu = 1.849*(10^-5);     % kg/ms
g = 9.81;               % m/s^2
L = 0.255;              % m
nu = mu/rho;            % m^2/s
alpha = k/(rho*c_p);    % m^2/s

% Dimensionless Parameters
% Assumes beta = 1/T, where T is the average of the two temperatures
Pr = nu/alpha;
Gr = (2*g*abs(T2-T1)*(L^3))/((T1+T2)*(nu^2));
Ra = Gr*Pr;

% Compute Nusselt Number: Note: differs from paper's computation, we use
% power law found in Sidebo's textbook (Table 9.2)

% Sidebotham's Computation
% Laminar Flow
if Ra < (4.545*10^9)
    Nu = 0.59*(Ra^0.25);
end
% Turbulent Flow
if Ra >= (4.545*10^9)
    Nu = 0.021*(Ra^0.4);
end

% % Churchill and Chu's Computation
% % Laminar Flow
% if Ra < (10^9)
%     Nu = 0.68 + ((0.67*(Ra^(1/4)))/((1+((0.492/Pr)^(9/16)))^(4/9)));
% end
% % Turbulent Flow
% if Ra >= (10^9)
%     Nu = (0.825 + ((0.387*(Ra^(1/4)))/((1+((0.492/Pr)^(9/16)))^(8/27))))^2;
% end

h = Nu*k/L;
    
end