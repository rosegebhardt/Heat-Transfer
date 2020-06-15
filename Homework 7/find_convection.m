function h = find_convection(entries)
   
    % Define Each Entry
    D = entries(1); V = entries(2); T_skin = entries(3); T_inf = entries(4);
    k = entries(5); rho = entries(6); c_p = entries(7); mu = entries(8);
    beta = entries(9); f_n = entries(10);
    
    % Compute Non-Dimensional Parameters
    Re = (rho*V*D)/mu;
    Pr = (mu*c_p)/k;
    Gr = (9.81*beta*abs(T_skin-T_inf)*(rho^2)*(D^3))/(mu^2);
    Ra = Gr*Pr;
    
    % f_n == 0 --> forced convection
    % f_n == 1 --> natural convection
    
    % Forced Convection
    if f_n == 0
        
        % Possible Error Statements
        if Re < 0.4 || Re > 400000
            error('reynolds number is out of forced convection range!!! :(')
        end
        if Pr < 0.5 || Pr > 50
            error('prandtl number is out of forced convection range!!! :(')
        end
        
        % Coefficients for Each Case
        if Re >= 0.4 && Re < 4
            C = 0.989; n = 0.330;
        end        
        if Re >= 4 && Re < 35
            C = 0.911; n = 0.385;
        end
        if Re >= 35 && Re < 4083
            C = 0.683; n = 0.466;
        end
        if Re >= 4083 && Re < 40045
            C = 0.193; n = 0.618;
        end
        if Re >= 40045 && Re <= 400000
            C = 0.0266; n = 0.805;
        end
        
        % Compute Nusselt Number
        Nu = C*(Re^n)*(Pr^(1/3));
        
    end
    
    % Natural Convection
    if f_n == 1
        
        % Possible Error Statements
        if Ra < 10^4 || Ra > 10^12
            error('rayleigh number is out of natural convection range!!! :(')
        end
        
        % Coefficients for Each Case
        if Ra >= 10^4 && Ra < 2.12*(10^7)
            C = 0.53; n = 0.25;
        end
        if Ra >= 2.12*(10^7) && Ra <= 10^12
            C = 0.13; n = 0.3333;
        end
        
        % Compute Nusselt Number
        Nu = C*(Ra^n);
    end
    
    h = Nu*k/D;
    
end