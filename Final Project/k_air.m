function k_air = k_air(T)
   
% This function is for computing the conduction coefficient of ideal air at
% 1 atm with varying temperatures. This is to be used to find conduction
% effects iteratively as temperature changes.

% This function take input temperatures in Celcuis!!!

% Interpolated Values of k
% Source for values can be found in the link below:
% https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm

if T<0
    error('the fabric is literally on fire, why is the temperture less than 0 degrees C?')
end
if T>400
    error('the firefighter is burnt and the material is useless')
end
if T>=0 && T<=5
    T_L = 0; T_H = 5; K_L = 0.02364; K_H = 0.02401;
end
if T>5 && T<=10
    T_L = 5; T_H = 10; K_L = 0.02401; K_H = 0.02439;
end
if T>10 && T<=15
    T_L = 10; T_H = 15; K_L = 0.02439; K_H = 0.02476;
end
if T>15 && T<=20
    T_L = 15; T_H = 20; K_L = 0.02476; K_H = 0.02514;
end
if T>20 && T<=25
    T_L = 20; T_H = 25; K_L = 0.02514; K_H = 0.02551;
end
if T>25 && T<=30
    T_L = 25; T_H = 30; K_L = 0.02551; K_H = 0.02588;
end
if T>30 && T<=35
    T_L = 30; T_H = 35; K_L = 0.02588; K_H = 0.02625;
end
if T>35 && T<=40
    T_L = 35; T_H = 40; K_L = 0.02625; K_H = 0.02662;
end
if T>40 && T<=45
    T_L = 40; T_H = 45; K_L = 0.02662; K_H = 0.02699;
end
if T>45 && T<=50
    T_L = 45; T_H = 50; K_L = 0.02699; K_H = 0.02735;
end
if T>50 && T<=60
    T_L = 50; T_H = 60; K_L = 0.02735; K_H = 0.02808;
end
if T>60 && T<=70
    T_L = 60; T_H = 70; K_L = 0.02808; K_H = 0.02881;
end
if T>70 && T<=80
    T_L = 70; T_H = 80; K_L = 0.02881; K_H = 0.02953;
end
if T>80 && T<=90
    T_L = 80; T_H = 90; K_L = 0.02953; K_H = 0.03024;
end
if T>90 && T<=100
    T_L = 90; T_H = 100; K_L = 0.03024; K_H = 0.03095;
end
if T>100 && T<=120
    T_L = 100; T_H = 120; K_L = 0.03095; K_H = 0.03235;
end
if T>120 && T<=140
    T_L = 120; T_H = 140; K_L = 0.03235; K_H = 0.03374;
end
if T>140 && T<=160
    T_L = 140; T_H = 160; K_L = 0.03374; K_H = 0.03511;
end
if T>160 && T<=180 
    T_L = 160; T_H = 180; K_L = 0.03511; K_H = 0.03646;
end
if T>180 && T<=200
    T_L = 180; T_H = 200; K_L = 0.03646; K_H = 0.03779;
end
if T>200 && T<=250
    T_L = 200; T_H = 250; K_L = 0.03779; K_H = 0.04104;
end
if T>250 && T<=300
    T_L = 250; T_H = 300; K_L = 0.04104; K_H = 0.04418;
end
if T>300 && T<=350
    T_L = 300; T_H = 350; K_L = 0.04418; K_H = 0.04721;
end
if T>350 && T<=400
    T_L = 350; T_H = 400; K_L = 0.04721; K_H = 0.05015;
end

k_air = ((K_H-K_L)/(T_H-T_L))*(T-T_L) + K_L;
    
end