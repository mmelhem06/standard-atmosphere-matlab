clc 
clear all
close all ;
g0 = 9.80665;
T0 = 288.16;
P0 = 101325;
density0 = 1.225;
R = 287.05;
lambda = 1.4;
h_base = [0, 11000, 25000, 47000, 53000, 79000, 90000, 105000];
T_base = [0, 0, 0, 0, 0, 0, 0, 0];
P_base = [0, 0, 0, 0, 0, 0, 0, 0]; 
lapse_rate = [-0.0065, 0, 0.003, 0, -0.0045, 0, 0.004, 0];

T_base(1) = T0;
P_base(1) = P0;

for i = 2:8
    T_base(i) = T_base(i-1) + lapse_rate(i-1)*(h_base(i)-h_base(i-1));
end


for i = 2:8
    dh = h_base(i) - h_base(i-1);
    a = lapse_rate(i-1);

    if a ~= 0
        T_ratio = T_base(i) / T_base(i-1);
        P_base(i) = P_base(i-1) * T_ratio^(-(g0)/(a*R));
    else
        P_base(i) = P_base(i-1) * exp(-(g0)*dh/(R*T_base(i-1)));
    end
end

LayerData = [h_base; T_base; P_base; lapse_rate]; %matrix

h_input = input('Enter altitude in m: ');

if h_input < 11000
    k=LayerData(:,1);
elseif h_input >= 11000 && h_input < 25000
    k=LayerData(:,2);
elseif h_input >= 25000 && h_input < 47000
    k=LayerData(:,3);
elseif h_input >= 47000 && h_input < 53000
    k=LayerData(:,4);
elseif h_input >= 53000 && h_input < 79000
    k=LayerData(:,5);
elseif h_input >= 79000 && h_input < 90000
    k=LayerData(:,6);
elseif h_input >= 90000 && h_input < 105000
    k=LayerData(:,7);
else
    error('Altitude out of range')
end

hb = k(1);
tb = k(2);
pb = k(3);
a = k(4);

T_output = (tb + a*(h_input - hb));

dh = h_input - hb;
if a ~= 0
        T_ratio2 = T_output / tb;
        P_output = pb * T_ratio2^(-(g0)/(a*R));
else
        P_output = pb * exp(-(g0)*dh/(R*tb));
end

density_output = P_output/(R*T_output);

a_sound_output = sqrt(lambda*R*T_output);

disp(['Temperature:  ', num2str(T_output), ' K']);

disp(['Pressure:  ', num2str(P_output), ' Pa']);

disp(['Density:  ', num2str(density_output), ' kg/m^3']);

disp(['Speed of Sound:  ', num2str(a_sound_output), ' m/s']);

figure 
plot(T_base, h_base/1000, 'LineWidth', 2)
ylim([0 105])
grid on
xlabel('Temperature (K)')
ylabel('Altitude (km)')
title('Temperature vs Altitude (Standard Atmosphere)')

text(240, 2, 'Troposphere')
text(208, 18, 'Tropopause')
text(250, 32, 'Stratosphere1')
text(270, 50, 'Stratosphere2')
text(220, 62, 'Mesosphere')
text(170, 85, 'Mesopause')
text(205, 105,'Thermosphere')

text(240, 10, 'a1 = -6.5*10^{-3} K/km')
text(240, 40, 'a2 = +3*10^{-3} K/km')
text(230, 68, 'a3 = -4.5*10^{-3} K/km')
text(200,95, 'a4 = 4*10^{-3} K/km')


h_plot = 0:1000:100000;

n= length(h_plot);
t_plot = zeros (1,n);
p_plot = zeros (1,n);
density_plot = zeros (1,n);
a_plot = zeros (1,n);

for j= 1:n 
   h= h_plot(j);
if h < 11000
    k=LayerData(:,1);
elseif h >= 11000 && h < 25000
    k=LayerData(:,2);
elseif h >= 25000 && h < 47000
    k=LayerData(:,3);
elseif h >= 47000 && h < 53000
    k=LayerData(:,4);
elseif h >= 53000 && h < 79000
    k=LayerData(:,5);
elseif h >= 79000 && h < 90000
    k=LayerData(:,6);
elseif h >= 90000 && h < 105000
    k=LayerData(:,7);
else
    error('Altitude out of range')
end

hb = k(1);
tb = k(2);
pb = k(3);
a = k(4);
dh = h - hb;

T = (tb + a*(h - hb));

if a ~= 0
        T_ratio2 = T / tb;
        P = pb * T_ratio2^(-(g0)/(a*R));
else
        P = pb * exp(-(g0)*dh/(R*tb));
end

t_plot(j)= T;
p_plot(j) = P;
density_plot(j)= P/(R*T);
a_plot(j) = sqrt(lambda*R*T);
end

AltitudeInKm = h_plot'/1000;

Standard_Atmospheric_Table= table (AltitudeInKm,t_plot',p_plot',density_plot',a_plot', ...
    'VariableNames',{'AltitudeInKm','Temperature in K','Pressure In P', ...
    'Density in Kg/m3','SpeedOfSound In m/s'});
disp(Standard_Atmospheric_Table)
writetable(Standard_Atmospheric_Table,'Standard_Atmospheric_Table.xlsx')


figure 
plot(t_plot, AltitudeInKm, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Temperature (K)');
ylabel('Altitude (km)');
title('Temperature vs Altitude');
hold on
plot(T_output, h_input/1000, 'x', 'MarkerSize', 8, 'LineWidth', 2)
hold off

figure
semilogx(p_plot, AltitudeInKm, 'y-', 'LineWidth', 1.5);
grid on;
xlabel('Pressure (Pa)');
ylabel('Altitude (km)');
title('Pressure vs Altitude');
hold on
plot(P_output, h_input/1000, 'x', 'MarkerSize', 8, 'LineWidth', 2)
hold off

figure
semilogx(density_plot, AltitudeInKm, 'g-', 'LineWidth', 1.5);
grid on;
xlabel('Density (kg/m^3)');
ylabel('Altitude (km)');
title('Density vs Altitude');
hold on
plot(density_output, h_input/1000, 'x', 'MarkerSize', 8, 'LineWidth', 2)
hold off

figure
plot(a_plot, AltitudeInKm, 'm-', 'LineWidth', 1.5);
grid on;
xlabel('Speed of Sound (m/s)');
ylabel('Altitude (km)');
title('Speed of Sound vs Altitude');
hold on
plot(a_sound_output, h_input/1000, 'x', 'MarkerSize', 8, 'LineWidth', 2)
hold off
