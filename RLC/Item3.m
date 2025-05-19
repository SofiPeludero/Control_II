clear all; clc; 

R=272.6843; L=0.1;C=9.8648e-6;

V0 = 12;       % Amplitud de la fuente de voltaje
f = 10;       % Frecuencia de la fuente (en Hz)
% Definir la función de voltaje 
V = @(t) V0 * square(2 * pi * f * t);  % Fuente de voltaje sinusoidal

% Definir el sistema de ecuaciones diferenciales
odefun = @(t, y) [y(2); (V(t) - R * y(2) - y(1) / C) / L];

% Resolver la ecuación diferencial
tspan = [0.05 0.2];  % Tiempo de simulación 
y0 = [0; 0];    % Condiciones iniciales (carga y corriente)
[t, y] = ode45(odefun, tspan, y0);

[z1]=xlsread('Curvas_Medidas_RLC_2024');
t0=z1(:,1);
I=z1(:,2);

hold on;
plot(t0,I,'r');title('Corriente Inductor I_L');grid;
plot(t, y(:, 2),'b');  % y(:, 2) es la corriente i(t)
xlim([0.049 0.2]);