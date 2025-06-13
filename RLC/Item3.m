clear all; clc; 

%Valores calculados en el item anterior: 
R=267.5231; L= 0.0981;C= 1.0055e-05;


V0 = 12;       % Amplitud de la fuente de voltaje
f = 10;       % Frecuencia de la fuente (en Hz)
% Definir la función de voltaje 
V = @(t) V0 * square(2 * pi * f * t);  % Fuente de voltaje sinusoidal

% Definir el sistema de ecuaciones diferenciales
odefun = @(t, y) [y(2); (V(t) - R * y(2) - y(1) / C) / L];

% Resolver la ecuación diferencial
tspan = [0.01 0.2];  % Tiempo de simulación 
y0 = [0; 0];    % Condiciones iniciales (carga y corriente)
[t, y] = ode45(odefun, tspan, y0);

[z1]=xlsread('Curvas_Medidas_RLC_2024');
t0=z1(:,1);
I=z1(:,2);

hold on;
plot(t0,I,'r');title('Corriente Inductor I_L');grid;
plot(t, y(:, 2),'b');  % y(:, 2) es la corriente i(t)
xlim([0.0 0.2]);
legend('Real','Modelo');
