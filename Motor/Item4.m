clear all; clc;
% Parámetros
LAA = 366e-6;
J = 5e-9;
RA = 55.6;
B = 0;
Ki = 6.49e-3;
Km = 6.53e-3;
TL = 0;            % Torque de carga 
Va = 12;          % Voltaje de entrada

% Tiempo de simulación
dt = 1e-7;        % Paso
Tf = 0.05;            
N = round(Tf / dt);

% Variables de estado
ia = zeros(1, N);
wr = zeros(1, N);
theta = zeros(1, N);
t = linspace(0, Tf, N);

% Simulación Euler
for k = 1:N-1
    dia = (-RA / LAA) * ia(k) - (Km / LAA) * wr(k) + (1 / LAA) * Va;
    dwr = (Ki / J) * ia(k) - (B / J) * wr(k) - (1 / J) * TL;
    dtheta = wr(k);

    ia(k+1) = ia(k) + dt * dia;
    wr(k+1) = wr(k) + dt * dwr;
    theta(k+1) = theta(k) + dt * dtheta;
end

% Gráficas
figure(1);
subplot(2,1,1)
plot(t, wr, 'b');
xlabel('Tiempo [s]'); ylabel('\omega_r [rad/s]');
title('Velocidad angular del motor');
hold on;
grid on;
subplot(2,1,2)
plot(t, ia, 'r');  legend('Corriente de armadura')
xlabel('Tiempo [s]');ylabel('i_a [A]');
title('Corriente de armadura');
grid on;

iamax= max(ia)
TLmax=Ki*iamax

% Simulación Euler
for k = 1:N-1
    dwr = (Ki / J) * ia(k) - (B / J) * wr(k) - (1 / J) * TLmax;
    wr(k+1) = wr(k) + dt * dwr;
end

subplot(2,1,1)
plot(t, wr, 'g');legend('Velocidad angular sin Torque','Velocidad angular con Torque Maximo');
xlabel('Tiempo [s]'); ylabel('\omega_r [rad/s]');
title('Velocidad angular del motor con TLmax');
grid on;



