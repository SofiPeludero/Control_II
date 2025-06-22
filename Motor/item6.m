clear all; clc;
%Parametros del motor item 5
Km = 1.0831e-02;   % Nm/A
Ki_m = 9.3023e-03;   % V*s/rad
Ra = 18.3275;      % Ohm
Laa = 4.0207e-04;   % H
Jm = 2.7398e-09;   % kg*m^2
B = 0;             % Fricción despreciable
TL_max = 0.0011;   % Torque perturbador máximo

% TIEMPO DE SIMULACIÓN
tau_el = Laa / Ra;                % Tiempo eléctrico (corriente)
tau_mec = Jm / (Km * Ki_m / Ra);    % Tiempo mecánico aproximado
Ts=min(tau_el, tau_mec) / 10;

%Ts = 1e-4;         % Paso de simulación
Tfinal = 0.1;      % Tiempo total
N = round(Tfinal / Ts);

% angulo de referencia
ref = 57.2958; %1pi[rad]=180 => 1rad=57.2958
%ref=pi/2;

Kp = 0.001;
Ki= 0.5;
Kd =0.11;

A1 = ((2*Kp*Ts)+(Ki*Ts^2)+(2*Kd))/(2*Ts);
B1 = ((-2*Kp*Ts)+(Ki*Ts^2)-(4*Kd))/(2*Ts);
C1 = Kd/Ts;

%Variables de est
X = [0; 0; 0]; % [ia; wr; theta]
%e = zeros(N+3, 1); 
u = 0; psi = 0;
e = zeros(N+2,1);

x1 = zeros(N,1); 
x2 = zeros(N,1); 
x3 = zeros(N,1);
acc = zeros(N,1); 
TL_hist = zeros(N,1); 
ref_hist = zeros(N,1);

%Del torque
TL_ap = zeros(N,1);
t = (0:N-1)*Ts;
TL_ap((t>=0.2)&(t<=0.33)) = TL_max;
TL_ap((t>=0.5)&(t<=0.63)) = TL_max;

% Controlador
for k = 3:N
    % Error con referencia
    e(k) = ref - X(3); % ángulo
    
    % PID discreto
    u = u + A1*e(k) + B1*e(k-1) + C1*e(k-2);

    % Simular motor con entrada 'u' y perturbación TL_ap(k)
    X = ModMotor(Ts, X, [u, TL_ap(k)]);

    % Guardar resultados
    x1(k) = X(1); % Corriente
    x2(k) = X(2); % Velocidad angular
    x3(k) = X(3); % Ángulo
    acc(k) = u;
    TL_hist(k) = TL_ap(k);
    ref_hist(k) = ref;
end

figure(1);
subplot(3,1,1);
plot(t, ref_hist, 'k--', t, x3, 'b');
title('Salida y, \theta '); ylabel('\theta (rad)'); grid on; legend('Ref','Salida');
%xlim([-0.000001 0.001]);

subplot(3,1,2);
plot(t, x1, 'r'); title('Corriente i_a'); ylabel('i_a (A)'); grid on;

subplot(3,1,3);
plot(t, acc, 'm'); title('Acción de Control V_a'); ylabel('V_a (V)'); xlabel('Tiempo (s)'); grid on;
