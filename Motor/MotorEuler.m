clear; % close all;
X = -[0; 0]; ii = 0;
t_etapa = 1e-7;   % Paso de simulación
tF = 0.001;       % Tiempo final
u = 12;           % Entrada constante: voltaje al motor

for t = 0:t_etapa:tF
    ii = ii + 1;
    X = ModMotor(t_etapa, X, u); 
    x1(ii) = X(1);  % Velocidad angular (omega)
    x2(ii) = X(2);  % Derivada de velocidad (aceleración)
    acc(ii) = u;    % Entrada aplicada
end

% Vector de tiempo para graficar
t = 0:t_etapa:tF;

% Gráficas
subplot(2,1,1); hold on;
plot(t, x1, 'b'); title('Salida y, \omega_t (sin PID)');
ylabel('Velocidad [rad/s]');
grid on;

subplot(2,1,2); hold on;
plot(t, acc, 'b'); title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');
ylabel('Voltaje aplicado [V]');
grid on;


function [X]=ModMotor(t_etapa, xant, accion)
Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3;
Va=accion;
h=1e-7;
omega= xant(1);
wp= xant(2);
for ii=1:t_etapa/h
wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
wp=wp+h*wpp;
omega = omega + h*wp;
end
X=[omega,wp];
end