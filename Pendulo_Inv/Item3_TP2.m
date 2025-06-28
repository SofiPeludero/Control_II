% Simulación de péndulo invertido con controlador y observador
clc; clear all;

% Parámetros físicos del sistema
m = 0.1;           % masa del péndulo [kg]
Fricc = 0.1;       % fricción del carro
long = 0.6;        % longitud del péndulo [m]
g = 9.8;           % gravedad [m/s^2]
M = 0.5;           % masa del carro [kg]
TamanioFuente = 12;

% Condición inicial del péndulo
alfa(1) = pi;      % equilibrio estable (hacia abajo)
colorc = '.r';
colorl = 'r';

% Modelo linealizado en el equilibrio estable
Mat_A = [0 1 0 0;
         0 -Fricc/M -m*g/M 0;
         0 0 0 1;
         0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_B = [0; 1/M; 0; 1/(long*M)];
xOP = [0; 0; pi; 0]; % punto de operación
Mat_C = [1 0 0 0;
         0 0 1 0];   % mide desplazamiento y ángulo

% Construcción del sistema ampliado (controlador con integrador)
Mat_Aa = [Mat_A zeros(4,1); -Mat_C(1,:) 0];
Mat_Ba = [Mat_B; 0];
Qa = diag([1e1 1e1 1e6 1e3 1e2]);
Ra = 1e4;

% Resolución del LQR ampliado
Ha = [Mat_Aa -Mat_Ba*inv(Ra)*Mat_Ba';
      -Qa    -Mat_Aa'];
[V, D] = eig(Ha);
MX1X2 = [];

for ii = 1:size(Ha,1)
    if real(D(ii,ii)) < 0
        MX1X2 = [MX1X2 V(:,ii)];
    end
end
MX1 = MX1X2(1:end/2,:);
MX2 = MX1X2(end/2+1:end,:);
Pa = real(MX2 / MX1);
Ka = inv(Ra) * Mat_Ba' * Pa;
K = Ka(1:4); KI = -Ka(5);

% Observador de estados (sistema dual)
Mat_Adual = Mat_A';
Mat_Bdual = Mat_C';
Mat_Cdual = Mat_B';
Qdual = diag([1e1 1 1e2 1]);
Rdual = diag([1e-2, 1e-1]);

Ha = [Mat_Adual -Mat_Bdual*inv(Rdual)*Mat_Bdual';
      -Qdual     -Mat_Adual'];
[V, D] = eig(Ha);
MX1X2 = [];
for ii = 1:size(Ha,1)
    if real(D(ii,ii)) < 0
        MX1X2 = [MX1X2 V(:,ii)];
    end
end
MX1 = MX1X2(1:end/2,:);
MX2 = MX1X2(end/2+1:end,:);
P = real(MX2 / MX1);
Ko = (inv(Rdual) * Mat_Bdual' * P)';

% Simulación
Ts = 1e-3;
T = 15;
KMAX = round(T / Ts);

% Condiciones iniciales
x = [0; 0; alfa(1); 0];
x_hat = [0; 0; 0; 0];
ref = 10;
h = Ts;
tita_pp = 0;
um = 0.1; % zona muerta
psi(1) = 0;

% Inicialización de variables
p(1) = x(1); p_p(1) = x(2); omega(1) = x(4);

for ki = 1:KMAX
    estado = [p(ki); p_p(ki); alfa(ki); omega(ki)];
    Y_ = Mat_C * estado;

    % Controlador con integrador
    psi_p = ref - Mat_C(1,:) * estado;
    psi(ki+1) = psi(ki) + psi_p * h;
    u(ki) = -K * (x_hat - xOP) + KI * psi(ki+1);

    % Zona muerta
    ui = u(ki);
    acci(ki) = ui;
    if abs(ui) < um
        u(ki) = 0;
    else
        u(ki) = ui - um * sign(ui);
    end

    % Modelo no lineal
    p_pp = (1 / (M + m)) * (u(ki) - m * long * tita_pp * cos(alfa(ki)) + m * long * omega(ki)^2 * sin(alfa(ki)) - Fricc * p_p(ki));
    tita_pp = (1 / long) * (g * sin(alfa(ki)) - p_pp * cos(alfa(ki)));

    % Integración
    p_p(ki+1) = p_p(ki) + h * p_pp;
    p(ki+1) = p(ki) + h * p_p(ki);
    omega(ki+1) = omega(ki) + h * tita_pp;
    alfa(ki+1) = alfa(ki) + h * omega(ki);

    % Observador
    x_hatp = Mat_A * (x_hat - xOP) + Mat_B * u(ki) + Ko * (Y_ - Mat_C * x_hat);
    x_hat = x_hat + h * x_hatp;
end

% Tiempo para graficar
t = 0:h:T;

% Gráficas
figure(1); hold on;
subplot(3,2,1); plot(t, alfa, '.m'); grid on; title('\phi_t','FontSize',TamanioFuente);
subplot(3,2,2); plot(t, omega, '.m'); grid on; title('$\dot{\phi_t}$','Interpreter','latex','FontSize',TamanioFuente);
subplot(3,2,3); plot(t, p, '.m'); grid on; title('\delta_t','FontSize',TamanioFuente);
subplot(3,2,4); plot(t, p_p, '.m'); grid on; title('$\dot{\delta_t}$','Interpreter','latex','FontSize',TamanioFuente);
subplot(3,1,3); plot(t(1:end-1), u, '.m', t(1:end-1), acci, '--k','LineWidth',1.3); grid on;
title('Acción de control','FontSize',TamanioFuente); xlabel('Tiempo en Seg.','FontSize',TamanioFuente);
legend('u','u\_inicial');

figure(2); hold on;
subplot(2,1,1); plot(alfa, omega, 'r'); grid on; xlabel('\phi_t','FontSize',TamanioFuente);
ylabel('$\dot{\phi_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);
subplot(2,1,2); plot(p, p_p, 'r'); grid on; xlabel('\delta_t','FontSize',TamanioFuente);
ylabel('$\dot{\delta_t}$','Interpreter','latex','Rotation',0,'FontSize',TamanioFuente);
