clc; clear all;

[z1]=xlsread('Curvas_Medidas_Motor_2024_v');
t0=z1(:,1);
w=z1(:,2);
Ia=z1(:,3);
Vi=z1(:,4);
Tl=z1(:,5);
StepAmplitude = 12;

%Visualizo los datos obtenido desde el archivo
figure (1);
subplot(2,2,1); hold on;
plot(t0,w,'r'); title('Velocidad angular w[rad/seg]'); grid;
subplot(2,2,2); plot(t0,Ia,'r');
title('Corriente de armadura Ia'); grid;
subplot(2,2,3); plot(t0,Vi,'r');title('Tension V'); grid;
xlabel('Tiempo [Seg.]');
subplot(2,2,4); plot(t0,Tl,'r');title('Torque'); grid;
xlabel('Tiempo [Seg.]');
hold off;

%METODO DE CHEN PARA w (velocidad angular)
tRet = 0.03502;
tc = 0.03507 - tRet;   %tiempo de chen (el elegido - el t retardo)
[val lugar] = min(abs(tc + tRet - t0));
y_t = w(lugar);
t = t0(lugar) - tRet;

[val lugar] = min(abs(2*t + tRet - t0));
y_t2 = w(lugar);
t2 = t0(lugar) - tRet;

[val lugar] = min(abs(3*t + tRet - t0));
y_t3 = w(lugar);
t3 = t0(lugar) - tRet;

k = 198.2 / StepAmplitude;

%CORROBORADOR DE PUNTOS extraidos
% figure(2);
% plot(t0,w);
% hold on;
% plot(t+tRet,y_t,'o');
% plot(t2+tRet,abs(y_t2),'o');
% plot(t3+tRet,y_t3,'o');

k1 = (1/StepAmplitude) * y_t / k - 1;
k2 = (1/StepAmplitude) * y_t2 / k - 1;
k3 = (1/StepAmplitude) * y_t3 / k - 1;

b = 4*(k1^3)*k3 - 3*(k1^2)*(k2^2) - 4*(k2^3) + (k3^2) + 6*k1*k2*k3;
alfa1 = (k1*k2 + k3 - sqrt(b)) / (2*((k1^2) + k2));
alfa2 = (k1*k2 + k3 + sqrt(b)) / (2*((k1^2) + k2));
beta = (2*(k1^3) + 3*k1*k2 + k3 - sqrt(b)) / sqrt(b);

T1 = real(-t / log(alfa1));
T2 = real(-t / log(alfa2));
T3 = real(beta * (T1 - T2) + T1);

sys_va = tf(k, conv([T1 1], [T2 1]));

%%%%%%%%%%%%%%%%%%
dt = 1e-5;
t_s = 0:dt:t0(end-1);

u1_Va = zeros(tRet/dt, 1);
u2_Va = 12 * ones(fix((0.6 - tRet)/dt), 1); 
u_Va = [u1_Va; u2_Va];

u1_T = zeros(fix(0.185/dt) + 1, 1);
u2_T = ones(fix((0.6 - 0.100)/dt), 1);
u_T = [u1_T; u2_T];

[y1, t1, ent] = lsim(sys_va, u_Va, t_s, [0 0]);

%METODO CHEN PARA Tl 
ret_tl = 0.18500;
tc_tl = 0.18505-ret_tl;

[val lugar] = min(abs(tc_tl + ret_tl - t0));
y_t_tl = w(lugar);
t_tl = t0(lugar) - ret_tl;

[val lugar] = min(abs(2*tc_tl + ret_tl - t0));
y_t2_tl = w(lugar);
t2_tl = t0(lugar) - ret_tl;

[val lugar] = min(abs(3*tc_tl + ret_tl - t0));
y_t3_tl = w(lugar);
t3_tl = t0(lugar) - ret_tl;

w1=198.2;                                  %Velocidad con TL=0
w2=164;                                %Velocidad angular con TL
TL=max(Tl)*0.93;
%k_tl=(w1-w2)/Tlmax;           %ganancia para Tl

%TL = 0.0011;
k_tl = -(w2 - w1) / TL;

yid_1 = -(y_t_tl - 198.2);
yid_2 = -(y_t2_tl - 198.2);
yid_3 = -(y_t3_tl - 198.2);

%CORROBORADOR DE PUNTOS extraidos
% figure(3);
% plot(t_tl+ret_tl,y_t_tl,'o');
% plot(t2_tl+ret_tl,y_t2_tl,'o');
% plot(t3_tl+ret_tl,y_t3_tl,'o');
% hold off;

k1_tl = (1/TL) * yid_1 / k_tl - 1;
k2_tl = (1/TL) * yid_2 / k_tl - 1;
k3_tl = (1/TL) * yid_3 / k_tl - 1;

b_tl = 4*(k1_tl^3)*k3_tl - 3*(k1_tl^2)*(k2_tl^2) - 4*(k2_tl^3 )+ (k3_tl^2) + 6*k1_tl*k2_tl*k3_tl;
alfa1_tl = (k1_tl*k2_tl + k3_tl - sqrt(b_tl)) / (2*((k1_tl^2) + k2_tl));
alfa2_tl = (k1_tl*k2_tl + k3_tl + sqrt(b_tl)) / (2*((k1_tl^2) + k2_tl));
beta_tl = (2*(k1_tl^3) + 3*k1_tl*k2_tl + k3_tl - sqrt(b_tl)) / sqrt(b_tl);

T1_tl = real(-t_tl / log(alfa1_tl));
T2_tl = real(-t_tl / log(alfa2_tl));
T3_tl =(real(beta_tl * (T1_tl - T2_tl) + T1_tl));

sys_T = tf(k_tl * [T3_tl 1], conv([T1_tl 1], [T2_tl 1]));

u2_T = TL * ones(fix((0.6 - 0.100)/dt), 1);
u_T = [u1_T; u2_T];

% Ajuste para evitar errores de longitud
N = length(t_s);
u_T = u_T(1:min(length(u_T), N));
u_T(N) = u_T(end);  

[y2, t2_, ent] = lsim(sys_T, u_T, t_s, [0 0]);

figure(4);
hold on;
plot(t0,w);
plot(t_s, y1 - y2, 'r');
legend('Datos', 'Modelado');
plot(t + tRet, y_t, 'o', t2 + tRet, y_t2, 'o', t3 + tRet, y_t3, 'o');
plot(t_tl + ret_tl, y_t_tl, 'o', t2_tl + ret_tl, y_t2_tl, 'o', t3_tl + ret_tl, y_t3_tl, 'o');

%Metodo de Chen para I, para obtener I/Vin
tRet_i = 0.0350;             % Tiempo de retardo estimado
I_aux = Ia(31:50860);    % Parte útil de la respuesta (corriente)
t_aux = t0(31:50860);    % Tiempo correspondiente
u_aux = Vi(31:50860);      % Entrada correspondiente

y1i = max(I_aux);         
uMax = max(Vi);            % Tensión máxima de entrada

% [val lugar] = min(abs(tc_tl + ret_tl - t0));
% y_t_tl = w(lugar);
% t_tl = t0(lugar) - ret_tl;
% 
% [val lugar] = min(abs(2*tc_tl + ret_tl - t0));
% y_t2_tl = w(lugar);
% t2_tl = t0(lugar) - ret_tl;
% 
% [val lugar] = min(abs(3*tc_tl + ret_tl - t0));
% y_t3_tl = w(lugar);
% t3_tl = t0(lugar) - ret_tl;
% Determinar los tiempos característicos
[val_i lugar] = min(abs(y1i - I_aux));
t1i = t_aux(lugar);
y1i = y1i / uMax;

[val_i lugar] = min(abs((t1i - tRet_i)*2 + tRet_i - t0));
t2i = t0(lugar);
y2i = Ia(lugar) / uMax;

[val_i lugar] = min(abs((t1i - tRet_i)*3 + tRet_i - t0));
t3i = t0(lugar);
y3i = Ia(lugar) / uMax;

alfa1_i = real((y2i - 4*y1i*y3i*sqrt(1/(4*y1i*y3i - 3*y2i^2)) + 3*(y2i^2)*sqrt(1/(4*y1i*y3i - 3*(y2i^2)))) / (2*y1i));
alfa2_i =real(y2i/y1i - alfa1);
T1_i = real(-(t1i - tRet_i) / log(alfa1_i));
T2_i = real(-(t1i - tRet_i) / log(alfa2_i));
beta_i = y1i / (alfa1_i - alfa2_i);
K_i = (y1i * (T1_i - T2_i) / (alfa1_i - alfa2_i));

sys_I=K_i*tf([1 0],conv([T1_i 1],[T2_i 1]))

[y_i,t_i,ent]=lsim(sys_I, Vi, t0, [0,0]);

T3_i=T3_tl;

% IDENTIFICACIÓN DE PARÁMETROS DEL MOTOR DC
La=(T1_i*T2_i)/K_i;
Ra=(T1_i + T2_i)/K_i;
Km=(T1_i*T2_i + (T2_i^2)-(T3_i*(T1_i+T2_i)))/(k*(T3_i^2));
%Km=1/k;
Ki=k*Ra/k_tl;
Jm=(T1_i + T2_i)/k_tl;
Bm=0;


% === Resultados ===
fprintf('\n---- PARÁMETROS IDENTIFICADOS DEL MOTOR DC ----\n');
fprintf('Constante de torque Km = %.4f Nm/A\n', Km);
fprintf('Constante del par motor Ki = %.4f Nm/A\n', Ki);
fprintf('Resistencia de armadura Ra = %.4f ohmios\n', Ra);
fprintf('Inductancia de armadura La = %.6f H\n', La);
fprintf('Momento de inercia Jm = %.6e kg·m²\n', Jm);
fprintf('Coef. de fricción viscosa Bm = %.6e N·m·s/rad\n', Bm);

figure(5);
plot(t0, Ia, 'b', 'DisplayName', 'Datos originales');
hold on;
plot(t_i, y_i, 'r--', 'DisplayName', 'Modelo estimado');
%plot([t1i t2i t3i], [y1i y2i y3i]*uMax, 'ko', 'DisplayName', 'Puntos medidos');
xlabel('Tiempo [s]'); ylabel('Corriente [A]');
%xlim([0 0.185]);
title('Sistema estimado (segundo orden con un cero)');
legend;
grid on;