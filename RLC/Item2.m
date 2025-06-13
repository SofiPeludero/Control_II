clear all; clc; 

[z1]=xlsread('Curvas_Medidas_RLC_2024');
t0=z1(:,1);
I=z1(:,2);
Vc=z1(:,3);
Vin=z1(:,4);

%Visualizo los datos obtenido desde el archivo
figure(1);
subplot(3,1,1);hold on;
plot(t0,I,'r');title('Corriente Inductor I_L');grid;
subplot(3,1,2);hold on;
plot(t0,Vin,'r');title('Tension Entrada V_i');grid;
subplot(3,1,3);hold on;
plot(t0,Vc,'r');title('Tension Capacitor V_C');grid;
xlabel('Tiempo [Seg.]');

%El sistema RLC un sistema de segundo orden
%por lo tanto lo podemos aproximar FT=K/[(T1s+1)(T2s+1)]

%Comienza a aplicar Chen
StepAmplitude=12; %12 V de entrada en Va

tRet=       0.01;   %tiempo de retardo del sistema sacado del gráfico de la función Vc (salida)
tc = 0.001      ;      %tiempo de chen (el de ret - el elegido)

[val lugar] =min(abs(tc+tRet-t0)); %Busco en ret+t1
y_t=Vc(lugar);
t=t0(lugar)-tRet; %t1

[val lugar] =min(abs(2*tc+tRet-t0));
y_t2=Vc(lugar);
t2=t0(lugar)-tRet;

[val lugar] =min(abs(3*tc+tRet-t0));
y_t3=Vc(lugar);
t3=t0(lugar)-tRet;

% K=y(inf)/U   => en este caso y(inf)= 12 y U=Vin=12 
k       =       1;     %ganancia de estado estacionario

%CORROBORADOR DE PUNTOS extraidos
figure(2);
plot(t0,Vc);
hold on
plot(t+tRet,y_t,'o')
plot(t2+tRet,y_t2,'o')
plot(t3+tRet,y_t3,'o')
 
 %METODO DE CHEN sobre la tensión
k1      =       (1/StepAmplitude)*y_t/k-1;
k2      =       (1/StepAmplitude)*y_t2/k-1;
k3      =       (1/StepAmplitude)*y_t3/k-1;
b       =       4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1   =       (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2   =       (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta    =       (2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1      =       (-t/log(alfa1));
T2      =       (-t/log(alfa2));
T3      =       (beta*(T1-T2))+T1;

T1=real(T1);T2=real(T2); T3=real(T3);%importa sólo la parte real

sys_va=tf(k*[0 1],conv([T1 1],[T2 1])) %saco la funcion de transf
% como tiene un retardo el num debería multiplicarse por e^-Ls
%s=sym('s');
%G=sys_va*exp(-tRet*s)
[y1,t1,ent]=lsim(sys_va, Vin, t0, [0,0]);

%METODO DE CHEN sobre la corriente

tRet_I=       0.01;   %tiempo de retardo del sistema sacado del gráfico de la función Vc (salida)
tc_I = 0.001      ;      %tiempo de chen (el de ret - el elegido)

[val lugar] =min(abs(tc_I+tRet_I-t0)); %Busco en ret+t1
y_t_I=I(lugar);
t_I=t0(lugar)-tRet; %t1

[val lugar] =min(abs(2*tc_I+tRet_I-t0));
y_t2_I=I(lugar);
t2_I=t0(lugar)-tRet_I;

[val lugar] =min(abs(3*tc_I+tRet_I-t0));
y_t3_I=I(lugar);
t3_I=t0(lugar)-tRet_I;

% K=y(inf)/U   => en este caso y(inf)= 12 y U=Vin=12 
t_aux=t0(3:500);
I_aux=I(3:500);
y_Imax=max(I_aux);
umax=max(Vin);
y1i=y_Imax/umax;
%CORROBORADOR DE PUNTOS extraidos
figure(3);
plot(t0,I);
hold on;
plot(t_I+tRet_I,y_t_I,'o')
plot(t2_I+tRet_I,y_t2_I,'o')
plot(t3_I+tRet_I,y_t3_I,'o')

K=y1i*(T1-T2)/(alfa1-alfa2);

k1_I      =       (1/StepAmplitude)*y_t_I/K-1;
k2_I      =       (1/StepAmplitude)*y_t2_I/K-1;
k3_I      =       (1/StepAmplitude)*y_t3_I/K-1;
b_I       =       4*k1_I^3*k3_I-3*k1_I^2*k2_I^2-4*k2_I^3+k3_I^2+6*k1_I*k2_I*k3_I;
alfa1_I   =       (k1_I*k2_I+k3_I-sqrt(b))/(2*(k1_I^2+k2_I));
alfa2_I   =       (k1_I*k2_I+k3_I+sqrt(b))/(2*(k1_I^2+k2_I));
beta_I    =       (2*k1_I^3+3*k1_I*k2_I+k3_I-sqrt(b_I))/(sqrt(b_I));
T1_I      =       (-t_I/log(alfa1_I));
T2_I      =       (-t_I/log(alfa2_I));
T3_I      =       (beta*(T1-T2))+T1;

T1_I=real(T1_I);T2_I=real(T2_I); T3_I=real(T3_I);%importa sólo la parte real

sys_I=tf(K*[0 1],conv([T1_I 1],[T2_I 1])) %saco la funcion de transf
% como tiene un retardo el num debería multiplicarse por e^-Ls
%s=sym('s');
%G=sys_va*exp(-tRet*s)
[y2,t2,ent]=lsim(sys_I, Vin, t0, [0,0]);
 


 %OBTENCIÓN DE VALORES DE R,L Y C

 %ft de RLC con salida en Vc = 1/(LCs^2 + RCs + 1)
 %LC=T1*T2  y RC=T1+T2 -- De acuerdo con la función aproximada
%Si valoroL => C=(T1*T2)/L y R=(T1+T2)/C

C=K;
L=(T1*T2)/C;
R=(T1+T2)/C;

fprintf('\n---- PARÁMETROS IDENTIFICADOS DEL SISTEMA RLC ----\n');
fprintf('Resistencia R= %.4f ohmios\n', R);
fprintf('Capacitor C = %.4f Nm/A\n', Ki);
fprintf('Inductancia L = %.4f \n', Ra);


% figure (4)
% subplot(2,1,1);hold on; plot(t0,Vc,'b'); subplot(2,1,1); plot(t1,y1,'k');
% legend('Identificada','Real');
% subplot(2,1,2); hold on; plot(t0,I,'b'); subplot(2,1,2); plot(t2,y2,'k');
% legend('Identificada','Real');
