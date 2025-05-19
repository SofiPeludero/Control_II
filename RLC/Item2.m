clear all; clc; 

[z1]=xlsread('Curvas_Medidas_RLC_2024');
t0=z1(:,1);
I=z1(:,2);
Vc=z1(:,3);
Vin=z1(:,4);

%Visualizo los datos obtenido desde el archivo
figure
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
% figure;
 %plot(t0,Vc);
 %hold on
 %plot(t+tRet,y_t,'o')
 %plot(t2+tRet,y_t2,'o')
 %plot(t3+tRet,y_t3,'o')
 
 %METODO DE CHEN
k1      =       (1/StepAmplitude)*y_t/k-1;
k2      =       (1/StepAmplitude)*y_t2/k-1;
k3      =       (1/StepAmplitude)*y_t3/k-1;
b       =       4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1   =       (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2   =       (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta    =       (2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1      =       (-t/log(alfa1))
T2      =       (-t/log(alfa2))
T1=real(T1);T2=real(T2);%importa sólo la parte real

sys_va=tf(k,conv([T1 1],[T2 1])) %saco la funcion de transf
% como tiene un retardo el num debería multiplicarse por e^-Ls
%s=sym('s');
%G=sys_va*exp(-tRet*s)

[y1,t1,ent]=lsim(sys_va, Vin, t0, [0,0]);

 figure
 plot(t0,Vc);
 hold on
 plot(t1,y1,'k')
 plot([t+tRet t2+tRet t3+tRet],[y_t y_t2 y_t3],'o'); legend('DatosExcel','Modelo','Mediciones');
 
 %OBTENCIÓN DE VALORES DE R,L Y C
 
 %ft de RLC con salida en Vc = 1/(LCs^2 + RCs + 1)
 %LC=T1*T2  y RC=T1+T2 -- De acuerdo con la función aproximada
%Si valoro L => C=(T1*T2)/L y R=(T1+T2)/C
L=0.1;
C=(T1*T2)/L
R=(T1+T2)/C

figure (4)
step(StepAmplitude*sys_va,'r',0.1),hold on
plot(t0,Vc,'b');
legend('Identificada','Real');