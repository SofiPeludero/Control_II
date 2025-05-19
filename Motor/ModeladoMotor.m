clear all; clc; 

[z1]=xlsread('Curvas_Medidas_Motor_2024_v');
t0=z1(:,1);
w=z1(:,2);
Ia=z1(:,3);
Vi=z1(:,4);
Tl=z1(:,5);

%Visualizo los datos obtenido desde el archivo
figure
subplot(2,2,1);hold on;
plot(t0,w,'r');title('Velocidad angular w[rad/seg]');grid;
subplot(2,2,2);hold on;
plot(t0,Ia,'r');title('Corriente de armadura Ia');grid;
subplot(2,2,3);hold on;
plot(t0,Vi,'r');title('Tension V');grid;
xlabel('Tiempo [Seg.]');
subplot(2,2,4);hold on;
plot(t0,Tl,'r');title('Torque');grid;
xlabel('Tiempo [Seg.]');

%Comienza a aplicar Chen
StepAmplitude=12; %12 V de entrada en Va

tRet=0.03509;   %tiempo de retardo del sistema sacado del gráfico de la función Vc (salida)
tc =0.03511-tRet;      %tiempo de chen (el elegido - el t retardo)

[val lugar] =min(abs(tc+tRet-t0)); %Busco en ret+t1
y_t=w(lugar);
t=t0(lugar)-tRet; %t1

[val lugar] =min(abs(2*tc+tRet-t0));
y_t2=w(lugar);
t2=t0(lugar)-tRet;

[val lugar] =min(abs(3*tc+tRet-t0));
y_t3=w(lugar);
t3=t0(lugar)-tRet;

% K=y(inf)/U   => en este caso y(inf)= 198.2 y U=Vin=12 
k       =     198.2/12;     %ganancia de estado estacionario
 
 %METODO DE CHEN PARA W[rad/seg]
k1      =       (1/StepAmplitude)*y_t/k-1;
k2      =       (1/StepAmplitude)*y_t2/k-1;
k3      =       (1/StepAmplitude)*y_t3/k-1;
b       =       4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1   =       (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2   =       (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta    =       (2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1      =       (-t/log(alfa1));
T2      =       (-t/log(alfa2));
T3      =       (b*(T1-T2))+T1;
T1=real(T1);T2=real(T2); T3=real(T3); %importa sólo la parte real

%saco la funcion de transf
sys_w=tf(k,conv([T1 1],[T2 1]))
% como tiene un retardo el num debería multiplicarse por e^-Ls

[y1,t1,ent]=lsim(sys_w, Vi, t0, [0,0]);
%figure
%plot(t0,Vi);
%hold on
%plot(t1,y1,'k')
%plot([t+tRet t2+tRet t3+tRet],[y_t y_t2 y_t3],'o'); legend('DatosExcel','Modelo','Mediciones'); 


%METODO DE CHEN EN TL
tRet_tl=0.18500;
tc_tl= 0.18505-tRet_tl;

[val lugar] =min(abs(tc_tl+tRet_tl-t0)); %Busco en ret+t1
y_t_tl=w(lugar);
t_tl=t0(lugar)-tRet_tl; %t1

[val lugar] =min(abs(2*tc_tl+tRet_tl-t0));
y_t2_tl=w(lugar);
t2_tl=t0(lugar)-tRet_tl;

[val lugar] =min(abs(3*tc_tl+tRet_tl-t0));
y_t3_tl=w(lugar);
t3_tl=t0(lugar)-tRet_tl;

w1=198.2;                                  %Velocidad con TL=0
w2=164;                                %Velocidad angular con TL
Tlmax=max(Tl)*.93;
k_tl=(w1-w2)/Tlmax;           %ganancia para Tl

%corrección de y para Tl
yid_1=(y_t_tl - w1);
yid_2=(y_t2_tl - w1);
yid_3=(y_t3_tl - w1);

k1_tl=(1/Tlmax)*yid_1/k_tl -1;
k2_tl=(1/Tlmax)*yid_2/k_tl -1;
k3_tl=(1/Tlmax)*yid_3/k_tl -1;

b_tl=4*k1_tl^3*k3_tl-3*k1_tl^2*k2_tl^2-4*k2_tl^3+k3_tl^2+6*k1_tl*k2_tl*k3_tl;

if b_tl > 0
    alfa1_tl = (k1_tl * k2_tl + k3_tl - sqrt(b_tl)) / (2 * (k1_tl^2 + k2_tl));
    alfa2_tl = (k1_tl * k2_tl + k3_tl + sqrt(b_tl)) / (2 * (k1_tl^2 + k2_tl));
else
    alfa1_tl = (k1_tl * k2_tl + k3_tl - sqrt(complex(b_tl))) / (2 * (k1_tl^2 + k2_tl));
    alfa2_tl = (k1_tl * k2_tl + k3_tl + sqrt(complex(b_tl))) / (2 * (k1_tl^2 + k2_tl));
end

beta_tl   =  (2*k1_tl^3+3*k1_tl*k2_tl+k3_tl-sqrt(b_tl))/(sqrt(b_tl));
T1_tl      =  (-t_tl/log(alfa1_tl));
T2_tl      =  (-t_tl/log(alfa2_tl));
T3_tl      =  (b_tl*(T1-T2))+T1;
T1_tl=real(T1);T2_tl=real(T2); T3_tl=real(T3);

%saco la funcion de transf de Tl
sys_T=tf(k_tl,conv([T1_tl 1],[T2_tl 1]))

%simualción
[y2,t2_tl,ent]=lsim(sys_T, Vi, t0, [0,0]);



%CORROBORADOR DE PUNTOS extraidos
figure;
plot(t0,w);
hold on
plot(t+tRet,y_t,'o');
plot(t2+tRet,abs(y_t2),'o');
plot(t3+tRet,y_t3,'o');
plot(t_tl+tRet_tl,y_t_tl,'o');
plot(t2_tl+tRet_tl,y_t2_tl,'o');
plot(t3_tl+tRet_tl,y_t3_tl,'o');
hold off;

%Resultados
 %figure;
% hold on;
% plot(t0, w,'b');
% plot(t1,y1, 'g');
% xlim([0,0.5]);
% ylim([0,300]);
% figure;
% plot(t2_tl,y2, 'r');

%plot([t + tRet, t2 + tRet, t3 + tRet], [y_t, y_t2, y_t3], 'o');
%plot([t_tl + tRet_tl, t2_tl + tRet_tl, t3_tl + tRet_tl], [y_t_tl, y_t2_tl, y_t3_tl], 'o');

%legend('Datos','Modelado','Puntos Va', 'Puntos Tl');


 


