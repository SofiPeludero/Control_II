clear all; clc; 

[z1]=xlsread('Curvas_Medidas_RLC_2024');
t0=z1(:,1);
I=z1(:,2);
Vc=z1(:,3);
Vin=z1(:,4);

subplot(3,1,1);hold on;
plot(t0,I,'r');title('Corriente Inductor I_L');grid;
subplot(3,1,2);hold on;
plot(t0,Vc,'r');title('Tension Capacitor V_C');grid;
plot(t0,Vin,'b')
xlabel('Tiempo [Seg.]');
subplot(3,1,3);hold on;
plot(t0,Vin,'r');title('Tension Entrada V_i');grid;