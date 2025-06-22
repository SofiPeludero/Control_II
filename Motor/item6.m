clear all; clc;

TL_Max =0.0011;

X=-[0; 0;0];ii=0;
%t_etapa=t_D(10)-t_D(9); 
t_etapa =   1e-05;
titaRef=57.2958;
tF=.30; 

%Constantes del PID
%Kp=20;Ki=250;Kd=0.001;color_='r';
Kp=10;Ki=250;Kd=0.00001;color_='r';

Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
N=round(tF/t_etapa);
e=zeros(N,1);u=0;k=2;TL_ap=0;

for jj=1:N
    ii=ii+1;k=k+1;
    if jj*t_etapa>.07 %Primeros 5seg con Torque nulo
        TL_ap=TL_Max;
    end
    if jj*t_etapa<=.15 %Primeros 5seg con Torque nulo
        titaRef=pi/2;        
    end
    if jj*t_etapa>.15
        titaRef=-pi/2;
        TL_ap=TL_Max;
    end
%     X=modmotor_identificado(t_etapa, X, [u,TL_ap]); %Con ésto se ajusta
%     el PID, desde cálculo
    X=modmotor_inicial(t_etapa, X, [u,TL_ap]); %Con ésto se prueba (simil planta)
    e(k)=titaRef-X(3); %ERROR tita
    u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
    x1(ii)=X(1);% 
    x2(ii)=X(2);% 
    x3(ii)=X(3);%X=[ia;wr;titar];
    acc(ii)=u;
    TL_D1(ii)=TL_ap;
    titaRef_(ii)=titaRef;
end
th=0:t_etapa:tF;
t=th(1:ii);
figure(1);
subplot(2,2,1);hold on;
plot(t,titaRef_,'--' ,t,x3,color_,'LineWidth',1.3);title('Salida y, \theta_t');legend('Ref','\theta');legend('boxoff');grid on;
subplot(2,2,2);hold on;grid on;
plot(t,x1,color_,'LineWidth',1.3);title('Corriente i_t');
xlabel('Tiempo [Seg.]');
subplot(2,2,3);hold on;grid on;
plot(t,TL_D1,color_,'LineWidth',1.3);title('Torque T_L');
xlabel('Tiempo [Seg.]');
subplot(2,2,4);hold on;grid on;
plot(t,acc,color_,'LineWidth',1.3);title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');grid on;
