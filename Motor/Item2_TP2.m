clear all; clc; 

TL_Max =0.0011; %Valor obtenido en TP1

X=-[0; 0;0]; % Estado inicial del motor: X = [ia; wr; tita] (corriente, velocidad, posición)
ii=0;%t_etapa=t_D(10)-t_D(9); 
t_etapa =   1.0000e-05;
titaRef=pi/2;
tF=.30; 

%Constantes del PID
Kp=20;Ki=250;Kd=0.0001;color_='b';
%Kp=7;Ki=250;Kd=0.00001;color_='b';

Ts=t_etapa;
%discretizacion del PID
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

N=round(tF/t_etapa);
e=zeros(N,1);
u=0;
k=2;
TL_ap=0;

for jj=1:N
    ii=ii+1;
    k=k+1;
    if jj*t_etapa>.07 
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
    X=ModMotor2(t_etapa, X, [u,TL_ap]); %Con ésto se prueba (simil planta)
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
subplot(2,2,1); hold on;
plot(t,titaRef_,'--' ,t,x3,color_);title('Salida y, \theta_t');legend('Ref','\theta');legend('boxoff');grid on;
subplot(2,2,2);hold on;grid on;
plot(t,x1,color_);title('Corriente i_t');
xlabel('Tiempo [Seg.]');
subplot(2,2,3);hold on;grid on;
plot(t,TL_D1,color_);title('Torque T_L');
xlabel('Tiempo [Seg.]');
subplot(2,2,4);hold on;grid on;
plot(t,acc,color_);title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');grid on;

%Control en el Estacio de estados

% Según el tiempo de muestreo que es 1e-3, lo máximo que puede ir a la
% izquierda algún polo es 3/1e-3
disp('Valor a la izquierda máximo es:')
1/(3*t_etapa)
% 3.3333e+04
% # Constantes Identificadas
%Km = 1.0831e-02; Ki = 9.3023e-03; Ra = 18.3275; Laa = 4.0207e-04; Jm = 2.7398e-09; Bm = 0;
Laa=500e-6; Jm=2.5e-9; Ra=20; Bm=0;Ki=10e-3; Km=60.53e-3;%Motor inicial

%Estados: X=[ia;wr;titar]^T;

%Determinando las matrices A, B, B_T, y C
Mat_A=[-Ra/Laa -Km/Laa 0; Ki/Jm -Bm/Jm 0; 0 1 0 ];
Mat_B=[1/Laa;0;0]; Mat_B_T=[0;-1/Jm;0];
Mat_C=[0 0 1;0 1 0]; %Sale tita, y omega.

[n_va,d]=ss2tf(Mat_A,Mat_B,Mat_C(1,:),0); %Para verificar
[n_TL,d2]=ss2tf(Mat_A,Mat_B_T,Mat_C(1,:),0); %Para verificar
disp('Respecto de Va:')
Num=n_va/d(3)
Den=d/d(3)
disp('Respecto de TL:')
NumTL=n_TL/d2(3)
Den=d2/d2(3)

% En estados, con matrices ampliadas
% Cálculo del Controlador %
eig(Mat_A)

% Valor máximo a la izquierda es 3.3333e+04
Polos_deseados=[-20e3;-20e3; -1e3;-50e1]; %Cuatro polos
%Polos_deseados=[-2e4;-2e4; -1e4;-3e3]; %Cuatro polos

% Sistema ampliado para referencia no nula
Aa=[Mat_A,[0;0;0];-Mat_C(1,:),0]; %ojo, el valor deseado es tita
Ba=[Mat_B;0];   
Mat_M=[Ba Aa*Ba Aa^2*Ba Aa^3*Ba]; %Matriz Controlabilidad
disp('Debe ser 4:')
rank(Mat_M) %Debe ser 4

%  Cálculo del controlador por asignación de polos
auto_val=eig(Aa);
c_ai=poly(auto_val);
Mat_W=[c_ai(4)  c_ai(3) c_ai(2) 1 
       c_ai(3)  c_ai(2) 1 0 
       c_ai(2)  1     0 0;
       1  0     0 0];
Mat_T=Mat_M*Mat_W;
alfa_ia=poly(Polos_deseados); % K=place(Aa,Ba,Polos_deseados);
%Ka=fliplr(alfa_ia(2:length(alfa_ia))-c_ai(2:length(alfa_ia)))*inv(Mat_T);
Ka=lqr(Aa,Ba,diag([2e-2 1e-5 1e-5  1e5]),1e-5);

eig(Aa-Ba*Ka) %Deben ser los Polos_deseados

%Cálculo del Observador
%Sistema Dual
hh=[0.5;0.5];
Ad=Mat_A';
Bd=Mat_C'*hh;
Cd=Mat_B';
Polos_deseados_D=[-3.3e4;-3.2e4; -2e4]; %Tres polos
% [-3e4;-2e3; -1e3;-1e2];
%Polos_deseados_D = [-3e4; -2.5e4; -2e4];
% Sistema ampliado para referencia no nula
 
Mat_Md=[Bd Ad*Bd Ad^2*Bd  ]; %Matriz Controlabilidad
disp('Debe ser 3:')
rank(Mat_Md) %Debe ser 3, y mide sólo tita

%  Cálculo del controlador por asignación de polos
auto_val=eig(Ad);
c_ai=poly(auto_val);
Mat_Wd=[ c_ai(3)  c_ai(2) 1  
       c_ai(2)  1     0 ;
       1  0     0 ];
Mat_Td=Mat_Md*Mat_Wd;
alfa_ia=poly(Polos_deseados_D); % K=place(Aa,Ba,Polos_deseados);
Kd=fliplr(alfa_ia(2:length(alfa_ia))-c_ai(2:length(alfa_ia)))*inv(Mat_Td);
eig(Ad-Bd*Kd) %Deben ser los Polos_deseados
Ko=(hh*Kd)';
Ko=(lqr(Mat_A',Mat_C',diag([2e-2 1e0 1e3]),diag([1e-5, 1e-5])))';
eig(Mat_A-Ko*Mat_C)

% Simulación del comportamiento dinámico
X=[0;0;0]; psi=0; 
x_hat=[0;0;0];
N=round(tF/t_etapa);
e=zeros(N,1);
u=0;
k=2;
TL_ap=0;
ii=0; 
um=0.1;

X_sinZM = [0; 0; 0];  % Estado inicial para la simulación sin zona muerta
x1i = zeros(N,1); x2i = zeros(N,1); x3i = zeros(N,1); acci = zeros(N,1);

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
    Y=Mat_C*X;
%     X=modmotor_identificado(t_etapa, X, [u,TL_ap]); %Con ésto se ajusta
%     el PID, desde cálculo
    X=ModMotor(t_etapa, X, [u,TL_ap]); %Con ésto se prueba (simil planta)

    e(k)=titaRef-Y(1); %ERROR tita

    %     u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID       
    u=-Ka *[X;psi];
    %u=-Ka *[x_hat;psi]; %Con observación de estados
    ui=u;
    
    X_sinZM = ModMotor2(t_etapa, X_sinZM, [ui, TL_ap]);  % sin zona muerta
    
    if abs(ui)<um %um es zona muerta
        u=0;
    else
        u=ui-um*sign(u);
    end   
    
    x_hat_p=Mat_A*x_hat+Mat_B*u + Ko*(Y-Mat_C*x_hat);
    x_hat=x_hat+x_hat_p*t_etapa;
    psi=psi+e(k)*t_etapa;
    x1(ii)=X(1);% 
    x2(ii)=X(2);% 
    x3(ii)=X(3);%X=[ia;wr;titar];
    acc(ii)=u;%sin zona muerta
   
    x1i(ii) = X_sinZM(1); % corriente sin zona muerta
    x2i(ii) = X_sinZM(2); % velocidad sin zona muerta
    x3i(ii) = X_sinZM(3); % ángulo sin zona muerta
    acci(ii) = ui;        % entrada sin zona muerta
    %acci(ii)=ui; %con zona muerta
    
    TL_D1(ii)=TL_ap;
    titaRef_(ii)=titaRef;
end
% t=0:t_etapa:tF;
th=0:t_etapa:tF;
t=th(1:ii);
color_='r';

figure(2);
subplot(2,2,1);hold on;
plot(t,titaRef_,'--' ,t,x3,'m');title('Salida y, \theta_t');legend('Ref','\theta');legend('boxoff');grid on;
subplot(2,2,2);hold on;grid on;
plot(t,x1,'m');title('Corriente i_t');
xlabel('Tiempo [Seg.]');
subplot(2,2,3);hold on;grid on;
plot(t,TL_D1,'m');title('Torque T_L');
xlabel('Tiempo [Seg.]');
subplot(2,2,4);hold on;grid on;
plot(t,acc,'m');title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');grid on;

% Gráfica comparativa: con y sin zona muerta
th = 0:t_etapa:tF;
t = th(1:ii);

figure(3); clf;

subplot(2,2,1); hold on; grid on;
plot(t, titaRef_, '--k', 'DisplayName','Referencia');
plot(t, x3i, 'b', 'DisplayName','\theta_t sin zona muerta');
plot(t, x3, 'k', 'DisplayName','\theta_t con zona muerta');
title('Ángulo \theta_t');
legend('Location','best');

subplot(2,2,2); hold on; grid on;
plot(t, x1i, 'b', 'DisplayName','i_a sin zona muerta');
plot(t, x1, 'k', 'DisplayName','i_a con zona muerta');
title('Corriente i_t');
legend('Location','best');

subplot(2,2,3); hold on; grid on;
plot(t, x2i, 'b', 'DisplayName','\omega sin zona muerta');
plot(t, x2, 'k', 'DisplayName','\omega con zona muerta');
title('Velocidad angular \omega');
legend('Location','best');

subplot(2,2,4); hold on; grid on;
plot(t, acci, 'b', 'DisplayName','u(t) sin zona muerta');
plot(t, acc, 'k', 'DisplayName','u(t) con zona muerta');
% yline(um, '--r', 'Zona muerta');
% yline(-um, '--r');
title('Entrada de control u_t');
legend('Location','best');
xlabel('Tiempo [s]');
