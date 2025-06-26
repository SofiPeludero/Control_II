clear all; clc;
%
TL_Max =0.0011;
% break
X=-[0; 0;0];ii=0;
%t_etapa=t_D(10)-t_D(9); 
t_etapa =   1e-05;
titaRef=pi/2;
% tF=t_D(end); %Tiempo final del Excel.
tF=.30; %10 segundos

%Constantes del PID
Kp=10;Ki=500;Kd=0.00001;color_='b';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
N=round(tF/t_etapa);
e=zeros(N,1);u=0;k=2;TL_ap=0;
% break
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
%PID Motor con valores identificados en TP1
figure(1);
subplot(4,1,1);hold on;
plot(t,titaRef_,'--' ,t,x3,color_);title('Salida y, \theta_t');legend('Ref','\theta');legend('boxoff');grid on;
subplot(4,1,2);hold on;grid on;
plot(t,x1(1:30000),color_);title('Corriente i_t');
xlabel('Tiempo [Seg.]');
subplot(4,1,3);hold on;grid on;
plot(t,TL_D1,color_);title('Torque T_L');
xlabel('Tiempo [Seg.]');
subplot(4,1,4);hold on;grid on;
plot(t,acc(1:30000),color_);title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');grid on;

%%Control en el Estacio de estados

% Según el tiempo de muestreo que es 1e-3, lo máximo que puede ir a la
% izquierda algún polo es 3/1e-3
disp('Valor a la izquierda máximo es:')
1/(3*t_etapa)
% 3.3333e+04
% # Constantes Identificadas
Km = 1.0831e-02;
Ki = 9.3023e-03;
Ra = 18.3275;
Laa = 4.0207e-04;
J = 2.7398e-09;
Bm = 0;

%Estados: X=[ia;wr;titar]^T;
%Determinando las matrices A, B, B_T, y C
Mat_A=[-Ra/Laa -Km/Laa 0;
    Ki/J -Bm/J 0;
    0 1 0 ];
Mat_B=[1/Laa;0;0]; Mat_B_T=[0;-1/J;0];
Mat_C=[0 0 1]; %Sale tita.
[n_va,d]=ss2tf(Mat_A,Mat_B,Mat_C,0); %Para verificar
[n_TL,d2]=ss2tf(Mat_A,Mat_B_T,Mat_C,0); %Para verificar
disp('Respecto de Va:')
Num=n_va/d(3)
Den=d/d(3)
disp('Respecto de TL:')
NumTL=n_TL/d2(3)
Den=d2/d2(3)

%En estados, con matrices ampliadas
%%%%%%%% Cálculo del Controlador %%%%%%%%%%%%%%%%%%%
%  eig(Mat_A)=
% 1.0e+04 *
%    0.0000 + 0.0000i
%   -2.0000 + 0.9178i
%   -2.0000 - 0.9178i
%Valor máximo a la izquierda es 3.3333e+04
Polos_deseados=[-20e3;-20e3; -1e3;-50e1]; %Cuatro polos
 % Sistema ampliado para referencia no nula
Aa=[Mat_A,[0;0;0];-Mat_C,0];
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
Ka=fliplr(alfa_ia(2:length(alfa_ia))-c_ai(2:length(alfa_ia)))*inv(Mat_T);
% Ka=lqr(Aa,Ba,diag([1e0 1e0 1e0  1e1]),1e-1);
eig(Aa-Ba*Ka) %Deben ser los Polos_deseados
%%%%%%%%%Cálculo del Observador
%Sistema Dual
Ad=Mat_A';Bd=Mat_C';Cd=Mat_B';

Polos_deseados_D=[-3e4;-2.2e4; -2e4]; %Tres polos
% [-3e4;-2e3; -1e3;-1e2];r
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
Ko=Kd';
eig(Mat_A-Ko*Mat_C)

% Simulación del comportamiento dinámico
X=[0;0;0]; psi=0; 
x_hat=[0;0;0];
N=round(tF/t_etapa);
e=zeros(N,1);u=0;k=2;TL_ap=0;ii=0;k=2;
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
    X=modmotor_inicial(t_etapa, X, [u,TL_ap]); %Con ésto se prueba (simil planta)
    e(k)=titaRef-Y; %ERROR tita
%     u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID       
%     u=-Ka *[X;psi];
    u=-Ka *[x_hat;psi]; %Con observación de estados
    x_hat_p=Mat_A*x_hat+Mat_B*u + Ko*(Y-Mat_C*x_hat);
    x_hat=x_hat+x_hat_p*t_etapa;
    psi=psi+e(k)*t_etapa;
    x1(ii)=X(1);% 
    x2(ii)=X(2);% 
    x3(ii)=X(3);%X=[ia;wr;titar];
    acc(ii)=u;
    TL_D1(ii)=TL_ap;
    titaRef_(ii)=titaRef;
end

% t=0:t_etapa:tF;
th=0:t_etapa:tF;
t=th(1:ii);

%con Espcio de estados
figure(2);
subplot(4,1,1);hold on;
plot(t,titaRef_,'--' ,t,x3,color_);title('Salida y, \theta_t');legend('Ref','\theta');legend('boxoff');grid on;
subplot(4,1,2);hold on;grid on;
plot(t,x1,color_);title('Corriente i_t');
xlabel('Tiempo [Seg.]');
subplot(4,1,3);hold on;grid on;
plot(t,TL_D1,color_);title('Torque T_L');
xlabel('Tiempo [Seg.]');
subplot(4,1,4);hold on;grid on;
plot(t,acc,color_);title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');grid on;