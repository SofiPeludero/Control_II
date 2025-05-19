clear all; clc;

% Valores de los componentes
R = 47; L = 1e-6; C = 100e-9;

% Matrices del sistema
Mat_A = [-R/L -1/L; 1/C 0]; 
Mat_B = [1/L; 0];
Mat_C = [R 0];

% Cálculo de constantes de tiempo
Constantes_Tiempo = eig(Mat_A);

At = 2.137379556099644e-09;  % Paso de tiempo (10 veces menor que la más rápida)
Tf = 1.5e-3;                   % Tiempo final de manera que se pueda observar 

N = round(Tf/At);  % Número total de pasos de simulación
x = [0; 0];        % Condiciones iniciales

% Vector de tiempo corregido
t = linspace(0, Tf, N);
Vo = zeros(1, N);  % Salida inicializada

% Simulación con entrada alternante cada 1 ms
for k = 1:N
    % Alternar la entrada cada 1 ms
    u = 12 * (-1)^(floor(t(k) / 1e-3));
    
    % Aplicar la ecuación de estado
    xp = Mat_A * x + Mat_B * u;
    x = x + At * xp;
    
    % Guardar la salida
    Vo(k) = Mat_C * x;
end

% Graficar la respuesta del sistema
figure;
plot(t, Vo, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Voltaje');
title('Respuesta del sistema RLC con entrada alternante');
grid on;
