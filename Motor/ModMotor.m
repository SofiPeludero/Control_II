function [X]=ModMotor(dt, xant, accion)
    Km = 1.0831e-02; 
    Ki = 9.3023e-03; 
    Ra = 18.3275;
    Laa = 4.0207e-04; 
    J = 2.7398e-09; 
    B = 0.0;
    Va = accion(1); 
    TL = accion(2);
    h = 1e-6; % Paso interno para mejor precisión

    ia = xant(1); wr = xant(2); theta = xant(3);
    for i = 1:round(dt/h)
        dia = (Va - Ra*ia - Ki*wr)/Laa;
        dwr = (Km*ia - TL - B*wr)/J;

        ia = ia + h*dia;
        wr = wr + h*dwr;
        theta = theta + h*wr;
    end
    X = [ia; wr; theta];
end
