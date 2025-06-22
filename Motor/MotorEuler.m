function [X]=modmotor_inicial(t_etapa, xant, accion)
h=5e-6;
%parametro motor item 5
Km = 1.0831e-02;
Ki = 9.3023e-03;
Ra = 18.3275;
Laa = 4.0207e-04;
J = 2.7398e-09;
Bm = 0;
%Laa=500e-6; J=2.5e-9; Ra=20; Bm=0;Ki=10e-3; Km=60.53e-3;
x=xant;
ia=xant(1);wr=xant(2);titar=xant(3);Va=accion(1);TL=accion(2);
for ii=1:t_etapa/h
    ia_p=-(Ra/Laa)*ia-(Km/Laa)*wr+(1/Laa)*Va;
    wr_p=(Ki/J)*ia-(Bm/J)*wr-(1/J)*TL;
    ia=ia+ia_p*h;
    wr=wr+wr_p*h;
    titar=titar+wr*h;
end
X=[ia;wr;titar];
