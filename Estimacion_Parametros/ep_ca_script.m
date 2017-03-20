%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA PRINCIPAL ESTIMACION DE PARAMETROS DE LA MI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Re Rr Le Lr VS Ler M Tr R L_1 J kl

%Parametros de la MI en por unidad
Re=0.03; Rr=0.03; Le=0.13; Lr=0.13; Ler=0.02; M=3*Ler/2; Tr=Lr/Rr; H=1;
kl=1; p=2; f=60; ws=2*pi*f/p; J=2*H*ws;

VS=sqrt(2/3)*[1 exp(1i*2*pi/3) exp(1i*4*pi/3)];

%Matrices para el modelo
R = [Re 0;0 Rr]; L=[Le M;M Lr]; L_1=inv(L);
x0 = [0,0,0,0,0]; t0=0; t=0:0.001*377:20*377;

% Llamada a la función de integración
[T,X]=ode45(@maquina_ep_ca,t,x0);

ie_delta = X(:,1);
ir_delta = X(:,2);
wr = X(:,3);
theta = X(:,4);
delta = X(:,5);

%Cálculo de Par
Te=M*imag(ie_delta.*conj(ir_delta));

% Paso a coordenadas primitivas
iae=sqrt(2/3)*real(ie_delta.*exp(1i*delta));
ibe=sqrt(2/3)*real(ie_delta.*exp(1i*delta)*exp(1i*4*pi/3));
ice=sqrt(2/3)*real(ie_delta.*exp(1i*delta)*exp(1i*2*pi/3));
iar=sqrt(2/3)*real(ir_delta.*exp(1i*(delta-theta)));
ibr=sqrt(2/3)*real(ir_delta.*exp(1i*(delta-theta)).*exp(1i*4*pi/3));
icr=sqrt(2/3)*real(ir_delta.*exp(1i*(delta-theta)).*exp(1i*2*pi/3));

% Variables de estado
Ve=exp(1i*T)/sqrt(2/3);
Ie=sqrt(2/3)*(iae+ibe*exp(1i*2*pi/3)+ice*exp(1i*4*pi/3));
Ir=sqrt(2/3)*(iar+ibr*exp(1i*2*pi/3)+icr*exp(1i*4*pi/3));

% Derivadas de las variables de estado
pVe = 1i*Ve;
pIe = 1i*Ie;
p2Ie = -Ie;

% pVe = diff(Ve)/(0.001*377);
% pIe = diff(Ie)/(0.001*377);
% p2Ie = diff(diff(Ie))/(0.001*377);

% Tres medidas para distintas cargas (puntos)
punto_op = [ 0.2 0.7 0.9];

% Inicialización de Variables Intermedias para la Regresión
wi=zeros(length(punto_op),3);
hi=zeros(length(punto_op),1);
NUM = zeros(3,1);
DEN = zeros(3);

% Calculo de coeficientes
for temp=1:length(punto_op)

    m=floor(punto_op(temp)*length(Ve));
    
    % Cuando son 5 coeficientes
    %wi(temp,:)=[(p2Ie(m)-1i*wr(m)*pIe(m)), pIe(m), -1i*wr(m)*Ir(m), -Ve(m), Ie(m)];
    %hi(temp)=pVe(m)-1i*wr(m)*Ve(m);    
    
    % Cuando son 3 coeficientes
    wi(temp,:)=[(p2Ie(m)-1i*wr(m)*pIe(m)), pIe(m), (Re*Ie(m)-Ve(m))];
    hi(temp)= pVe(m)-1i*wr(m)*Ve(m)+1i*Re*wr(m)*Ie(m); 
    
    NUM = wi(temp,:)'*hi(temp) + NUM;
    DEN = wi(temp,:)'*wi(temp,:) + DEN;
    
end

K = DEN\NUM;

Rr_est = real((K(2)-Re));
Lr_est = real((K(2)-Re)/K(3));
Le_est = Lr_est;
Ler_est= real(sqrt((Le_est-K(1))*Lr_est));

var = [Rr_est,Lr_est,Ler_est];

disp(var)