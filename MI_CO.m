% Modelo en campo orientado de la máquina de inducción

global Re Rr Le Lr VS Ler M Tr J

%Parametros de la MI en por unidad
Re=0.03; Rr=0.03; Le=0.13; Lr=0.13; Ler=0.02; M=3*Ler/2; Tr=Lr/Rr; p=2;
J=2*2*pi*60/p;

VS=sqrt(2/3)*[1 exp(1i*2*pi/3) exp(1i*4*pi/3)];

x0 = [0,0,0.05,0,0]; t0=0; t=0:0.001*377:7*377;

[T,X]=ode45(@maquina_co,t,x0);

Te = (M^2)*X(:,3).*X(:,2)/Lr;

% Transformacion a Coordenadas Primitivas

Ie_co = X(:,1)+1i*X(:,2);
Im_co = X(:,3).*exp(1i*X(:,4));

Iae = sqrt(2/3)*real(Ie_co.*exp(1i*X(:,4)));
Ibe = sqrt(2/3)*real(Ie_co.*exp(1i*X(:,4))*exp(1i*4*pi/3));
Ice = sqrt(2/3)*real(Ie_co.*exp(1i*X(:,4))*exp(1i*2*pi/3));

Ir_co = M*(X(:,3)-Ie_co)/Lr;
 
Iar = sqrt(2/3)*real(Ir_co.*exp(1i*X(:,4)));
Ibr = sqrt(2/3)*real(Ir_co.*exp(1i*X(:,4))*exp(1i*4*pi/3));
Icr = sqrt(2/3)*real(Ir_co.*exp(1i*X(:,4))*exp(1i*2*pi/3));