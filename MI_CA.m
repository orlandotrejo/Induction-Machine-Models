%Modelo en coordenadas arbitrarias de la máquina de inducción

global Re Rr Le Lr VS Ler M Tr R L_1 Tm J

%Parametros de la MI en por unidad
Re=0.03; Rr=0.03; Le=0.13; Lr=0.13; Ler=0.02; M=3*Ler/2; Tr=Lr/Rr; Tm=0;
H=1; p=2; f=60; ws=2*pi*f/p; J=2*H*ws;

VS=sqrt(2/3)*[1 exp(1i*2*pi/3) exp(1i*4*pi/3)];

%Matrices para el modelo
R = [Re 0;0 Rr]; L=[Le M;M Lr]; L_1=inv(L);
x0 = [0,0,0,0,0]; t0=0; t=0:0.001*377:7*377;

% Llamada a la función de integración
[T,X]=ode45(@maquina_ca,t,x0);

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