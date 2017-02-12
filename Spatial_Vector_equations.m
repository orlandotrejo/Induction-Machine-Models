%Modelo dinamico de la maquina de induccion

global R L L_1 G Ve VS Ler J k

VS=sqrt(2/3)*[1 exp(1i*2*pi/3) exp(1i*4*pi/3)];

%Parametros de la MI
Re=0.03; Rr=0.03; Le=2.10; Lr=2.10; Ler=2.0; M=Ler; 
Ve=1;
%Parámetros de la carga
J=200; k=1.0

%Matrices para el modelo
R = [Re 0;0 Rr]; L=[Le M;M Lr]; G=[0 0;M Lr]; L_1=inv(L);

x0=[0,0,0,0];
[T,X]=ode45(@maquina,[0 1*377],x0);

Te=Ler*imag(X(:,1).*conj(X(:,2)));

figure(1)
plot(T,Te)
grid

figure(2)
plot(T,X(:,3))
grid

figure(3)
plot(real(X(:,1)),imag(X(:,1)))
axis equal
grid

iae=sqrt(2/3)*real(X(:,1))
ibe=sqrt(2/3)*real(X(:,1)*exp(1j*4*pi/3))

figure(4)
plot(T,iae,T,ibe,T,-iae-ibe)
grid

