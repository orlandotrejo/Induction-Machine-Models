function px=maquina_ca(t,x)

global Re Rr Le Lr VS Ler M Tr R L_1 Tm J

ie_delta=x(1); ir_delta=x(2); wr=x(3); theta=x(4); delta=x(5);

ve = VS*[cos(t);cos(t-2*pi/3);cos(t-4*pi/3)];
vr = 0;

% TRANSFORMACION A COORDENADAS DELTA

ve_delta = ve*exp(-1i*delta);
vr_delta = vr*exp(-1i*delta);

p_delta = wr; %Alineado con la fase a del rotor

% MODELO
G = [p_delta*Le, p_delta*M; (p_delta-wr)*M, (p_delta-wr)*Lr];

pi_delta = L_1 * ([ve_delta;vr_delta] - R*[ie_delta;ir_delta] - ...
           1i*G*[ie_delta;ir_delta]);

pie_delta = pi_delta(1,1);
pir_delta = pi_delta(2,1);

pwr=(M*imag(ie_delta*conj(ir_delta))-Tm)/J;
       
px = [pie_delta;pir_delta;pwr;wr;p_delta];

end
