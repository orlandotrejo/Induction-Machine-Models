function px=maquina(t,x)

    global R L_1 G Ve VS Ler J k 
    ie=x(1); ir=x(2); wm=x(3); theta=x(4);
        
    ve = Ve*VS*[cos(t);cos(t-2*pi/3);cos(t-4*pi/3)];
    vr = 0+1i*0;
    
    pii = L_1*([ve;vr]-(R-1i*wm*G)*[ie;ir]); 
    pwm = (Ler*imag(ie*conj(ir))-k*wm^2)/J;
    px=[pii(1);pii(2);pwm;wm];
           
    
end
    
