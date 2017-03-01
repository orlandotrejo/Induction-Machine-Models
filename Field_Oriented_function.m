function px=Field_Oriented_Function(t,x)

global Re Rr Le Lr VS M Tr J

ide=x(1); iqe=x(2); im=x(3); delta=x(4); wm=x(5); 

if im < 0.05 
    im = 0.05;
end

ve = VS*[cos(t);cos(t-2*pi/3);cos(t-4*pi/3)];

T = [cos(delta), sin(delta); -sin(delta), cos(delta)];
               
vedq = T * [real(ve);imag(ve)];

vde = vedq(1);
vqe = vedq(2);

px(1) = ((Le - M*M/Lr)^(-1)) * (vde-(Re+Rr*M*M/(Lr*Lr))*ide) + ...
         wm*iqe + iqe*iqe/Tr/im + Rr*M*M*im/Lr/Lr;
    
px(2) = -wm*ide - ide*iqe/Tr/im - (Le - M*M/Lr)^(-1) * ...
        ((Re+Rr*M*M/Lr/Lr)*iqe - M*M/Lr*wm*im - vqe);
        
px(3) = (ide-im)/Tr;

px(4) = wm + iqe/(Tr*im);

px(5) = (M*M/Lr*im*iqe)/J;

px=px';

end
