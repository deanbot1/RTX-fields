function f = reactions(A0,R0,C0,kon20,kon16,koff20,koff16,tf,nr_time)

    y0 = [A0,R0,C0,0,0,0];
    
    tspan = linspace(0,tf,nr_time);
    [T,Y] = ode45(@(t,y) model1(t,y,kon20,kon16,koff20,koff16),tspan,y0);
    
    f = Y(end,6);

end

function dydt = model1(t,y,kon20,kon16,koff20,koff16)

A=y(1);
B=y(2);
C=y(3);
AB=y(4);
BC=y(5);
ABC=y(6);

dydt(1) = -kon20*B.*A+koff20*AB-kon20*A.*BC+koff20*ABC;
dydt(2) = koff20*AB+koff16*BC-kon20*A.*B-kon16*B.*C;
dydt(3) = -kon16*B.*C+koff16*B.*C-kon16*C.*AB + koff16*ABC;
dydt(4) = kon20*B.*A-kon16*C.*AB-koff20*AB+koff16*ABC;
dydt(5) = kon16*B.*C - kon20*A.*BC-koff16*BC+koff20*ABC;
dydt(6) = kon16*C.*AB + kon20*A.*BC-koff16*ABC-koff20*ABC;
dydt    = [dydt(1);dydt(2);dydt(3);dydt(4);dydt(5);dydt(6)];

end