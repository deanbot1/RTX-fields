function f = et_code_rhs(t,y,pars)
    
    T     = y(1); 
    E     = y(2);
    D     = y(3);

    f     = zeros(3,1); % need to return a column vector
    g     = pars(1);
    r     = pars(2);
    kexp  = pars(3);
    NE    = pars(4);
    gamma = pars(5);
    f_rate= pars(6);
   
    f(1)  = g*T - fcn_f(f_rate,gamma)*E^r*T;
    f(2)  = kexp*(NE-E) - r*fcn_f(f_rate,gamma)*E^r*T;
    f(3)  = fcn_f(f_rate,gamma)*E^r*T;

end

function g = fcn_f(f_rate,gamma)
    
    g     = gamma*f_rate;

end