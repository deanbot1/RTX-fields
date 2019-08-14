function f = et_code_rhs(t,y,pars)
    
    T     = y(1); 
    E     = y(2);
    D     = y(3);

    f     = zeros(3,1); % need to return a column vector
    g     = pars(1);
    r     = pars(2);
    kexp  = pars(3);
    NE    = pars(4);
    
    pars_f = pars(5:end);
   
    f(1)  = g*T - fcn_f(pars_f)*E^r*T;
    f(2)  = kexp*(NE-E) - r*fcn_f(pars_f)*E^r*T;
    f(3)  = fcn_f(pars_f)*E^r*T;

end

function f = f_molecular(pars_f)
    
    kdRCD20 = pars_f(1);
    kdRCD16 = pars_f(2);
    CD20    = pars_f(4);
    CD16    = pars_f(5);
    R       = pars_f(6);
    f = R.*(CD20/(kdRCD20+CD20) + CD16/(kdRCD16+CD16));

end

function g = fcn_f(pars_f)
    
    gamma = pars_f(3);
    g = gamma*f_molecular(pars_f);

end