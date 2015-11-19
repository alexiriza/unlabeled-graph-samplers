#Basic Functions.
randomize();
random:=proc() (rand(1..10^20)()-1)/10.^20; end;

generalScheme:=proc(k_init, init, pNext, lambda)
local U, S, k, p;
    U:=random();
    p:=evalf(init);
    S:=p;
    k:=k_init;
    while U>S do
        p:=pNext(p,lambda,k);
        S:=S+p;
        k:=k+1
    end do;
    k
end;

#Bernoulli Law
bern:=proc(lambda)
local r;
    r:=random();
    if r < lambda
    then 1
    else 0
    end if
end;

#Geometric Law
nextGeom:=proc(p,lambda,k)
    lambda*p
end;

geom:=proc(lambda)
    generalScheme(0,(1-lambda),nextGeom,lambda)
end;

fast_geom:=proc(lambda)
     local U;
         U:=random();
         floor(log(U)/log(lambda));
     end;

#Poisson Law
nextPoiss:=proc(p,lambda,k)
    lambda*p/(k+1)
end;


poiss:=proc(lambda)
    generalScheme(0,evalf(exp(-lambda)),nextPoiss,lambda)
end;

non_zero_poiss:=proc(lambda)
    generalScheme(1,evalf(lambda/(exp(lambda)-1)),nextPoiss,lambda)
end;


#Logarithmic Law
nextLoga:=proc(p,lambda,k)
    lambda*p*(k)/(k+1)
end;

loga:=proc(lambda)
    generalScheme(1,-lambda/log(1-lambda),nextLoga,lambda)
end;
