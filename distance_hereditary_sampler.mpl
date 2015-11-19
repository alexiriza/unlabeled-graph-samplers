# Unbiased sampler for (unlabeled) distance-hereditary graphs
# Given by the function: gen_dh_cycle_pointed(z) (0 < z <= 0.137935)

with(combstruct):

# Generating functions (formal power series) - see Sections 5.4.1 & 5.5.1 of paper for explanation

DistanceHereditary:={SX=Prod(Union(Z,K,SC),Set(Union(Z,K,SX),card>=1)),SC=Set(Union(Z,K,SX),card>=2),K=Set(Union(Z,SC,SX),card>=2)}:

acc:=2000:
kpoly:=add(count([K,DistanceHereditary,unlabeled],size=i)*z^i,i=0..acc):
dkpoly:=add(i*count([K,DistanceHereditary,unlabeled],size=i)*z^(i-1),i=0..acc):
d2kpoly:=add(i*(i-1)*count([K,DistanceHereditary,unlabeled],size=i)*z^(i-2),i=0..acc):
sxpoly:=add(count([SX,DistanceHereditary,unlabeled],size=i)*z^i,i=0..acc):
dsxpoly:=add(i*count([SX,DistanceHereditary,unlabeled],size=i)*z^(i-1),i=0..acc):
xdsxpoly:=expand(z*dsxpoly):
d2sxpoly:=add(i*(i-1)*count([SX,DistanceHereditary,unlabeled],size=i)*z^(i-2),i=0..acc):
scpoly:=add(count([SC,DistanceHereditary,unlabeled],size=i)*z^i,i=0..acc):
dscpoly:=add(i*count([SC,DistanceHereditary,unlabeled],size=i)*z^(i-1),i=0..acc):
d2scpoly:=add(i*(i-1)*count([SC,DistanceHereditary,unlabeled],size=i)*z^(i-2),i=0..acc):
xddhpoly:=expand(z+2*z^2+z*(sxpoly+scpoly+kpoly)+z^2*subs(z=z^2,dsxpoly)+z^2*subs(z=z^2,dscpoly)+z*(1+dkpoly+dscpoly+dsxpoly)-z*(1+dscpoly+dsxpoly)*(1+z+kpoly+scpoly+sxpoly)-z^2*(1+subs(z=z^2,dsxpoly)+subs(z=z^2,dscpoly))+(z+kpoly+scpoly)*(z*(1+dkpoly+dscpoly+dsxpoly)-z*(1+dkpoly+dsxpoly)*(1+z+kpoly+scpoly+sxpoly))):
xddhpoly:=add(coeff(xddhpoly,z,i)*z^i,i=0..acc): # Truncation
dhpoly:=add(coeff(xddhpoly, z, i)/i*z^i,i=1..acc):


# Generating functions (oracles)

sx:=proc(x)
  eval(sxpoly,z=x);
end proc:

dsx:=proc(x)
  eval(dsxpoly,z=x);
end proc:

xdsx:=proc(x)
  x*dsx(x);
end proc:

d2sx:=proc(x)
  eval(d2sxpoly,z=x);
end proc:

sc:=proc(x)
  eval(scpoly,z=x);
end proc:

dsc:=proc(x)
  eval(dscpoly,z=x);
end proc:

d2sc:=proc(x)
  eval(d2scpoly,z=x);
end proc:

k:=proc(x)
  eval(kpoly,z=x);
end proc:

dk:=proc(x)
  eval(dkpoly,z=x);
end proc:

d2k:=proc(x)
  eval(d2kpoly,z=x);
end proc:

xksx:=proc(x)
  x+k(x)+sx(x);
end proc:

xscsx:=proc(x)
  x+sc(x)+sx(x);
end proc:

dxksx:=proc(x)
  1+dk(x)+dsx(x);
end proc:

dxscsx:=proc(x)
  1+dsc(x)+dsx(x);
end proc:

xddh:=proc(x)
  eval(xddhpoly,z=x);
end proc:

dh:=proc(x)
  eval(dhpoly,z=x);
end proc:


# Samplers

read "proba_laws.mpl": # By Carine Pivoteau: contains implementations of
                       # standard probability laws (bernoulli, geometric,
                       # poisson, and logarithmic)

read "sampler_tools.mpl": # Contains Boltzmann samplers for sets and
                          # cycle-pointed sets, as well as functions to estimate
                          # radius of convergence of power series and to count
                          # the number of nodes of a certain type in a split
                          # tree string

gen_xksc:=proc(x)
  local n;
  n:=bern_option(x, k(x), sc(x));
  return piecewise(n=0, Z, n=1, gen_k(x), n=2, gen_sc(x));
end proc:

gen_xksx:=proc(x)
  local n;
  n:=bern_option(x, k(x), sx(x));
  return piecewise(n=0, Z, n=1, gen_k(x), n=2, gen_sx(x));
end proc:

gen_xscsx:=proc(x)
  local n;
  n:=bern_option(x, sc(x), sx(x));
  return piecewise(n=0, Z, n=1, gen_sc(x), n=2, gen_sx(x));
end proc:

gen_sx:=proc(x)
  local gamma,n;
  gamma:=gen_xksc(x);
  return SX(gamma, gen_set_of(x, xksx, gen_xksx, 1));
end proc:

gen_sc:=proc(x)
  return SC(gen_set_of(x, xksx, gen_xksx, 2));
end proc:

gen_k:=proc(x)
  return K(gen_set_of(x, xscsx, gen_xscsx, 2));
end proc:

gen_xksc_cycle_pointed:=proc(x)
  local n;
  n:=bern_option(x, x*dk(x), x*dsc(x));
  return piecewise(n=0, Z, n=1, gen_k_cycle_pointed(x), n=2, gen_sc_cycle_pointed(x));
end proc:

gen_xksx_cycle_pointed:=proc(x)
  local n;
  n:=bern_option(x, x*dk(x), x*dsx(x));
  return piecewise(n=0, Z, n=1, gen_k_cycle_pointed(x), n=2, gen_sx_cycle_pointed(x));
end proc:

gen_xscsx_cycle_pointed:=proc(x)
  local n;
  n:=bern_option(x, x*dsc(x), x*dsx(x));
  return piecewise(n=0, Z, n=1, gen_sc_cycle_pointed(x), n=2, gen_sx_cycle_pointed(x));
end proc:

gen_sx_cycle_pointed:=proc(x)
  local gamma,S,j;
  if bern(((x+x*dk(x)+x*dsc(x))*(evalf(exp(Sum((x^j + k(x^j) + sx(x^j))/j,j=1..infinity))) - 1))/(x*dsx(x)))=1 then
    gamma:=gen_xksc_cycle_pointed(x);
    S:=gen_set_of(x, xksx, gen_xksx, 1);
    return SX(gamma, S);
  else
    gamma:=gen_xksc(x);
    S:=gen_cycle_pointed_set_of(x, xksx, dxksx, gen_xksx, gen_xksx_cycle_pointed, 1, false);
    return SX(gamma, S);
  end if;
end proc:

gen_sc_cycle_pointed:=proc(x)
  return SC(gen_cycle_pointed_set_of(x, xksx, dxksx, gen_xksx, gen_xksx_cycle_pointed, 2, false));
end proc:

gen_k_cycle_pointed:=proc(x)
  return K(gen_cycle_pointed_set_of(x, xscsx, dxscsx, gen_xscsx, gen_xscsx_cycle_pointed, 2, false));
end proc:

gen_kscsx:=proc(x)
  local n;
  n:=bern_option(k(x), sc(x), sx(x));
  return piecewise(n=0, gen_k(x), n=1, gen_sc(x), n=2, gen_sx(x));
end proc:

gen_dh_cycle_pointed:=proc(x) # Boltzmann sampler for the cycle-pointed class of split trees of distance-hereditary graphs (hence an unbiased sampler for split trees of distance-hereditary graphs)
  local tau;
  if bern((x*(sx(x)+sc(x)+k(x)))/(xddh(x)))=1 then
    return Z(gen_kscsx(x))
  elif bern((x^2*dsx(x^2))/(xddh(x) - x*(sx(x)+sc(x)+k(x))))=1 then
    tau:=gen_sx_cycle_pointed(x^2);
    return e(tau, tau);
  elif bern((x^2*dsc(x^2))/(xddh(x) - x*(sx(x)+sc(x)+k(x)) - x^2*dsx(x^2)))=1 then
    tau:=gen_sc_cycle_pointed(x^2);
    return e(tau, tau);
  elif bern((x*(1+dk(x)+dsc(x)+dsx(x)) - x*(1+dsc(x)+dsx(x))*(1+x+k(x)+sc(x)+sx(x)) - x^2*(1+dsx(x^2)+dsc(x^2)))/(xddh(x) - x*(sx(x)+sc(x)+k(x)) - x^2*dsx(x^2) - x^2*dsc(x^2)))=1 then
    return KR(gen_cycle_pointed_set_of(x, xscsx, dxscsx, gen_xscsx, gen_xscsx_cycle_pointed, 3, true));
  else
    return SR(gen_xksc(x), gen_cycle_pointed_set_of(x, xksx, dxksx, gen_xksx, gen_xksx_cycle_pointed, 2, true));
  end if;
end proc:
