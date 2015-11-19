# Unbiased sampler for (unlabeled) three-leaf power graphs
# Given by the function: gen_3lp_cycle_pointed(z) (0 < z <= 0.259845)

with(combstruct):

# Generating functions (formal power series) - see Sections 5.4.2 & 5.5.2 of paper for explanation

ThreeLeafPower:={SC=Set(Union(A,SX),card>=2), SX=Prod(A, Set(Union(A, SX),card>=1)), A=Set(Z, card>=1)}:

acc:=2000:
apoly:=add(count([A,ThreeLeafPower,unlabeled],size=i)*z^i,i=0..acc):
dapoly:=add(i*count([A,ThreeLeafPower,unlabeled],size=i)*z^(i-1),i=0..acc):
xdapoly:=expand(z*dapoly):
d2apoly:=add(i*(i-1)*count([A,ThreeLeafPower,unlabeled],size=i)*z^(i-2),i=0..acc):
sxpoly:=add(count([SX,ThreeLeafPower,unlabeled],size=i)*z^i,i=0..acc):
dsxpoly:=add(i*count([SX,ThreeLeafPower,unlabeled],size=i)*z^(i-1),i=0..acc):
xdsxpoly:=expand(z*dsxpoly):
d2sxpoly:=add(i*(i-1)*count([SX,ThreeLeafPower,unlabeled],size=i)*z^(i-2),i=0..acc):
scpoly:=add(count([SC,ThreeLeafPower,unlabeled],size=i)*z^i,i=0..acc):
dscpoly:=add(i*count([SC,ThreeLeafPower,unlabeled],size=i)*z^(i-1),i=0..acc):
xdscpoly:=expand(z*dscpoly):
d2scpoly:=add(i*(i-1)*count([SC,ThreeLeafPower,unlabeled],size=i)*z^(i-2),i=0..acc):
xd3lppoly:=expand(z+2*z^2+z*(sxpoly+scpoly+apoly*(sxpoly+scpoly)+apoly-z)+z^2*subs(z=z^2,dsxpoly)+apoly*(z*(dapoly+dscpoly+dsxpoly)-z*(dapoly+dsxpoly)*(1+apoly+scpoly+sxpoly))+(z*apoly)*(apoly+1)-z^2+(scpoly+sxpoly)*(z*apoly)*(apoly+1)):
xd3lppoly:=add(coeff(xd3lppoly, z, i)*z^i, i=0..acc):
threelppoly:=add(coeff(xd3lppoly, z, i)/i*z^i, i=1..acc):


# Generating functions (oracles)

a:=proc(x)
  eval(apoly,z=x);
end proc:

da:=proc(x)
  eval(dapoly,z=x);
end proc:

xda:=proc(x)
  eval(xdapoly,z=x);
end proc:

d2a:=proc(x)
  eval(d2apoly,z=x);
end proc:

sx:=proc(x)
  eval(sxpoly,z=x);
end proc:

dsx:=proc(x)
  eval(dsxpoly,z=x);
end proc:

xdsx:=proc(x)
  eval(xdsxpoly, z=x);
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

xdsc:=proc(x)
  eval(xdscpoly,z=x);
end proc:

d2sc:=proc(x)
  eval(d2scpoly,z=x);
end proc:

asx:=proc(x)
  a(x)+sx(x);
end proc:

dasx:=proc(x)
  da(x)+dsx(x);
end proc:

id:=proc(x)
  x;
end proc:

did:=proc(x)
  1;
end proc:

xd3lp:=proc(x)
  eval(xd3lppoly,z=x);
end proc:

threelp:=proc(x)
  eval(threelppoly,z=x);
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

gen_x:=proc(x)
  return Z;
end proc:

gen_asx:=proc(x)
  if bern(a(x)/(a(x)+sx(x)))=1 then
    return gen_a(x);
  else
    return gen_sx(x);
  end if;
end proc:

gen_sc:=proc(x)
  return SC(gen_set_of(x, asx, gen_asx, 2));
end proc:

gen_sx:=proc(x)
  return SX(gen_a(x), gen_set_of(x, asx, gen_asx, 1));
end proc:

gen_a:=proc(x)
  local S;
  S:=gen_set_of(x, id, gen_x, 1);
  if S=Z then return Z
  else return K(S);
  end if;
end proc:

gen_x_cycle_pointed:=proc(x)
  return Z;
end proc:

gen_asx_cycle_pointed:=proc(x)
  if bern((x*da(x))/(x*da(x)+x*dsx(x)))=1 then
    return gen_a_cycle_pointed(x);
  else
    return gen_sx_cycle_pointed(x);
  end if;
end proc:

gen_sc_cycle_pointed:=proc(x)
  return SC(gen_cycle_pointed_set_of(x, asx, dasx, gen_asx, gen_asx_cycle_pointed, 2, false));
end proc:

gen_sx_cycle_pointed:=proc(x)
  local gamma, S, j;
  if bern((x*da(x)*(evalf(exp(Sum((a(x^j) + sx(x^j))/j,j=1..infinity))) - 1))/(x*dsx(x)))=1 then
    gamma:=gen_a_cycle_pointed(x);
    S:=gen_set_of(x, asx, gen_asx, 1);
    return SX(gamma, S);
  else
    gamma:=gen_a(x);
    S:=gen_cycle_pointed_set_of(x, asx, dasx, gen_asx, gen_asx_cycle_pointed, 1, false);
    return SX(gamma, S);
  end if;
end proc:

gen_a_cycle_pointed:=proc(x)
  local S;
  S:=gen_cycle_pointed_set_of(x, id, did, gen_x, gen_x_cycle_pointed, 1, false);
  if S=Z then return Z;
  else return K(S);
  end if;
end proc:

gen_sxsc:=proc(x)
  if bern(sx(x)/(sx(x)+sc(x)))=1 then
    return gen_sx(x);
  else
    return gen_sc(x);
  end if;
end proc:

gen_3lp_cycle_pointed:=proc(x) # Boltzmann sampler for the cycle-pointed class of split trees of three-leaf power graphs (hence an unbiased sampler for split trees of three-leaf power graphs)
  local A, S, part, build, tau;
  if bern((x*(sx(x)+sc(x)))/(xd3lp(x)))=1 then
    return Z(gen_sxsc(x));
  elif bern((x*a(x)*(sx(x)+sx(x)))/(xd3lp(x) - (x*(sx(x)+sc(x)))))=1 then
    A:=gen_a(x);
    S:=gen_sxsc(x);
    if A=Z then return Z(K(Z, S));
    else
      build:=NULL;
      for part in A do
        build:=build,part;
      od:
      return Z(K(build, S));
    end if;
  elif bern((x*(a(x) - x))/(xd3lp(x) - (x*(sx(x)+sc(x))) - (x*a(x)*(sx(x)+sx(x)))))=1 then
    return Z(K(gen_set_of(x, id, gen_x, 2)));
  elif bern((x^2*dsx(x^2))/(xd3lp(x) - (x*(sx(x)+sc(x))) - (x*a(x)*(sx(x)+sx(x))) - (x*(a(x) - x))))=1 then
    tau:=gen_sx_cycle_pointed(x^2);
    return e(tau, tau);
  elif bern((a(x)*(x*(da(x)+dsc(x)+dsx(x))-x*(da(x)+dsx(x))*(1+a(x)+sc(x)+sx(x))))/(xd3lp(x) - (x*(sx(x)+sc(x))) - (x*a(x)*(sx(x)+sx(x))) - (x*(a(x) - x)) - (x^2*dsx(x^2))))=1 then
    return SR(gen_a(x), gen_cycle_pointed_set_of(x, asx, dasx, gen_asx, gen_asx_cycle_pointed, 2, true));
  elif bern(((x*a(x))*(a(x)+1)-x^2)/(xd3lp(x) - (x*(sx(x)+sc(x))) - (x*a(x)*(sx(x)+sx(x))) - (x*(a(x) - x)) - (x^2*dsx(x^2)) - (a(x)*(x*(da(x)+dsc(x)+dsx(x))-x*(da(x)+dsx(x))*(1+a(x)+sc(x)+sx(x))))))=1 then
    return KR(gen_cycle_pointed_set_of(x, id, did, gen_x, gen_x_cycle_pointed, 3, true));
  else
    return KR(gen_sxsc(x), gen_cycle_pointed_set_of(x, id, did, gen_x, gen_x_cycle_pointed, 2, true));
  end if;
end proc:
