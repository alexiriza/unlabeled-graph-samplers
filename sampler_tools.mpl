# Radius of convergence

radius:=proc(poly) # Computes an estimate of the radius of convergence of the specified generating function power series, using the Domb-Sykes method (See p. 37 of paper)
  local X, Y, deg;
  deg:=degree(poly, z);
  X:=Vector([seq](1/n,n=(acc/2-1)..(acc-1)));
  Y:=Vector([seq](evalf(coeff(poly, z, n+1)/coeff(poly, z, n)),n=(acc/2-1)..(acc-1)));
  return 1/Statistics[LinearFit]([t->1, t->t], X, Y)[1];
end proc:


# Counting parameters (cliques, stars, leaves)

count_param:=proc(tree, param) # Returns the number of occurrences of the specified param in the specified tree
  local param_string, s, subtree, header;
  if tree=Z then if param=Z then return 1; else return 0; end if; end if;
  s:=0;
  for subtree in tree do
    s:=s+count_param(subtree, param);
  od;
  param_string:=convert(param, string);
  header:=StringTools[Split](convert(tree, string), "(")[1];
  if header=param_string then return s+1; else return s; end if;
end proc:

num_cliques:=proc(tree) # Returns the number of clique nodes in the specified tree
  return count_param(tree, K) + count_param(tree, KR);
end proc:

num_stars:=proc(tree) # Returns the number of star nodes in the specified tree
  return count_param(tree, SX) + count_param(tree, SC) + count_param(tree, SR);
end proc:

num_Z:=proc(tree) # Returns the number of leaves in the specified tree (i.e. the size of the tree)
  return count_param(tree, Z);
end proc:


# Boltzmann samplers for sets and cycle-pointed sets

read "proba_laws.mpl": # By Carine Pivoteau: contains implementations of
                       # standard probability laws (bernoulli, geometric,
                       # poisson, and logarithmic)

max_index:=proc(A, x) # MAX_INDEX generator, from 'Boltzmann Samplers, Polya Theory, and Cycle Pointing' p. 32
  local U,v,k,j;
  U := random();
  k := 0;
  v := 1/evalf(exp(Sum(A(x^j)/j,j=1..infinity)));
  while v < U do
    k := k+1;
    v := v * (evalf(exp(A(x^k)/k)));
  od;
  return k;
end proc:

max_index_geq_1:=proc(A, x) # MAX_INDEX generator, conditioned on the output being >= 1
  local U,v,c,k,j;
  U:=random();
  k:=1;
  v:=evalf(exp(A(x)));
  c:=1/evalf(exp(Sum(A(x^j)/j,j=1..infinity))-1);
  while (v-1)*c < U do
    k:=k+1;
    v:=v*evalf(exp(A(x^k)/k));
  od;
  return k;
end proc:

gen_set_of:=proc(x, inner_oracle, inner_gen, length_lower_bound) # Boltzmann sampler for an unlabeled set, with option to specify lower bound on the size of the set
  local S,i,j,J,k,kJ,t,copy;
  S:=NULL;
  if length_lower_bound <= 0 then
    J:=max_index(inner_oracle, x);
    if J <> 0 then
      kJ:=non_zero_poiss(inner_oracle(x^J)/J);
    end if;
  elif length_lower_bound = 1 then
    J:=max_index_geq_1(inner_oracle, x);
    kJ:=non_zero_poiss(inner_oracle(x^J)/J);
  elif length_lower_bound = 2 then
    do
      J:=max_index_geq_1(inner_oracle,x);
      kJ:=non_zero_poiss(inner_oracle(x^J)/J);
      if J>1 or kJ>1 then
        break;
      end if;
    od;
  else
    error "not implemented"
  end if;
  for j from 1 to J do
    if j < J then
       k:=poiss(inner_oracle(x^j)/j);
    else
      k:=kJ;
    end if;
    for i from 1 to k do
      t:=inner_gen(x^j);
      for copy from 1 to j do
        S:=S,t;
      od;
    od;
  od;
  return S
end proc:

bern_option:=proc(a, b, c)
  if bern(a/(a+b+c))=1 then 0;
  elif bern(b/(b+c))=1 then 1;
  else 2;
  end if;
end proc:

root_cycle_size:=proc(x, d_inner_oracle) # ROOT_CYCLE_SIZE generator, from 'Boltzmann Samplers, Polya Theory, and Cycle Pointing' p. 36
  local U,c,k,v,j;
  U:=random();
  k:=0;
  v:=0;
  c:=evalf(Sum(x^j*d_inner_oracle(x^j),j=1..infinity));
  while v < c*U do
    k:=k+1;
    v:=v+x^k*d_inner_oracle(x^k);
  end do;
  return k;
end proc:

root_cycle_size_geq_2:=proc(x, d_inner_oracle) # ROOT_CYCLE_SIZE generator, conditioned on the output being >= 2
  local U,c,k,v,j;
  U:=random();
  k:=1;
  v:=0;
  c:=evalf(Sum(x^j*d_inner_oracle(x^j),j=2..infinity));
  while v < c*U do
    k:=k+1;
    v:=v+x^k*d_inner_oracle(x^k);
  end do;
  return k;
end proc:

gen_cycle_pointed_set_of:=proc(x, inner_oracle, d_inner_oracle, inner_gen, inner_gen_cycle_pointed, length_lower_bound, symmetric) # Boltzmann sampler for an unlabeled cycle-pointed set, with options to specify lower bound on the size of the set and to restrict to a symmetric cycle-pointed set (i.e. one with marked cycle of length >= 2)
  local K,gamma,S,copy,j;
  S:=NULL;
  if symmetric then
    if length_lower_bound <= 2 then
      K:=root_cycle_size_geq_2(x, d_inner_oracle);
      S:=gen_set_of(x, inner_oracle, inner_gen, 0);
    elif length_lower_bound = 3 then
      do
        K:=root_cycle_size_geq_2(x, d_inner_oracle);
        S:=gen_set_of(x, inner_oracle, inner_gen, 0);
        if K>2 or S<>NULL then
          break;
        end if;
      od;
    else
      error "not implemented"
    end if;
  else
    if length_lower_bound <= 1 then
      K:=root_cycle_size(x, d_inner_oracle);
      S:=gen_set_of(x, inner_oracle, inner_gen, 0);
    elif length_lower_bound = 2 then
      do
        K:=root_cycle_size(x, d_inner_oracle);
        S:=gen_set_of(x, inner_oracle, inner_gen, 0);
        if K>1 or S<>NULL then
          break;
        end if;
      od;
    else
      error "not implemented"
    end if;
  end if;
  gamma:=inner_gen_cycle_pointed(x^K);
  for copy from 1 to K do
    S:=S,gamma;
  od;
  return S
end proc:
