# Copyright © 2024, Inria (Fabrice Rouillier Fabrice.Rouillier@inria.fr)
# Permission is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the “Software”), to deal in 
# the Software without restriction, including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is furnished to do so, subject 
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind, express or implied, 
# including but not limited to the warranties of merchantability, fitness for a particular 
# purpose and noninfringement. In no event shall the authors or copyright holders X be 
# liable for any claim, damages or other liability, whether in an action of contract, 
# tort or otherwise, arising from, out of or in connection with the software or the 
# use or other dealings in the Software.
#
# Except as contained in this notice, the name of  Inria (Fabrice Rouillier)
# shall not be used in advertising or otherwise to promote the sale, use or 
# other dealings in this Software without prior written authorization from 
# the Inria (Fabrice Rouillier).

zds := module()

option package;
export
     rur,isolate:
local
     tri_sep_restrict_zp,biv_test_sep_gblex_zp,test_sep_gblex_zp,lin_extract_gblex_zp,
     param_shape_lemma_radicalize_zp,find_sep_elt,param_general_zp,pack_compute_seq,
     param_zerodim_internal_mm,param_zerodim,isolate_zerodim,
     gb_drl_zp,gb_lex_2_biv,gb_lex_zp,tri_split_zp,tri_lin_extract_zp,biv_lin_extract_gblex_zp:

# ******************************************************************************
# Groebner Engine: direct call to the limgb package 
# ******************************************************************************

gb_drl_zp:=proc(sys,vars,pr)
    return(libmgb:-gbasis(pr,0,vars,sys));
end:
gb_lex_zp:=proc(sys,vars,pr)
    local m,gro:
    gro:=libmgb:-gbasis(pr,0,vars,sys);
    return(libmgb:-fbasis(pr,0,vars,gro,vars));
end:
gb_lex_2_biv:=proc(gro,vars,v,pr)
    local m,vars2,p:
    vars2:=[v,vars[nops(vars)]]:
    return(libmgb:-fbasis(pr,nops(vars),vars,gro,vars2));
end:

# ******************************************************************************
# From a triangular set
# ******************************************************************************

# univ(x),tri(x,y) is a triangular set 
# lazy partial factorization of univ
# in a part above which y is separating
# and a part above which it is not

tri_sep_restrict_zp:=proc(univ,tri,pr)
    local vu,vb,pu,n,aa,a,i:
    vb:=op(indets(tri) minus indets(univ)):
    vu:=op(indets(univ)):
    pu:=univ:
    n:=degree(tri,vb):
    aa:=coeff(tri,vb,n-1):
    a:=-aa/n mod pr:
    for i from 2 to n do 
       aa:=-((n-i+1)/i)*aa*a mod pr;
       pu:=Gcd(aa-coeff(tri,vb,n-i),pu) mod pr:
    od:
    return(Quo(univ,pu,vu) mod pr,[pu,a]):
end:

# supposes that univ is squarefree
# univ(x),tri(x,y) is a triangular set 
# lazy partial factorization of univ
# in a part above which the leading coefficient
# of tri is invertible 
# and a part above which it is not

tri_split_zp:=proc(univ,tri,pr)
   local vu,vb,removed_fact,remaining_fact,reduced_tri,i,inv1,inv2:
   vu:=op(indets(univ)):
   vb:=op(indets(tri) minus indets(univ)):
   removed_fact:=Gcd(univ,lcoeff(tri,vb)) mod pr:
   remaining_fact:=Quo(univ,removed_fact,vu) mod pr:
   if (remaining_fact=1) then
      return(removed_fact,[]); 
   fi:
   reduced_tri:=vb^(degree(tri,vb))*(Gcdex(remaining_fact,lcoeff(tri,vb),vu,'inv1','inv2') mod pr):
   for i from 0 to degree(tri,vb)-1 do 
         reduced_tri:=(reduced_tri+(Rem(coeff(tri,vb,i)*inv2,remaining_fact,vu) mod pr)*vb^i) mod pr:
   od:
   return(removed_fact,[remaining_fact,reduced_tri]):
end:

#supposes that univ is squarefree
#univ(x),tri(x,y) is a triangular set 
#the function extracts a RUR from a triangular system
#supposing that univ separates the second coordinate

tri_lin_extract_zp:=proc(univ,tri,pr)
   local vu,vb,removed_fact,remaining_fact,reduced_tri,i,n:
   vu:=op(indets(univ)):
   vb:=op(indets(tri) minus indets(univ)):
   removed_fact:=Gcd(univ,lcoeff(tri,vb)) mod pr:
   remaining_fact:=Quo(univ,removed_fact,vu) mod pr:
   if (remaining_fact=1) then
      return(removed_fact,[]);
   fi:
   n:=degree(tri,vb):
   reduced_tri:=n*lcoeff(tri,vb)*vb+coeff(tri,vb,n-1);
   return(removed_fact,[remaining_fact,reduced_tri]):
end:

# ******************************************************************************
# From a bivariate gb lex
# ******************************************************************************

biv_test_sep_gblex_zp:=proc(gl,vars,pr)
  local pp,i,tr,uu,vu:
  vu:=op(indets(gl[1])):
  pp:=gl[1]:
  pp:=Quo(pp,(Gcd(pp,diff(pp,vu)) mod pr),vu) mod pr:
  for i from 2 to nops(gl) do
      pp,tr:=tri_split_zp(pp,gl[i],pr);
      if (nops(tr)>0) then 
         uu,tr:=tri_sep_restrict_zp(tr[1],tr[2],pr);
         if(uu<>1) then return(false): end:
      fi:
  od:
  return(true);
end:

biv_lin_extract_gblex_zp:=proc(gl,vars,pr)
  local pp,i,tr,uu,vu,prod,vb,bp,pp0,inv1,inv2:
  vu:=vars[2]:
  vb:=vars[1]:
  pp:=gl[1]:
  prod:=1:
  pp0:=Quo(pp,(Gcd(pp,diff(pp,vu)) mod pr),vu) mod pr:
  pp:=pp0:                        
  bp:=0:
  for i from 2 to nops(gl) do
      pp,tr:=tri_lin_extract_zp(pp,gl[i],pr);
      if (nops(tr)>0) then 
         bp:=bp+(expand(tr[2]*prod) mod pr) mod pr:
	 prod:=prod*tr[1]:
      fi:
  od:
  pp:=Gcdex(pp0,(Rem(coeff(bp,vb,1),pp0,vu) mod pr),vu,'inv1','inv2') mod pr:
  pp:=Rem(diff(pp0,vu)*Rem((Rem(coeff(bp,vb,0),pp0,vu) mod pr)*inv2,pp0,vu) mod pr,pp0,vu) mod pr:
  return([pp0,-pp mod pr]):
end:

# ******************************************************************************
# From a general gb lex
# ******************************************************************************

test_sep_gblex_zp:=proc(glex,vars,pr)
  local i, test_ok,vars2,gg3:
  i:=1:test_ok:=true:
  while ((i<nops(vars)) and (test_ok)) do
    vars2:=[vars[i],vars[nops(vars)]]:
    gg3:=gb_lex_2_biv(glex,vars,vars[i],pr):
    test_ok:=biv_test_sep_gblex_zp(gg3,vars2,pr):
    if (test_ok) then i:=i+1: fi:
  od:
  return(evalb(i=nops(vars))):
end:

lin_extract_gblex_zp:=proc(glex,vars,pr)
  local i, test_ok,vars2,gg3,tr,tr2:
  vars2:=[vars[1],vars[nops(vars)]]:
  gg3:=gb_lex_2_biv(glex,vars,vars[1],pr):
  tr:=biv_lin_extract_gblex_zp(gg3,vars2,pr);
  for i from 2 to (nops(vars)-1) do
    vars2:=[vars[i],vars[nops(vars)]]:
    gg3:=gb_lex_2_biv(glex,vars,vars[i],pr):
    tr2:=biv_lin_extract_gblex_zp(gg3,vars2,pr);
    tr:=[op(tr),tr2[2]]:
  od:
  return(tr):
end:

# ******************************************************************************
# From a system known to be shape lemma
# ******************************************************************************

param_shape_lemma_radicalize_zp:=proc(sys,vars,pr)
    local gg, gl, min_p,i:
    gl:=gb_lex_zp(sys,vars,pr):
    min_p:=Quo(gl[1],Gcd(gl[1],diff(gl[1],vars[nops(vars)])) mod pr,vars[nops(vars)]) mod pr:
    if (nops(gl)>nops(vars)) then error("Not Shape Lemma"): fi:
    return([min_p,seq(Rem((-coeff(gl[nops(vars)-i+1],vars[i],0))*diff(min_p,vars[nops(vars)]) mod pr,
            min_p,vars[nops(vars)]) mod pr,i=1..(nops(vars)-1))]):
end:

param_general_zp:=proc(sys,vars,pr)
    local gg, gl:
    gl:=gb_lex_zp(sys,vars,pr):
    return(lin_extract_gblex_zp(gl,vars,pr));
end:

pack_compute_seq:=proc(sys,vars,pri,my_function,pack)
   local pr,k,p_pr,l_pr,l_min_p,lift_p:
   pr:=pri:
   l_pr:=[]:
   p_pr:=1:
   l_min_p:=[]:
   for k from 1 to pack do
       pr:=prevprime(pr-1):
       l_min_p:=[op(l_min_p),my_function(sys,vars,pr)]:
       l_pr:=[op(l_pr),pr]:
       p_pr:=p_pr*pr:
   od:
   lift_p:=chrem(l_min_p,l_pr):
   return(pr,p_pr,lift_p):
end:

param_zerodim_internal_mm:=proc(sys,vars,param_engine_zp)
 local pr, prod_pr,i,reco_p,min_p,lift_p,lift_0,k,pack2,p_pr:
 pr:=prevprime(2^31):
 lift_0:=param_engine_zp(sys,vars,pr):
 prod_pr:=pr:
 reco_p:=FAIL:
 i:=1:
 while evalb(reco_p = FAIL) do
 #  pack2:=max(8,2^(ceil(log[2](i*10/100))));
   pack2:=max(2,ceil(i*10/100));
   print(i,pack2);
   pr,p_pr,lift_p:=pack_compute_seq(sys,vars,pr,param_engine_zp,pack2):
   lift_0:=chrem([lift_0,lift_p],[prod_pr,p_pr]):
   prod_pr:=prod_pr*p_pr:
   reco_p:=iratrecon(lift_0,prod_pr,scaled):
   i:=i+pack2:
 od:
 return(reco_p):
end:

# ******************************************************************************
# From a general system 
# ******************************************************************************

find_sep_elt:=proc(sys,vars)
  local pr,gdrl,glex,sys2,vars2,p,stop_comp,roll,i,is_cyclic,ii:
  roll:=rand(-10..10):
  pr:=prevprime(2^31);
  vars2:=vars:
  glex:=gb_lex_zp(sys,vars,pr):
  stop_comp:=test_sep_gblex_zp(glex,vars,pr):
  p:=vars[nops(vars)]:
  while (not(stop_comp)) do
    p:=_Z-(vars[1]+add(roll()*vars[i],i=2..nops(vars))):
    sys2:=[op(sys),p]:
    vars2:=[op(vars),_Z]:
    glex:=gb_lex_zp(sys2,vars2,pr):
    stop_comp:=test_sep_gblex_zp(glex,vars2,pr):
  end:
  if (nops(glex)=nops(vars)) then
     is_cyclic:=evalb(0=add(degree(glex[ii],vars2[nops(vars2)-ii+1])-1,ii=2..(nops(vars2))));
  else
     is_cyclic:=false:
  fi:
  return(p,is_cyclic);
end:

# returns ext,[coord(1),...,coord(n)] in the generl case
# and ext,[coord(1),...,coord(n-1)] whenever xn separates

param_zerodim:=proc(sys,vars)
     local p,is_cyclic:
     p,is_cyclic:=find_sep_elt(sys,vars):
     if (is_cyclic) then
        if (p=vars[nops(vars)]) then
	   print("Shape Lemma");
	   return(param_zerodim_internal_mm(sys,vars,param_shape_lemma_radicalize_zp));
	else
	   print("Cyclic quotient algebra");
	   return(param_zerodim_internal_mm([op(sys),p],[op(vars),_Z],param_shape_lemma_radicalize_zp)):
	fi:
     else
     	print("General case",vars,p);
	if (vars[nops(vars)]=p) then
	  return(param_zerodim_internal_mm(sys,vars,param_general_zp)):	  
	else
	  return(param_zerodim_internal_mm([op(sys),p],[op(vars),_Z],param_general_zp)):
	fi:
     fi:
end:

isolate_zerodim:=proc(sys,vars,prec)
     local rr,m,constr,i,de,iso:
     rr:=param_zerodim(sys,vars):
     m:=lcm(op(map(u->denom(u),rr))):
     constr:=numer([m*diff(rr[1],op(indets(rr[1]))),seq(m*rr[i],i=2..nops(vars))]):
     de:=expand(m*diff(rr[1],op(indets(rr[1])))):
     if (has(vars,op(indets(rr[1])))) then
     	return(fgbrs:-rs_isolate_rur(numer(rr[1]),de,[seq(expand(m*rr[i]),i=2..nops(rr)),expand(op(indets(rr[1]))*de)],op(indets(rr[1])),precision=prec));
     else
        return(fgbrs:-rs_isolate_rur(numer(rr[1]),de,[seq(expand(m*rr[i]),i=2..nops(rr))],op(indets(rr[1])),precision=prec));
     fi:
end:

# ******************************************************************************
# User Interface
# ******************************************************************************

rur:=proc(sys::depends(list(polynom(integer))),
          vars::list(name):=[op(indets(sys))])
     local rr,m,de,i:
     rr:=param_zerodim(sys,vars):
     m:=lcm(op(map(u->denom(u),rr))):
     de:=expand(m*diff(rr[1],op(indets(rr[1])))):
     if (has(vars,op(indets(rr[1])))) then 
        return(numer(rr[1])=0,[seq(vars[i]=expand(m*rr[i+1])/de,i=1..(nops(vars)-1)),vars[nops(vars)]=vars[nops(vars)]*de/de]):
     else
        return(numer(rr[1])=0,[seq(vars[i]=expand(m*rr[i+1])/de,i=1..(nops(vars)))]):        
     end:
end:

isolate:=proc(sys::depends(list(polynom(integer))),
              vars::list(name):=[op(indets(sys))],{precision::nonnegint:=Digits},$)
     local liso,prec,res,interv,i:
     liso:=isolate_zerodim(sys,vars,precision):
     res:=[]:
     for interv in liso do
         res:=[op(res),{seq(vars[i]=interv[i],i=1..nops(interv))}]:
     od:
     return(res):
end:

end module:

#savelib('zds',"/tmp");
