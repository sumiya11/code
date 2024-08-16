# Copyright © 2024, Alexander Demin (asdemin_2@edu.hse.ru) and Fabrice Rouillier (Fabrice.Rouillier@inria.fr)
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
# Except as contained in this notice, the names of  Alexander Demin  and  Fabrice Rouillier
# shall not be used in advertising or otherwise to promote the sale, use or 
# other dealings in this Software without prior written authorization from 
# Alexander Demin  and  Fabrice Rouillier.

module RationalUnivariateRepresentation

import AbstractAlgebra,Groebner,Nemo,Primes,Random
import AbstractAlgebra: QQ, polynomial_ring
import Base.Threads: nthreads
import Distributed
import ThreadPools: @tspawnat

export zdim_parameterization, QQ, polynomial_ring, to_file

#####################################################
# Logging
#####################################################

const _verbose = Ref{Bool}(true)

rur_print(xs...) = rur_print(stdout, xs...)
function rur_print(io::IO, xs...)
    _verbose[] && print(io, xs...)
    Base.flush(io)
    nothing
end

#####################################################
# Add on to AbstractAlgebra
#####################################################

function convert_to_mpol_UInt32(sys_z::Vector{AbstractAlgebra.Generic.MPoly{BigInt}},pr::UInt32)
   sys_Int32=AbstractAlgebra.Generic.MPoly{AbstractAlgebra.GFElem{Int64}}[]
   for p in sys_z push!(sys_Int32,AbstractAlgebra.change_coefficient_ring(AbstractAlgebra.GF(Int64(pr)),p)) end
   return(sys_Int32)
end

function extract_raw_data(sys_z::Vector{AbstractAlgebra.Generic.MPoly{BigInt}})
    expvecs=Vector{Vector{PP}}()
    cfs_zz=Vector{Vector{BigInt}}()
    for i in 1:length(sys_z)
        push!(expvecs, map(PP,collect(AbstractAlgebra.exponent_vectors(sys_z[i]))))
        push!(cfs_zz, map(BigInt,collect(AbstractAlgebra.coefficients(sys_z[i]))))
    end
    return expvecs,cfs_zz
 end

function coeff_mod_p(x::AbstractAlgebra.GFElem{Int32},pr::UInt32)::UInt32
    return(UInt32(AbstractAlgebra.data(x)))
end

function coeff_mod_p(x::AbstractAlgebra.GFElem{Int64},pr::UInt32)::UInt32
    return(UInt32(AbstractAlgebra.data(x)))
end

function sys_mod_p(sys::Vector{AbstractAlgebra.Generic.MPoly{T}},pr::UInt32)::Vector{PolUInt32} where {T}
    res=PolUInt32[]
    for p in sys
        push!(res,PolUInt32(map(u->PP(u),collect(AbstractAlgebra.exponent_vectors(p))),map(u->coeff_mod_p(u,pr),collect(AbstractAlgebra.coefficients(p)))))
    end
    return(res)
end

#####################################################
# Add on to Groebner
#####################################################

function groebner_linform(sys_Int32,linform)
    gro=Groebner.groebner(sys_Int32,ordering=Groebner.DegRevLex(),threaded=:no);
    return gro
end

function groebner_learn_linform(sys_Int32,linform)
    graph,gro=Groebner.groebner_learn(sys_Int32,ordering=Groebner.DegRevLex(),threaded=:no);
    return graph,gro
end

function kbase_linform(gro,pr,linform)
    g=sys_mod_p(gro,pr);
    ltg=map(u->u.exp[1],g);
    R=AbstractAlgebra.parent(gro[1])
    q1=compute_quotient_basis(ltg);
    quo=sort(map(u->R([1],[Vector{Int64}(u.data)]),q1))
    return quo,g
end

#####################################################
# internal polynomials
#####################################################

struct PP
    data::Vector{Int32}
    function PP(v::Vector{T}) where {T}
        new(Vector{Int32}(v))
    end
end

function pp_isless_drl(ea::PP, eb::PP)
    if sum(ea.data) < sum(eb.data)
        return true
    elseif sum(ea.data) != sum(eb.data)
        return false
    end
    i = length(a)
    @inbounds while i > 1 && a[i] == b[i]
        i -= 1
    end
    @inbounds return a[i] > b[i]
end

struct PolUInt32
   exp::Vector{PP}
   co::Vector{UInt32}
end

function monomial(p::PolUInt32,i::Int32)::PP
    return(p.exp[i])
end

function coeff(p::PolUInt32,i::Int32)::UInt32
    return(p.co[i])
end

mutable struct StackVect
   pos::Int32
   mon::PP
   prev::Int32
   var::Int32
end

function coeff_mod_p(x::BigInt,pr::UInt32)
    if (x<0) return(pr-UInt32((-x)%pr))
    elseif return(UInt32(x%pr)) end
end

function total_degree(pp::PP)::Int64
    sum(pp.data)
end

# tests if pp1 is divisible by pp2
# then returns the greatest variable
# that appear in the divisor
@inline function divides(pp1::PP,pp2::PP)
    a, b = pp1.data, pp2.data
    var = 0
    @inbounds for j in 1:length(a)
        if a[j] > b[j]
            var = j
        end
        if a[j] < b[j]
            return false, var
        end
    end
    return true, var
end

function find_divisor(m::PP,lm::Vector{PP})
    pos=length(lm);
    @inbounds for j in pos:-1:1
       te,v=divides(m,lm[j])
       if te return(j) end
    end
    return(0)
end
    
function find_in_list(m::PP,lm::Vector{PP})
    pos=length(lm);
    @inbounds for j in pos:-1:1
       if m.data==lm[j].data return(j) end
    end
    return(0)
end

function compute_quotient_basis(ltg::Vector{PP})
    if length(ltg)==0 return([]) end
    nbv=length(ltg[1].data)
    quo=Vector{PP}()
    todo=Vector{PP}([PP([Int32(0) for i=1:nbv])])
    inspected=Vector{PP}([PP([Int32(0) for i=1:nbv])])
    pos=1
    while (length(todo)>0)
        m=popfirst!(todo)
        pos=find_divisor(m,ltg)
        if pos==0 
            push!(quo,m)
            @inbounds for i in 1:nbv
                pp=PP(m.data)
                pp.data[i]+=1
                if (find_in_list(pp,inspected)==0)
                   push!(todo,pp)
                   push!(inspected,pp)
                end           
            end
        end
    end
    rur_print("(",length(quo),",",length(inspected),")")
    return(quo)
end

#####################################################
# Multiplication Matrices
#####################################################

# test if a monomial is already in the border and, if not, try to find a divisor
# of total degree -1. the function supposes that it does not belong to the
# quotient variant where we look for the smallest divisor
function find_in_border(m::PP,
                        t_xw::Vector{StackVect})
   pos=length(t_xw)
   res_flag=0
   res_dd=0
   res_pos=0
   # On y va par degrÃ©s dÃ©croissants pour limiter
   # le nombre de tests
   tm=total_degree(m)
   while (pos>0)
      mm=t_xw[pos].mon
      tmm=total_degree(mm)
      if (tmm==tm)
         if (m==mm) return(1,m,pos) end
      elseif (tmm==tm-1)
            # test if mm is not in the quotient basis
            if (!((t_xw[pos].prev>0) && (t_xw[pos].var==0))) 
               flag,dd=divides(m,mm)
               if (flag)
                  res_flag=2 
                  res_dd=dd
                  res_pos=pos
               end
            end
      else
            # the total degree is to high to find a predecessor
            if (res_flag>0) 
               return(res_flag,res_dd,res_pos)
            else
               rur_print("\nError find in border : ",tm," ",tmm)
               return(0,0,0)
            end
      end
      pos=pos-1
   end
   if (res_flag>0)
      return(res_flag,res_dd,res_pos)
   else 
      rur_print("\n\nError in finding a predecessor")
      return(0,0,0)
   end
end

function mul_pp_by_var!(m::PP,v_index::Int64)
    mm=PP(m.data)
    mm.data[v_index]=mm.data[v_index]+1
    return(mm)
end

function prepare_table_mxi(ltg::Vector{PP},
                           kb::Vector{PP})
   nbv=length(ltg[1].data)
   tablex=[[Int32(0) for i=1:length(kb)] for j=1:nbv]
   general_stack=Vector{StackVect}()
   nb_stack=0    
   for j in eachindex(kb)
          m=PP(kb[j].data)
          for ii in 1:nbv
              i=nbv-ii+1   
	          nm=mul_pp_by_var!(m,i)
	          pos=findfirst(item -> item.data == nm.data, kb)
	          if isnothing(pos)
	  	         pos=findfirst(item -> item.data == nm.data, ltg)
     	         if isnothing(pos)
	               flag,v,prec=find_in_border(nm,general_stack)
	               if (flag==1) tablex[i][j]=Int32(prec)
		           elseif (flag==2)
		             nb_stack=nb_stack+1
	                 tablex[i][j]=Int32(nb_stack)
                         push!(general_stack,StackVect(nb_stack,PP(nm.data),Int32(prec),Int32(v)))
	             else rur_print("\n*** Error search table ***\n") end
	         else
		     # nm is a leading monomial of an element of the Gb
		     # we insert it with the flags prev=0 and var=pos in the gb
                     prev=findlast(item-> item.mon.data==nm.data,general_stack)
                     if isnothing(prev)
		                 nb_stack=nb_stack+1
		                 tablex[i][j]=Int32(nb_stack)
		                push!(general_stack,StackVect(nb_stack,PP(ltg[pos].data),Int32(0),Int32(pos)))
                     else tablex[i][j]=Int32(prev) end
		     end
	         else
		          # nm is a element of the quotient basis
		          # we insert it with the flags prev=pos and var=0
                  prev=findlast(item-> item.mon.data==nm.data,general_stack)
                  if isnothing(prev)
                        nb_stack=nb_stack+1
        		        tablex[i][j]=Int32(nb_stack)
        		        push!(general_stack,StackVect(nb_stack,PP(kb[pos].data),Int32(pos),Int32(0)))
                  else tablex[i][j]=Int32(prev) end
             end
          end
   end
   return(tablex,general_stack)
end

########################################
# vectors add-ons
########################################

function convert_coeffs_to_coeffs_z(coeffs::Vector{Vector{Rational{BigInt}}})
    coeffs_z=Vector{Vector{BigInt}}[]
   for c in coeffs
	 push!(coeffs_z,map(numerator, (c .* lcm(map(denominator,c)))))
   end
   return coeffs_z
end

function reduce_mod_p(
    cfs_zz::Vector{Vector{BigInt}},
    prime::UInt32
)
    @assert prime < typemax(UInt32)
    prime_zz = BigInt(prime)
    cfs_zp=[Vector{UInt32}(undef, length(cfs_zz[i])) for i in 1:length(cfs_zz)]
    @inbounds for i in 1:length(cfs_zz)
        cfs_zz_i = cfs_zz[i]
        for j in 1:length(cfs_zz_i)
            cfs_zp[i][j] = UInt32(mod(cfs_zz_i[j], prime_zz))
        end
        # if lead vanishes
        iszero(cfs_zp[i][1]) && return false, cfs_zp
    end
    true, cfs_zp
end

@inline function add_mul!(vres::Vector{UInt64},a::UInt32,v::Vector{UInt32})
    @fastmath @inbounds @simd ivdep for i in eachindex(vres)
        vres[i]+=UInt64(a)*UInt64(v[i])
    end
end

function reduce_mod!(vres::Vector{UInt64},pr::UInt32,arithm)
   @fastmath @inbounds @simd ivdep for i in eachindex(vres)
      vres[i]=Groebner.mod_p(vres[i],arithm)
   end
end

function vectorize_pol_gro!(p::PolUInt32,
                           kb::Vector{PP},
                           pr::UInt32,
                           arithm,
                           res::Vector{UInt32})
   pos=length(kb)
   @inbounds for i in 2:length(p.co)
      m=monomial(p,Int32(i))
      while ((pos>0) && (kb[pos].data!=m.data)) pos=pos-1 end
      if (pos<1) rur_print("\nError vectorisation ",pos," ",m)
      else
         res[pos]=pr-coeff(p,Int32(i))
      end
   end
   return res
end

function compute_fill_quo_gb!(t_xw::Vector{StackVect},
                              gro::Vector{PolUInt32},
                              quo::Vector{PP},
                              pr::UInt32,
                              arithm)
   t_v=[Vector{UInt32}() for i=1:length(t_xw)]
   @inbounds for i in eachindex(t_xw)
      if ((t_xw[i].var>0) && (t_xw[i].prev==0))
         # var>0 and prev=0 => groebner basis element at position var in the list gro
         t_v[i]=zeros(UInt32,length(quo))
         vectorize_pol_gro!(gro[t_xw[i].var],quo,pr,arithm,t_v[i])
      elseif ((t_xw[i].var==0) && (t_xw[i].prev>0))
         # var=0 (and prev>0 )=> quotient element at position prev in the list quo
          t_v[i]=[t_xw[i].prev]
        end
   end
   return t_v
end

function _mul_var_quo_UInt32!(v::Vector{UInt32},
                        ii::Int32,
                        t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        pr::UInt32,
                        arithm,
                        vv::Vector{UInt64},       # buffer
                        vres::Vector{UInt32})     # result
    dim=length(v); @assert dim==length(vres)
    pack=Int32(2^(floor(63-2*log(pr-1)/log(2))))
    resize!(vv, dim)
    @inbounds for i in 1:dim
        vv[i] = zero(UInt64)
    end
    continuer=true
    @inbounds for j=1:dim
        iszero(v[j]) && continue
        if (length(t_v[i_xw[ii][j]])>1)
            add_mul!(vv,v[j],t_v[i_xw[ii][j]])
            if j%pack==0
                reduce_mod!(vv,pr,arithm)
            end
        elseif (length(t_v[i_xw[ii][j]])==1)
            kk=t_v[i_xw[ii][j]][1]
            vv[kk]=Groebner.mod_p(vv[kk]+UInt64(v[j]), arithm)
        else
            continuer=false
            break;
        end
    end
    return continuer,reduce_mod!(vv,pr,arithm,vres)
end

function mul_var_quo_UInt32!(v::Vector{UInt32},
                        ii::Int32,
                        t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        pr::UInt32,
                        arithm,
                        vv::Vector{UInt64},
                        vres::Vector{UInt32})
    f,r= _mul_var_quo_UInt32!(v,ii,t_v,i_xw,pr,arithm,vv,vres)
    return r
end

function mul_var_quo_UInt32(v::Vector{UInt32},
                        ii::Int32,
                        t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        pr::UInt32,
                        arithm,
                        vv::Vector{UInt64})
    vres=Vector{UInt32}(undef,length(v))
    f,r= _mul_var_quo_UInt32!(v,ii,t_v,i_xw,pr,arithm,vv,vres)
    return r
end

function learn_compute_table!(t_v::Vector{Vector{UInt32}},
                              t_xw::Vector{StackVect},
                              i_xw::Vector{Vector{Int32}},
                              quo::Vector{PP},
                              pr::UInt32,
                              arithm)
   nb=1
   t_learn=Int32[]
   buf = Vector{UInt64}()
   while (nb>0)
      nb=0
      for i in eachindex(t_xw)
         # test if already computed
         if (length(t_v[i])==0)
            continuer=false
            # test if the ancestor is computed
            if (length(t_v[t_xw[i].prev])>0)
                vv=Vector{UInt32}(undef,length(t_v[t_xw[i].prev]))
                continuer,vv= _mul_var_quo_UInt32!(t_v[t_xw[i].prev],t_xw[i].var,t_v,i_xw,pr,arithm,buf,vv)
                if (continuer) 
                    t_v[i]=vv
                    push!(t_learn,i)
                    nb=nb+1
                end
            end
         end
      end
   end
   return(t_learn)
end

function learn_compute_table_sl!(t_v::Vector{Vector{UInt32}},
                                 t_xw::Vector{StackVect},
                                 i_xw::Vector{Vector{Int32}},
                                 quo::Vector{PP},
                                 pr::UInt32,
                                 arithm)

   nb=1
   t_learn=Int32[]
   buf = Vector{UInt64}()
   to_compute=copy(i_xw[end])
   for i in 1:(length(i_xw)-1)
        push!(to_compute,i_xw[i][1])
   end
   while (nb>0)
      nb=0
      for i in eachindex(to_compute)
         # test if already computed
         if (length(t_v[i])==0)
            continuer=false
            # test if the ancestor is computed
            if (length(t_v[t_xw[i].prev])>0)
                vv=Vector{UInt32}(undef,length(t_v[t_xw[i].prev]))
                continuer,vv= _mul_var_quo_UInt32!(t_v[t_xw[i].prev],t_xw[i].var,t_v,i_xw,pr,arithm,buf,vv)
                if (continuer) 
                    t_v[i]=vv
                    push!(t_learn,i)
                    nb=nb+1
                end
            end
         end
      end
   end
   return(t_learn)
end

function apply_compute_table!(t_v::Vector{Vector{UInt32}},
                              t_learn::Vector{Int32},
                              t_xw::Vector{StackVect},
                              i_xw::Vector{Vector{Int32}},
                              quo::Vector{PP},
                              pr::UInt32,
                              arithm)
    buf = Vector{UInt64}()
    @inbounds for i in t_learn
      if (length(t_v[i])==0)
          t_v[i]=mul_var_quo_UInt32(t_v[t_xw[i].prev],t_xw[i].var,t_v,i_xw,pr,arithm,buf)
      end
   end
   return nothing
end

function reduce_mod(v::Vector{UInt64},pr::UInt32,arithm)
   vres=Vector{UInt32}(undef, length(v))
    @fastmath @inbounds @simd for i in eachindex(vres)
      vres[i]=Groebner.mod_p(v[i], arithm)%UInt32
   end
   return vres
end

function reduce_mod!(v::Vector{UInt64},pr::UInt32,arithm,vres::Vector{UInt32})
     @fastmath @inbounds @simd for i in eachindex(vres)
       vres[i]=Groebner.mod_p(v[i], arithm)%UInt32
    end
    return vres
 end

############################################
# learn_zdim
############################################

@noinline __throw_not_generic(pr) = throw("Basis modulo $pr may be not generic enough. ")

function learn_zdim_quo(sys::Vector{AbstractAlgebra.Generic.MPoly{BigInt}},pr::UInt32,arithm,linform,cyclic,dd)
   sys_Int32=convert_to_mpol_UInt32(sys,pr)
   # Sasha: groebner_learn can use multi-threading itself (but let's not use it for now).
   graph,gro=groebner_learn_linform(sys_Int32,linform);
   for pr2 in Primes.nextprimes(UInt32(2^30), 2)
    sys_Int32_2=convert_to_mpol_UInt32(sys,pr2)
    # Sasha: groebner can use multi-threading (-//-)
    gb_2 = Groebner.groebner(sys_Int32_2,threaded=:no)
    length(gro) != length(gb_2) && __throw_not_generic(pr)
    for j in 1:length(gro)
        length(gro[j]) != length(gb_2[j]) && __throw_not_generic(pr)
    end
   end
   quo,g=kbase_linform(gro,pr,linform);
   gb_expvecs=map(poly->poly.exp,g)
   ltg=map(u->u.exp[1],g);
   q=map(u->u.exp[1],sys_mod_p(map(u->AbstractAlgebra.leading_monomial(u),quo),pr));
   i_xw,t_xw=prepare_table_mxi(ltg,q);
   t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
   if (cyclic)
      t_learn=learn_compute_table_sl!(t_v,t_xw,i_xw,q,pr,arithm);
      success,zp_param=zdim_parameterization(t_v,i_xw,pr,Int32(dd),Int32(1),arithm)
      if (success) 
        rur_print("\nApply cyclic optimization ")
      else
        rur_print("\nSwitch off cyclic optimization ")
        i_xw,t_xw=prepare_table_mxi(ltg,q);
        t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
        t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm);
      end
   else
      t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm);
   end
   return(graph,t_learn,t_v,q,i_xw,t_xw,pr,gb_expvecs)
end

function apply_zdim_quo!(graph,
                         t_learn::Vector{Int32},
                         q::Vector{PP},
                         i_xw::Vector{Vector{Int32}},
                         t_xw::Vector{StackVect},
                         pr::UInt32,
                         arithm,
                         gb_expvecs::Vector{Vector{PP}},
                         cfs_zp::Vector{Vector{UInt32}},
                         sys,
                         linform::Bool)
    success,gro=Groebner.groebner_applyX!(graph,cfs_zp,pr);
    if (success)
        g = [PolUInt32(gb_expvecs[i],gro[i]) for i in 1:length(gb_expvecs)]  # :^)
        t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
        apply_compute_table!(t_v,t_learn,t_xw,i_xw,q,pr,arithm);
        return success,t_v
    else
        rur_print("\n*** Bad prime detected in apply_zdim ***\n")
        return success,nothing
    end
end

function zdim_mx(t_v::Vector{Vector{Int32}},
                 i_xw::Vector{Vector{Int32}},
                 i::Int32,d::Int32)
    res=Vector{Int32}[]
    for j=1:d
        push!(res,t_v[i_xw[i][j]])
    end
    return(res)
end

############################################
# Rational Parameterization
############################################
# first non null element of a reduced row stored at position i in the gred is at
# position i+1 to simulate that gred[0]is not stored this element is the
# opposite element of the true pivot rows are normalized

function gauss_reduct(
        v::Vector{UInt32},
        gred::Vector{Vector{UInt32}},
        dg::Int32,
        dv::Int32,
        paquet::Int32,
        pr::UInt32,
        arithm,
        b::Vector{UInt64})   # buffer
    i=Int32(0)
    j=Int32(0)
    last_nn=Int32(-1)
    pivot=UInt32(0)
    resize!(b, dv)
    @inbounds for ell in 1:dv
        b[ell] = UInt64(v[ell])
    end
    @inbounds for i in 1:(dg-1)
        b[i+1]=Groebner.mod_p(b[i+1], arithm)
        iszero(b[i+1]) && continue
        if (length(gred[i])==0)
            last_nn=i%Int32
            break
        end
        pivot=(UInt64(pr)-b[i+1])%UInt32;
        add_mul!(b,pivot,gred[i]);
        b[i+1]=pivot;
        if (j<paquet)
            j+=Int32(1);
        else
            reduce_mod!(b,pr,arithm);
            j=Int32(0);
        end
    end
    if (j>0) 
        @inbounds for k in 1:dv v[k]=Groebner.mod_p(b[k], arithm)%UInt32; end
    else 
        @inbounds for k in 1:dv v[k]=b[k]%UInt32; end; 
    end
    if (last_nn==-1)
        @inbounds for ii in dg:(dv-1)
            if (v[ii+1]!=0)
                last_nn=ii%Int32;
                break;
            end
        end
    end
    if last_nn==-1
        return dv
    else
        return last_nn
    end
end

function make_monic_row!(v::Vector{UInt32},tmp::UInt32,pr::UInt32,arithm)
     @fastmath @inbounds @simd for i in eachindex(v)
        v[i]=Groebner.mod_p((v[i]%UInt64)*(tmp%UInt64), arithm)%UInt32
     end
end

function first_variable(t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        ii::Int32,pr::UInt32,
                        arithm)
    pack=Int32(2^(floor(63-2*log(pr-1)/log(2))))
    d=Int32(length(i_xw[ii])) 
    free_set=[append!([UInt32(1)],[UInt32(0) for i=2:d])]
    if (length(t_v[i_xw[ii][1]])>1)
        v=t_v[i_xw[ii][1]]
    else
        v = zeros(UInt32, d)
        v[t_v[i_xw[ii][1]][1]]=UInt32(1)
    end
    push!(free_set,v)
    # index[deg]=pos (pos=position in the gred with the convention that 0 is not stored)
    index=[Int32(i) for i=1:d];
    # normalization values
    hom=[UInt32(1) for i=1:d];
    gred=[Vector{UInt32}() for i=1:d];

    i=Int32(2);continuer=1;dg=Int32(1);new_i=Int32(1);deg=Int32(0);
    while ((i<d) && (v[i]==0)) i+=1; end
    gred[i-1]=Vector{UInt32}(v)
    if (gred[i-1][1]!=0) gred[i-1][1]=pr-gred[i-1][1]; end
    if (gred[i-1][i]!=1)
        hom[i-1]=invmod(gred[i-1][i],pr)%UInt32
        make_monic_row!(gred[i-1],hom[i-1],pr,arithm)
    end
    # T^0 is stored at gred[0] virtually
    index[1]=i-1;
    # gred[dg+i] not used i=0..d-dg
    dg=Int32(i);
    # we are now going to compute T^2
    deg=Int32(2);
    w=Vector{UInt32}(undef,d)
    buf1=Vector{UInt64}(undef,d)
    buf2=Vector{UInt64}(undef,d)
    @inbounds while (continuer==1)
        v=mul_var_quo_UInt32(v,ii,t_v,i_xw,pr,arithm,buf1)
        resize!(w, length(v))
        for ell in 1:length(v)
            w[ell] = v[ell]
        end
        # reduce with gred[0]
        if (w[1]!=0) w[1]=pr-w[1]; end
        # reduce with gred[i] i=1..dg-1
        # new_i = possible position to store the reduced vector in gred
        new_i=gauss_reduct(w,gred,dg,d,pack,pr,arithm,buf2);
        if (new_i<d)
            push!(free_set,v)
            gred[new_i]=copy(w);
            if (!(new_i<dg)) dg=Int32(new_i+1); end;
            index[deg]=Int32(new_i);
            hom[new_i]=invmod(gred[new_i][new_i+1],pr)%UInt32
            make_monic_row!(gred[new_i],hom[new_i],pr,arithm)
            deg+=Int32(1);
        else
            continuer=0;
        end;
    end
    # set v[i] cofficient of T^(i-1) in the min poly
    v=Vector{UInt32}(undef,deg)
    v[1]=w[1];
    @inbounds for i in 2:deg v[i]=Groebner.mod_p((w[index[i-1]+1]%UInt64)*(hom[index[i-1]]%UInt64), arithm)%UInt32 ; end;
    return(v,gred,index,dg,hom,free_set)
end

function biv_lex(t_v::Vector{Vector{UInt32}},
                 i_xw::Vector{Vector{Int32}},
                 gred::Vector{Vector{UInt32}},
                 index::Vector{Int32},
                 dg::Int32,
                 hom::Vector{UInt32},
                 free_set::Vector{Vector{UInt32}},
                 ii::Int32,
                 pr::UInt32,
                 arithm)
    pack=Int32(2^(floor(63-2*log(pr-1)/log(2))))
    d=Int32(length(i_xw[length(i_xw)])) 
    new_free_set=copy(free_set);
    deg=length(free_set);
    new_generators=Vector{Vector{Int32}}()
    new_monomial_free_set=[(Int32(i-1),Int32(0)) for i=1:deg]
    new_monomial_basis=copy(new_monomial_free_set);
    new_leading_monomials=Vector{Tuple{Int32,Int32}}()
    buf1=Vector{UInt64}()
    buf2=Vector{UInt64}()
    while (length(new_free_set)>0)
        tmp_set=Vector{Vector{UInt32}}()
        tmp_mon_set=Vector{Tuple{Int32,Int32}}()
        @inbounds for j in eachindex(new_free_set)
            curr_mon=new_monomial_free_set[j]
            curr_mon=(curr_mon[1],curr_mon[2]+Int32(1))
            v=mul_var_quo_UInt32(new_free_set[j],ii,t_v,i_xw,pr,arithm,buf1) 
            w=Vector{UInt32}(v)
            if (w[1]!=0) w[1]=pr-w[1]; end
            new_i=gauss_reduct(w,gred,dg,d,pack,pr,arithm,buf2);
            if (new_i<d)
                push!(tmp_set,v)
                push!(tmp_mon_set,curr_mon);
                gred[new_i]=Vector{UInt32}(w);
                if (!(new_i<dg)) dg=Int32(new_i+1); end;
                index[deg]=Int32(new_i);
                hom[new_i]=invmod(gred[new_i][new_i+1],pr)%UInt32
                make_monic_row!(gred[new_i],hom[new_i],pr,arithm)
                deg+=Int32(1);
            else
                v=Vector{UInt32}(undef,deg)
                v[1]=w[1];
                @inbounds for i in 2:deg v[i]=Groebner.mod_p((w[index[i-1]+1]%UInt64)*(hom[index[i-1]]%UInt64), arithm)%UInt32 ; end;
                push!(new_generators,copy(v))
                push!(new_leading_monomials,curr_mon);
                break;
            end;
        end;
        new_free_set=copy(tmp_set)
        new_monomial_free_set=copy(tmp_mon_set)
        append!(new_monomial_basis,copy(tmp_mon_set))
    end
    return(new_monomial_basis,new_leading_monomials,new_generators);
end

function convert_biv_lex_2_biv_pol(n_g,m_b,lt_g,pr)
    l_base=Vector{Vector{Vector{Int32}}}()
    for kk in eachindex(n_g)
         pp=n_g[kk]
         p_mat=Vector{Vector{UInt32}}()
         ldeg2=Int32(0)
         lco=[UInt32(0) for i=1:length(m_b)]
         for i in eachindex(pp)
                deg2=m_b[i][2]
                if deg2>ldeg2
                    push!(p_mat,lco)
                    for j in (ldeg2+1):(deg2-1)
                        push!(p_math,Vector{UInt32}())
                    end
                    lco=[UInt32(0) for j=1:length(m_b)]
                    ldeg2=deg2
                end
                lco[m_b[i][1]+1]=pp[i]
         end
         deg2=lt_g[kk][2]
         if (deg2==ldeg2) 
            lco[lt_g[kk][1]+1]=UInt32(1)
         else 
            push!(p_mat,lco);
            for i in (ldeg2+1):(deg2-1)
                        push!(p_math,Vector{UInt32}())
            end
            lco=[UInt32(0) for i=1:length(m_b)]
            lco[lt_g[kk][1]+1]=UInt32(1)
         end
         push!(p_mat,lco);
         push!(l_base,p_mat)
    end
    return(l_base)
end

function coeff_mod_p(x::Nemo.fpFieldElem,pr::UInt32)
    return(UInt32(Nemo.data(x)))
end

function check_separation_biv(bli,f,C)
    k=length(bli)-1;
    akk=C(bli[k+1])
    invakk=Nemo.invmod(k*akk,f)
    b=C(bli[k])*invakk
    for i in (k-1):-1:1
        tmp=Nemo.mod((k-i+1)*C(bli[i])-(i)*b*C(bli[i+1]),f)
        if (tmp!=0) 
            return(false);
        end
    end
    return(true);
end

function zdim_parameterization(t_v::Vector{Vector{UInt32}},
                               i_xw::Vector{Vector{Int32}},
                               pr::UInt32,dd::Int32,check::Int32,
                               arithm)
    res=Vector{Vector{UInt32}}()
    ii=Int32(length(i_xw))
    v,gred,index,dg,hom,free_set=first_variable(t_v,i_xw,ii,pr,arithm)
    C, _Z = Nemo.polynomial_ring(Nemo.Native.GF(Int64(pr),cached=false),cached=false)
    f=C(v)+_Z^(length(v));
    ifp=Nemo.derivative(f)
    f=f/Nemo.gcd(f,ifp)
    ifp=Nemo.derivative(f)
    push!(res,map(u->coeff_mod_p(u,pr),collect(Nemo.coefficients(f))));
    if (dd<0) flag=true 
    else  flag=(Int32(Nemo.degree(f))==dd)  end
    if (flag)    
        @inbounds for j in (ii-1):-1:1
            m_b,lt_b,n_g=biv_lex(t_v,i_xw,copy(gred),copy(index),copy(dg),copy(hom),copy(free_set),Int32(j),pr,arithm);
            bl=convert_biv_lex_2_biv_pol(n_g,m_b,lt_b,pr)
            s1=C([Int32(0)])
            s0=C([Int32(0)])
            pro=C([Int32(1)])
            ft=C(f)
            @inbounds for i in 1:length(bl)
                d1=length(bl[i])-1
                lc1=C(bl[i][d1+1])
                co0=C(bl[i][d1])
                f2=Nemo.gcd(ft,lc1)
                f1=ft/f2           
                # in some non radical case
                if (Nemo.degree(f1)>0)
                  if (check>0)
                        flag=check_separation_biv(bl[i],f1,C)
                        if (!flag) return(false,[],j) end
                  end
                  s1+=Nemo.mulmod(d1*lc1,pro,f)
                  s0+=Nemo.mulmod(co0,pro,f)
                  pro=pro*f1
                end
                ft=C(f2)
            end
            is1=Nemo.invmod(s1,f)
            s0=Nemo.mulmod(s0,is1,f)
            s0=-Nemo.mulmod(s0,ifp,f)
            push!(res,map(u->coeff_mod_p(u,pr),collect(Nemo.coefficients(s0))))
        end
    else 
        rur_print("Bad prime number for parameterization ",Int32(Nemo.degree(f)),"(",dd,")")
    end
    return(flag,res,ii)
end

function zz_mat_same_dims(zpm::Vector{Vector{UInt32}})::Vector{Vector{BigInt}}
    zzm=Vector{Vector{BigInt}}(undef,length(zpm))
    @inbounds for i in eachindex(zpm)
        zzm[i]=[BigInt(0) for j=1:length(zpm[i])]
    end
    return(zzm)
end

function qq_mat_same_dims(zpm::Vector{Vector{UInt32}})::Vector{Vector{Rational{BigInt}}}
    zzm=Vector{Vector{Rational{BigInt}}}(undef,length(zpm))
    @inbounds for i in eachindex(zpm)
        zzm[i]=[Rational{BigInt}(0) for j=1:length(zpm[i])]
    end
    return(zzm)
end

function param_modular_image(
        graph,
        t_learn,q,i_xw::Vector{Vector{Int32}},t_xw::Vector{StackVect},
        gb_expvecs,sys_z,cfs_zz,
        pr::UInt32,dd,linform::Bool)
    arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
    redflag,cfs_zp=reduce_mod_p(cfs_zz,pr)
    if !redflag
        rur_print("\n*** bad prime for lead detected ***\n")
        return false,nothing
    end
    success,t_v=apply_zdim_quo!(graph,t_learn,q,i_xw,t_xw,pr,arithm,gb_expvecs,cfs_zp,sys_z,linform);
    if !success
        rur_print("\n*** bad prime for Gbasis detected ***\n")
        return false,nothing
    end
    success,zp_param=zdim_parameterization(t_v,i_xw,pr,Int32(dd),Int32(0),arithm);
    if !success
        rur_print("\n*** bad prime for RUR detected ***\n")
        return false,nothing
    end
    return true,zp_param
end

function param_many_modular_images(
        graph,
        t_learn,q,i_xw::Vector{Vector{Int32}},t_xw::Vector{StackVect},
        gb_expvecs,sys_z,cfs_zz,
        prms::Vector{UInt32},dd,linform::Bool,
        prms_range,t_globl_success,t_globl_zp_params)
    for prm_id in prms_range
        pr=prms[prm_id]
        success,zp_param=param_modular_image(graph,t_learn,q,i_xw,t_xw,gb_expvecs,sys_z,cfs_zz,pr,dd,linform)
        !success && return nothing
        t_globl_success[prm_id]=success
        t_globl_zp_params[prm_id]=zp_param
    end
    return nothing
end

function general_param_serial(sys_z, nn, dd, linform::Bool,cyclic::Bool)::Vector{Vector{Rational{BigInt}}}
    qq_m=Vector{Vector{Rational{BigInt}}}()
    t_pr=Vector{UInt32}();
    t_param=Vector{Vector{Vector{UInt32}}}();
    pr=UInt32(Primes.prevprime(2^nn-1));
    arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
    rur_print("\nLearn ");
    graph,t_learn,t_v,q,i_xw,t_xw,pr,gb_expvecs=learn_zdim_quo(sys_z,pr,arithm,linform,cyclic,dd);
    backup=deepcopy(graph)
    expvecs,cfs_zz=extract_raw_data(sys_z)
    continuer=true
    bloc_p=Int32(2)
    lift_level=0
    kk=Int32(0)
    while(continuer)
        kk+=1
        pr=UInt32(Primes.prevprime(pr-1))
        success,zp_param=param_modular_image(graph,t_learn,q,i_xw,t_xw,gb_expvecs,sys_z,cfs_zz,pr,dd,linform)
        if !success
            kk-=1
            # The object may be corrupted after the failure. Revive it.
            graph=backup
            backup=deepcopy(backup)
            continue
        end
        length(t_param) > 0 && @assert map(length, t_param[end]) == map(length, zp_param)
        # From this point, we assume that the prime is OK
        push!(t_pr,UInt32(pr));
        push!(t_param,zp_param);
        kk < bloc_p && continue
        # Attempt reconstruction
        rur_print(length(t_pr));
        kk=0
        zz_p=BigInt(1);
        if (lift_level==0)
            rur_print("-");
            zz_m=zz_mat_same_dims([t_param[1][1]]);
            qq_m=qq_mat_same_dims([t_param[1][1]]);
            tt=[[t_param[ij][1]] for ij=1:length(t_pr)]
            Groebner.crt_vec_full!(zz_m,zz_p,tt,t_pr);
            aa=Groebner.ratrec_vec_full!(qq_m,zz_m,zz_p)
            if (aa) lift_level=1 end
        end
        if (lift_level==1)
            rur_print("+");
            zz_m=zz_mat_same_dims(t_param[1]);
            qq_m=qq_mat_same_dims(t_param[1]);
            Groebner.crt_vec_full!(zz_m,zz_p,t_param,t_pr);
            continuer=!Groebner.ratrec_vec_full!(qq_m,zz_m,zz_p);
        end
        bloc_p=Int32(max(floor(length(t_pr)/10),2))
    end;
    return qq_m;
end

function general_param(sys_z, nn, dd, linform::Bool,cyclic::Bool,parallelism::Symbol,threads,procs)
    if parallelism==:serial
        general_param_serial(sys_z,nn,dd,linform,cyclic)
    elseif parallelism==:multithreading
        # TODO(Sasha)
        # general_param_multithreading(sys_z,nn,dd,linform,cyclic,threads)
    else
        # general_param_multiprocessing(sys_z,nn,dd,linform,cyclic,procs)
    end
end

function swap_elts!(v,i,j)
    t=v[i]
    v[i]=v[j]
    v[j]=t
    return(v)
end

function isuniv(v)
    pos=0
    for i in eachindex(v)
        if (v[i]!=0) 
            if (pos>0) return(0)
            else pos=i end
        end
    end
    return(1)
end
function count_univ(l)
    nb=0
    for i in eachindex(l)
        nb+=l[i]
    end
    return(nb)
end

function prepare_system(ring, symbols, monoms, coeffs, nn, use_block)
    pr=UInt32(Primes.prevprime(2^nn-1));
    arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
    sys=(copy(monoms), copy(coeffs))
    
    flag, sys_Int32 = reduce_mod_p(coeffs,pr)
    @assert flag

    rur_print("\nTest the shape position")
    ring = Groebner.PolyRing(ring.nvars, ring.ord, pr)
    gro=Groebner.groebner(ring, monoms, sys_Int32, threaded=:no);
    if !(count_univ(map(w->isuniv(w),map(u -> u[1], gro[1]))) == ring.nvars)
        rur_print("\nError :  System is not zero-dimensional \n")
        return(-1,nothing,nothing,nothing,nothing,nothing)
    end

    # g=sys_mod_p(gro,pr);
    g = [PolUInt32(PP.(gro[1]), gro[2]) for i in 1:length(gro[1])]

    ltg=map(u->u.exp[1],g);
    # RU=parent(sys_Int32[1])
    q1=compute_quotient_basis(ltg);
    q=sort(q1,lt=pp_isless_drl)
    
    i_xw,t_xw=prepare_table_mxi(ltg,q);
    t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
    t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm); 
    ii=Int32(length(i_xw))
    
    flag,zp_param,uu=zdim_parameterization(t_v,i_xw,pr,Int32(-1),Int32(1),arithm);

    dd0=0
    if (flag) 
        dd0=length(zp_param[1])-1
        i_max=ii
        rur_print("\nSystem has a separating variable (",ii,")")
        return(dd0,length(q),(monoms,coeffs),false,dd0==length(q));
    end

    rur_print("\nTest if variables are separating")
    
    for ii in (length(i_xw)-1):-1:1
      ls=copy(symbols)
      ls=swap_elts!(ls,length(ls),ii)
      rur_print("(",ls[length(ls)],")")
    #   C,ls2=AbstractAlgebra.polynomial_ring(AbstractAlgebra.ZZ,push!(ls),internal_ordering=:degrevlex);
    #   lls=AbstractAlgebra.gens(C)
      monoms0 = map(u -> swap_elts!(u,length(ls),ii), monoms)
      # sys0=map(u->C(collect(AbstractAlgebra.coefficients(u)),map(u->swap_elts!(u,length(ls),ii),collect(AbstractAlgebra.exponent_vectors(u)))),sys_z);
      # sys_Int32=convert_to_mpol_UInt32(sys0,pr)
      flag, sys_Int32 = reduce_mod_p(coeffs,pr)

      gro=Groebner.groebner(ring, monoms0, sys_Int32,threaded=:no);

      g=sys_mod_p(gro,pr);

      ltg=map(u->u.exp[1],g);

      RU=parent(sys_Int32[1])
      q1=compute_quotient_basis(ltg);
      q=map(u->u.exp[1],sys_mod_p(sort(map(u->RU([1],[Vector{Int64}(u.data)]),q1)),pr))

        
      i_xw,t_xw=prepare_table_mxi(ltg,q);
      t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
      t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm); 
      ii=Int32(length(i_xw))
    
      flag,zp_param,uu=zdim_parameterization(t_v,i_xw,pr,Int32(-1),Int32(1),arithm);
    
      if (flag) 
        if ((length(zp_param[1])-1)>dd0)
            dd0=(length(zp_param[1])-1)
            i_max=ii
        end
        rur_print("\nSystem has a separating variable (",ls[length(ls)],")")
        return(dd0,length(q),sys0,AbstractAlgebra.symbols(C),false,dd0==length(q));
      end
    end

    rur_print("\nFind a separating vector")

    dd0=0
    
    ls=Vector{Symbol}(AbstractAlgebra.symbols(R))
    C,ls2=AbstractAlgebra.polynomial_ring(AbstractAlgebra.ZZ,push!(ls,:_Z),internal_ordering=:degrevlex);
    lls=AbstractAlgebra.gens(C)
    sys0=map(u->C(collect(AbstractAlgebra.coefficients(u)),map(u->push!(u,0),collect(AbstractAlgebra.exponent_vectors(u)))),sys_z);

    sep=[ BigInt(0) for i=1:length(lls)]
    sep[length(lls)]=1
    vv=length(lls)-1
    cc=BigInt(1)
    sep[length(lls)-1]=-BigInt(1)
    vv-=1
    dd=-1
    while (!flag) 
      sys=copy(sys0)
      lf=C(BigInt(0))
      for j in eachindex(lls)
            lf+=sep[j]*lls[j]
      end
      rur_print("\nTry ",lf)
      push!(sys,lf) 

      sys_Int32=convert_to_mpol_UInt32(sys,pr)
      linform=false 
      if (use_block)
        linform=true    # if a linear form was appended
        gro=groebner_linform(sys_Int32,linform)
        quo,g=kbase_linform(gro,pr,linform)
      else 
        gro=Groebner.groebner(sys_Int32,ordering=Groebner.DegRevLex(),threaded=:no);
        g=sys_mod_p(gro,pr);
      end
      ltg=map(u->u.exp[1],g);

    RU=parent(sys_Int32[1])
    q1=compute_quotient_basis(ltg);
    q=map(u->u.exp[1],sys_mod_p(sort(map(u->RU([1],[Vector{Int64}(u.data)]),q1)),pr))

      i_xw,t_xw=prepare_table_mxi(ltg,q);
      t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
      t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm); 
      ii=Int32(length(i_xw))
      flag,zp_param,uu=zdim_parameterization(t_v,i_xw,pr,Int32(-1),Int32(1),arithm);
      if (!flag)
        rur_print(" (",uu,",",vv,")")
        if (sep[uu]<0) sep[uu]=-sep[uu]
        else sep[uu]=-sep[uu]-1 end
        vv=uu
        if (vv<1)
            rur_print("\n Error in the choice of the separating form \n") 
            return(-1,-1,sys,AbstractAlgebra.symbols(C),linform)
        end
      else
        dd=length(zp_param[1])-1
      end  
    end

    if (dd==dd0)
         rur_print("\nSystem is in fact a separating variable (",ii,")")
        return(dd0,length(q),sys_z,AbstractAlgebra.symbols(R),false,dd==length(q))
    end    
    rur_print("\nSeparating form : ",sys[length(sys)],"\n")
    return(dd,length(q),sys,AbstractAlgebra.symbols(C),false,dd==length(q))
end

"""
    zdim_parameterization(sys; options...)
    zdim_parameterization(monoms, coeffs; options...)

Computes a RUR of a zero-dimensional system.

## Optional Arguments

- `nn`: an integer `<= 32`, use `nn`-bit prime numbers. 
    Default is `28`.
- `verbose`: a bool, whether to print progress statements or not. 
    Default is `true`.
- `parallelism`: a symbol, parallelism mode. 
    The options are `:multithreading`, `:multiprocessing`, and `:serial`.
    With `:multithreading`, the number of Julia threads can be specified via `threads` option.
    With `:multiprocessing`, the worker processes must be configured manually and specified via `procs` option.
    Default is `:serial`.
- `threads`: an integer, the number of threads. 
    Can only be used together with `parallelism=:multithreading`. 
    Default is `nthreads()`.
- `procs`: an array, the IDs of worker processes.
    Uses the process with `pid == 1` as the master process.
    Can only be used together with `parallelism=:multiprocessing`. 
    Default is not specified.

## Returns

- Returns the arrays of coefficients of `f(x), g_1(x), ..., g_n(x)`, 
    where `n` is the number of variables and `f(x)` is an annihilating polynomial of the separating element.
"""
function zdim_parameterization(
        sys;
        nn::Int32=Int32(28),
        use_block::Bool=false,
        verbose::Bool=true,
        parallelism=:serial,
        threads=nothing,
        procs=nothing)
    @assert AbstractAlgebra.base_ring(AbstractAlgebra.parent(sys[1])) == AbstractAlgebra.QQ
    @assert AbstractAlgebra.internal_ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
    monoms = map(f -> collect(AbstractAlgebra.exponent_vectors(f)), sys)
    coeffs = map(f -> collect(AbstractAlgebra.coefficients(f)), sys)
    ring = Groebner.PolyRing(
        AbstractAlgebra.nvars(AbstractAlgebra.parent(sys[1]))
        Groebner.DegRevLex(),
        0
    )
    symbols = AbstractAlgebra.symbols(AbstractAlgebra.parent(sys[1]))
    _zdim_parameterization(ring, symbols, monoms, coeffs, nn=nn, verbose=verbose)
end

function _zdim_parameterization(
        ring, symbols, monoms, coeffs; 
        nn::Int32=Int32(28),
        use_block::Bool=false,
        verbose::Bool=true,
        parallelism=:serial,
        threads=nothing,
        procs=nothing)
    _verbose[]=verbose
    @assert 1 <= nn <= 32 
    @assert parallelism in (:serial, :multithreading, :multiprocessing)
    if parallelism == :serial
        if !isnothing(threads) || !isnothing(procs)
            @warn "Options `threads` and `procs` are not supported in `:serial` mode and were ignored"
        end
        rur_print("Using serial mode.\n")
    end
    coeffs_z = convert_coeffs_to_coeffs_z(coeffs)
    rur_print("\nSeparation step");
    dm,Dq,sys_T,_vars,linform,cyclic=prepare_system(ring, symbols, monoms, coeffs_z,nn,use_block);
    if !(dm>0)
        rur_print("\nSomething went wrong");
        return []
    end
    rur_print("\nStart the computation (cyclic = ",cyclic,")");
    qq_m=general_param(sys_T,nn,dm,linform,cyclic,parallelism,threads,procs);
    rur_print("\n");
    return qq_m
end

function to_file(f_name,rur)
   cof=lcm(map(v->lcm(map(u->denominator(u),v)),rur));
   C, _Z = Nemo.polynomial_ring(Nemo.ZZ)
   rur_z=[C(map(u->Nemo.ZZ(numerator(u*cof)),rur[i])) for i=1:length(rur)];
   open(f_name, "w") do f
           write(f,"rur:=[",string(rur_z[1]))
           for i in 2:length(rur)
              write(f,",\n",string(rur_z[i]))
           end
           write(f,"]:\n")
       end
    return(nothing);
end

using PrecompileTools
include("precompile.jl")

end # module
