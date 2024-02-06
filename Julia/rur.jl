
import AbstractAlgebra,Nemo,Primes,Random

#####################################################
# Add on to AbstractAlgebra
#####################################################

function convert_sys_to_sys_z(sys::Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}})
   sys_z=AbstractAlgebra.Generic.MPoly{BigInt}[]
   for p in sys
	 push!(sys_z,AbstractAlgebra.change_coefficient_ring(AbstractAlgebra.ZZ,p*lcm(map(denominator,collect(AbstractAlgebra.coefficients(p))))))
   end
   return(sys_z)
end

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

#####################################################
# Add on to Groebner
#####################################################

const _enable_product_ord = Ref{Bool}(true)

# If linform == true and _enable_product_ord[] == true, then the basis of
#   [sys, f(_Z)]
# is computed by computing the basis of sys and appending f(_Z) to the result.

# We also tweak manually the order of monomials in f(_Z), since
# AbstractAlgebra.leading_monomial and friends operate in :degrevlex.

# If linform == true, then we assume that
# - the last element of the input array is the linear form,
# - the last variable in the AbstractAlgebra ring is the new variable.

function groebner_apply_linform!(graph,cfs_zp,pr,linform)
    if linform && _enable_product_ord[]
        success,gro=Groebner.groebner_applyX!(graph,cfs_zp[1:end-1],pr);
        gro=vcat(gro,[cfs_zp[end]])
        # Move the monomial with _Z up front as the largest
        _form=gro[end]
        _form=vcat(_form[end],_form[1:end-1])
        gro[end]=_form
        return success,gro
    else
        success,gro=Groebner.groebner_applyX!(graph,cfs_zp,pr);
        return success,gro
    end
end

function groebner_linform(sys_Int32,linform)
    if linform && _enable_product_ord[]
        gro=Groebner.groebner(sys_Int32[1:end-1],ordering=Groebner.DegRevLex());
        if any(f -> AbstractAlgebra.total_degree(AbstractAlgebra.leading_monomial(f)) == 1, gro)
            @warn "The basis is possibly not inter-reduced. Switching back to degrevlex."
            _enable_product_ord[] = false  # globally change the setting
            gro=Groebner.groebner(sys_Int32,ordering=Groebner.DegRevLex());
            return gro
        end
        gro=vcat(gro,sys_Int32[end])  
        return gro
    else
        gro=Groebner.groebner(sys_Int32,ordering=Groebner.DegRevLex());
        return gro
    end
end

function groebner_learn_linform(sys_Int32,linform)
    if linform && _enable_product_ord[]
        graph,gro=Groebner.groebner_learn(sys_Int32[1:end-1],ordering=Groebner.DegRevLex());
        @assert !any(f -> AbstractAlgebra.total_degree(AbstractAlgebra.leading_monomial(f)) == 1, gro)
        gro=vcat(gro,sys_Int32[end])  
        return graph,gro
    else
        graph,gro=Groebner.groebner_learn(sys_Int32,ordering=Groebner.DegRevLex());
        return graph,gro
    end
end

function kbase_linform(gro,pr,linform)
    if linform && _enable_product_ord[]
        R=AbstractAlgebra.parent(gro[1])
        vars=AbstractAlgebra.gens(R)
        @assert length(unique(map(string,vars))) == length(vars)   # variable names are unique
        ord=Groebner.Lex(vars[end])*Groebner.DegRevLex(vars[1:end-1])
        quo=Groebner.kbase(gro,ordering=ord);
        g=sys_mod_p(gro,pr)
        # Move the monomial with _Z up front as the largest
        _f=g[end]
        @assert length(_f.exp) == findfirst(e -> e.data[end] == 1, _f.exp)
        _f=PolUInt32(vcat(_f.exp[end], _f.exp[1:end-1]), vcat(_f.co[end], _f.co[1:end-1]))
        g[end] = _f
        return quo,g
    else
        quo=Groebner.kbase(gro,ordering=Groebner.DegRevLex());
        g=sys_mod_p(gro,pr);
        return quo,g
    end
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

## tests if pp1 is divisible by pp2
## then returns the greatest variable
## that appear in the divisor
@inline function divides(pp1::PP,pp2::PP)
    # tmp=pp1.data-pp2.data
    # i=1
    # var=0
    # while ((i<=length(tmp)) && (tmp[i]>=0))
    #     if (tmp[i]>0) var=i end
    #     i=i+1
    # end
    # return((i>length(tmp)),var)
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
# Interface with AbstractAlgebra
#####################################################
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
# Multiplication Matrices
#####################################################

#test if a monomial is already in the
#border and, if not, try to find a
#divisor of total degree -1
#the function supposes that it does not belong to
#the quotient
#variant where we look for the smallest divisor
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
            #test if mm is not in the quotient basis
            if (!((t_xw[pos].prev>0) && (t_xw[pos].var==0))) 
               flag,dd=divides(m,mm)
               if (flag)
                  res_flag=2 
                  res_dd=dd
                  res_pos=pos
               end
            end
      else
            #the total degree is to high to find a predecessor
            if (res_flag>0) 
               return(res_flag,res_dd,res_pos)
            else
               print("\nError find in border : ",tm," ",tmm)
               return(0,0,0)
            end
      end
      pos=pos-1
   end
   if (res_flag>0)
      return(res_flag,res_dd,res_pos)
   else 
      print("\n\nError in finding a predecessor")
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
	             else print("\n*** Error search table ***\n") end
	         else
		     #nm is a leading monomial of an element of the Gb
		     #we insert it with the flags prev=0 and var=pos in the gb
                     prev=findlast(item-> item.mon.data==nm.data,general_stack)
                     if isnothing(prev)
		                 nb_stack=nb_stack+1
		                 tablex[i][j]=Int32(nb_stack)
		                push!(general_stack,StackVect(nb_stack,PP(ltg[pos].data),Int32(0),Int32(pos)))
                     else tablex[i][j]=Int32(prev) end
		     end
	         else
		          #nm is a element of the quotient basis
		          #we insert it with the flags prev=pos and var=0
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

# stuff
########################################
function reduce_mod_p(
    cfs_zz::Vector{Vector{BigInt}},
    prime::UInt32
)
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

# vectors add-ons
########################################
@inline function add_mul!(vres::Vector{UInt64},a::UInt32,v::Vector{UInt32})
    @fastmath @inbounds @simd ivdep for i in eachindex(vres)
        vres[i]+=UInt64(a)*UInt64(v[i])
    end
end

@inline function add_mul_jam2!(
        vres::Vector{UInt64},
        a1::UInt32,v1::Vector{UInt32},
        a2::UInt32,v2::Vector{UInt32})
    @fastmath @inbounds @simd ivdep for i in 1:length(vres)
        s=vres[i]
        s=s+UInt64(a1)*UInt64(v1[i])
        s=s+UInt64(a2)*UInt64(v2[i])
        vres[i]=s
    end
end

@inline function add_mul_jam4!(
        vres::Vector{UInt64},
        a1::UInt32,v1::Vector{UInt32},
        a2::UInt32,v2::Vector{UInt32},
        a3::UInt32,v3::Vector{UInt32},
        a4::UInt32,v4::Vector{UInt32})
    @fastmath @inbounds @simd ivdep for i in 1:length(vres)
        s=vres[i]
        s=s+UInt64(a1)*UInt64(v1[i])
        s=s+UInt64(a2)*UInt64(v2[i])
        s=s+UInt64(a3)*UInt64(v3[i])
        s=s+UInt64(a4)*UInt64(v4[i])
        vres[i]=s
    end
end

function reduce_mod!(vres::Vector{UInt64},pr::UInt32,arithm)
   @fastmath @inbounds @simd ivdep for i in eachindex(vres)
      vres[i]=Groebner.mod_p(vres[i],arithm)
   end
end

function vectorize_pol_gro(p::PolUInt32,
                           kb::Vector{PP},
                           pr::UInt32,
                           arithm)
   res=zeros(UInt32, length(kb))
   pos=length(kb)
   @inbounds for i in 2:length(p.co)
      m=monomial(p,Int32(i))
      while ((pos>0) && (kb[pos].data!=m.data)) pos=pos-1 end
      if (pos<1) print("\nError vectorisation ",pos," ",m)
      else
         res[pos]=pr-coeff(p,Int32(i))
      end
   end
   return(res)
end

# elements of the quotient basis
# are represented by a vector of length 1
# whose entry is the place of the element in the quotient basis
function compute_fill_quo_gb!(t_xw::Vector{StackVect},
                              gro::Vector{PolUInt32},
                              quo::Vector{PP},
                              pr::UInt32,
                              arithm)
   t_v=[Vector{UInt32}() for i=1:length(t_xw)]
   @inbounds for i in eachindex(t_xw)
      #print("\n",t_xw[i])
      if ((t_xw[i].var>0) && (t_xw[i].prev==0))
         #var>0 and prev=0 => groebner basis element at position var in the list gro
         t_v[i]=vectorize_pol_gro(gro[t_xw[i].var],quo,pr,arithm)
      elseif ((t_xw[i].var==0) && (t_xw[i].prev>0))
         #var=0 (and prev>0 )=> quotient element at position prev in the list quo
          t_v[i]=[t_xw[i].prev]
      end
   end
   return(t_v)
end

function _mul_var_quo_UInt32(v::Vector{UInt32},
                        ii::Int32,
                        t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        pr::UInt32,
                        arithm,
                        vv::Vector{UInt64})  # buffer
    dim=length(v)
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
    return continuer,reduce_mod(vv,pr,arithm)
end

function mul_var_quo_UInt32(v::Vector{UInt32},
                        ii::Int32,
                        t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        pr::UInt32,
                        arithm,
                        vv::Vector{UInt64})
    f,r= _mul_var_quo_UInt32(v,ii,t_v,i_xw,pr,arithm,vv)
    return(r)
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
         #test if already computed
         if (length(t_v[i])==0)
            continuer=false
            #test if the ancestor is computed
            if (length(t_v[t_xw[i].prev])>0)
                continuer,vv= _mul_var_quo_UInt32(t_v[t_xw[i].prev],t_xw[i].var,t_v,i_xw,pr,arithm,buf)
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
end

function reduce_mod(v::Vector{UInt64},pr::UInt32,arithm)
   vres=Vector{UInt32}(undef, length(v))
    @fastmath @inbounds @simd for i in eachindex(vres)
      vres[i]=Groebner.mod_p(v[i], arithm)%UInt32
   end
   return(vres)
end

############################################
#learn_zdim
############################################

function learn_zdim_quo(sys::Vector{AbstractAlgebra.Generic.MPoly{BigInt}},pr::UInt32,arithm,linform)
#   pr=UInt32(Primes.prevprime(2^27-1));
#   @timeit tmr "convert" 
    sys_Int32=convert_to_mpol_UInt32(sys,pr)
# @timeit tmr "groebner learn"  
   # graph,gro=Groebner.groebner_learn(sys_Int32,ordering=Groebner.DegRevLex());
   graph,gro=groebner_learn_linform(sys_Int32,linform);
   # @timeit tmr "kbase" 
   # quo=Groebner.kbase(gro00,ordering=Groebner.DegRevLex());
   quo,g=kbase_linform(gro,pr,linform);

#    g=sys_mod_p(gro,pr);
   gb_expvecs=map(poly->poly.exp,g)
   ltg=map(u->u.exp[1],g);
   q=map(u->u.exp[1],sys_mod_p(map(u->AbstractAlgebra.leading_monomial(u),quo),pr));

   # @timeit tmr "prepare table" 
   i_xw,t_xw=prepare_table_mxi(ltg,q);
   t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
# @timeit tmr "learn table"  
   t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm); 
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
                         linform::Bool)
#   @timeit tmr "convert" 
    # sys_Int32=convert_to_mpol_UInt32(sys,pr);
#   @timeit tmr "groebner" 
    # f,gro=Groebner.groebner_apply!(graph,sys_Int32);
    success,gro=groebner_apply_linform!(graph,cfs_zp,pr,linform)
    #success,gro=Groebner.groebner_applyX!(graph,cfs_zp,pr);
    if (success)
#   @timeit tmr "gro mod p" 
    # g=sys_mod_p(gro,pr);
        g = [PolUInt32(gb_expvecs[i],gro[i]) for i in 1:length(gb_expvecs)]  # :^)
#   @timeit tmr "fill"
        t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
#  @timeit tmr "table" 
        apply_compute_table!(t_v,t_learn,t_xw,i_xw,q,pr,arithm);
        return(success,t_v)
    else
        print("\n*** Bad prime detected in apply_zdim ***\n")
        return(success,nothing)
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
#first non null element of a reduced row stored at position i 
#in the gred is at position i+1 to simulate that
#gred[0]is not stored
#this element is the opposite element of the true pivot
#rows are normalized

# Run "_gauss_reduct_jam[] = 1" to disable jamming
const _gauss_reduct_jam = Ref{Int}(2)

# @timeit tmr "gred" 
function gauss_reduct(
        v::Vector{UInt32},
        gred::Vector{Vector{UInt32}},
        dg::Int32,
        dv::Int32,
        paquet::Int32,
        pr::UInt32,
        arithm,
        b::Vector{UInt64})   # buffer
    if paquet > 3 && _gauss_reduct_jam[] == 4
        return gauss_reduct_jam4(v,gred,dg,dv,paquet,pr,arithm,b)
    elseif  paquet > 1 && _gauss_reduct_jam[] == 2
        return gauss_reduct_jam2(v,gred,dg,dv,paquet,pr,arithm,b)
    end
    i=Int32(0)
    j=Int32(0)
    last_nn=Int32(-1)
    pivot=UInt32(0)
    resize!(b, dv)
    @inbounds for ell in 1:dv
        b[ell] = UInt64(v[ell])
    end
    # b=[UInt64(v[i]) for i=1:dv]
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

function gauss_reduct_jam2(
        v::Vector{UInt32},
        gred::Vector{Vector{UInt32}},
        dg::Int32,
        dv::Int32,
        paquet::Int32,
        pr::UInt32,
        arithm,
        b::Vector{UInt64})
    @assert paquet > 1
    j=Int32(0)
    last_nn=Int32(-1)
    pivot=UInt32(0)
    resize!(b, dv)
    @inbounds for ell in 1:dv
        b[ell] = UInt64(v[ell])
    end
    # b=[UInt64(v[i]) for i=1:dv]
    ub2 = (dg-1) & xor(typemax(Int), 1)
    k = 1
    @inbounds while k <= ub2
        w1, w2 = gred[k], gred[k+1]
        b[k+1]=Groebner.mod_p(b[k+1], arithm)
        if isempty(w1)
            if !iszero(b[k+1])
                last_nn=k%Int32
                break
            end
            # Go for a single iteration with w2
            b[k+2]=Groebner.mod_p(b[k+2], arithm)
            if iszero(b[k+2])
                k += 2
                continue
            end
            if isempty(w2)
                last_nn=(k+1)%Int32
                break
            end
            pivot=(UInt64(pr)-b[k+2])%UInt32;
            add_mul!(b,pivot,w2);
            b[k+2]=pivot;
            if (j<paquet)
                j+=Int32(1);
            else
                reduce_mod!(b,pr,arithm);
                j=Int32(0);
            end
            k += 2
            continue
        end
        # from this point, w1 is non-empty
        if isempty(w2)
            # Go for a single iteration with w1
            pivot=(UInt64(pr)-b[k+1])%UInt32;
            add_mul!(b,pivot,w1);
            b[k+1]=pivot;
            if (j<paquet)
                j+=Int32(1);
            else
                reduce_mod!(b,pr,arithm);
                j=Int32(0);
            end
            b[k+2]=Groebner.mod_p(b[k+2], arithm)
            if !iszero(b[k+2])
                last_nn=(k+1)%Int32
                break
            end
            k+=2
            continue
        end
        # from this point, w1 and w2 are non-empty
        pivot1=(UInt64(pr)-b[k+1])%UInt32;
        b[k+2]=Groebner.mod_p(b[k+2], arithm)
        pivot2=(UInt64(pr)-Groebner.mod_p(b[k+2]+UInt64(pivot1)*w1[k+2], arithm))%UInt32;
        add_mul_jam2!(b, pivot1, w1, pivot2, w2) 
        b[k+1]=Groebner.mod_p(UInt64(pivot1)+UInt64(pivot2)*UInt64(w2[k+1]), arithm)%UInt32;
        b[k+2]=pivot2;
        if (j<paquet)
            j+=Int32(2);
        else
            reduce_mod!(b,pr,arithm);
            j=Int32(0);
        end
        k+=2
    end

    if last_nn == -1
        @inbounds for i in (k):(dg-1)
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

function gauss_reduct_jam4(
        v::Vector{UInt32},
        gred::Vector{Vector{UInt32}},
        dg::Int32,
        dv::Int32,
        paquet::Int32,
        pr::UInt32,
        arithm,
        b::Vector{UInt64})      # buffer
    @assert paquet > 3
    j=Int32(0)
    last_nn=Int32(-1)
    pivot=UInt32(0)
    resize!(b, dv)
    @inbounds for ell in 1:dv
        b[ell] = UInt64(v[ell])
    end
    # b=[UInt64(v[i]) for i=1:dv]
    # the largest number not greater than (dg-1) and divisible by 4
    ub4 = (dg-1) & xor(typemax(Int), 3)
    k = 1
    # k repeatedly increases by exactly 4
    @inbounds while k <= ub4
        k1,k2,k3,k4 = k,k+1,k+2,k+3
        w1,w2,w3,w4 = gred[k1],gred[k2],gred[k3],gred[k4]
        # If any of the reducers are empty, then default to scalar iteration
        if isempty(w1) || isempty(w2) || isempty(w3) || isempty(w4)
            for ki in k:(k+3)
                b[ki+1]=Groebner.mod_p(b[ki+1], arithm)
                iszero(b[ki+1]) && continue
                if (length(gred[ki])==0)
                    last_nn=ki%Int32;
                    break
                end
                pivot=Groebner.mod_p((UInt64(pr)-b[ki+1]),arithm)%UInt32;
                add_mul!(b,pivot,gred[ki]);
                b[ki+1]=pivot;
                if (j<paquet)
                    j+=Int32(1);
                else
                    reduce_mod!(b,pr,arithm);
                    j=Int32(0);
                end
            end
            (last_nn != -1) && break    # new pivot found, break out
            k=k+4
            continue
        end
        # From this point, w1, w2, w3, w4 are non-empty.
        # Jam!
        b[k1+1]=Groebner.mod_p(b[k1+1], arithm)
        b[k2+1]=Groebner.mod_p(b[k2+1], arithm)
        b[k3+1]=Groebner.mod_p(b[k3+1], arithm)
        b[k4+1]=Groebner.mod_p(b[k4+1], arithm)
        pivot1=(UInt64(pr)-b[k1+1])%UInt32;
        pivot2=(UInt64(pr)-Groebner.mod_p(
            b[k2+1] + UInt64(pivot1)*w1[k2+1], arithm))%UInt32;
        pivot3=(UInt64(pr)-Groebner.mod_p(
            b[k3+1] + UInt64(pivot1)*w1[k3+1] + UInt64(pivot2)*w2[k3+1], arithm))%UInt32;
        pivot4=(UInt64(pr)-Groebner.mod_p(
            b[k4+1] + UInt64(pivot1)*w1[k4+1] + UInt64(pivot2)*w2[k4+1] + UInt64(pivot3)*w3[k4+1], arithm))%UInt32;
        add_mul_jam4!(b, pivot1, w1, pivot2, w2, pivot3, w3, pivot4, w4)
        b[k1+1]=Groebner.mod_p(
            UInt64(pivot1) + UInt64(pivot4)*UInt64(w4[k1+1]) + UInt64(pivot3)*UInt64(w3[k1+1]) + UInt64(pivot2)*UInt64(w2[k1+1]), arithm)%UInt32;
        b[k2+1]=Groebner.mod_p(
            UInt64(pivot2) + UInt64(pivot4)*UInt64(w4[k2+1]) + UInt64(pivot3)*UInt64(w3[k2+1]), arithm)%UInt32
        b[k3+1]=Groebner.mod_p(
            UInt64(pivot3) + UInt64(pivot4)*UInt64(w4[k3+1]), arithm)%UInt32
        b[k4+1]=pivot4
        if (j<paquet)
            j+=Int32(4);
        else
            reduce_mod!(b,pr,arithm);
            j=Int32(0);
        end
        k=k+4
    end

    # the scalar tail. At most 3 iterations.
    if last_nn == -1
        @inbounds for i in (k):(dg-1)
            b[i+1]=Groebner.mod_p(b[i+1], arithm)
            iszero(b[i+1]) && continue
            if (length(gred[i])==0)
                last_nn=i%Int32;
                break
            end
            pivot=Groebner.mod_p((UInt64(pr)-b[i+1]),arithm)%UInt32;
            add_mul!(b,pivot,gred[i]);
            b[i+1]=pivot;
            if (j<paquet)
                j+=Int32(1);
            else
                reduce_mod!(b,pr,arithm);
                j=Int32(0);
            end
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

function normalize_row!(v::Vector{UInt32},tmp::UInt32,pr::UInt32,arithm)
     @fastmath @inbounds @simd for i in eachindex(v)
        v[i]=Groebner.mod_p((v[i]%UInt64)*(tmp%UInt64), arithm)%UInt32
     end
end

function first_variable(t_v::Vector{Vector{UInt32}},
                        i_xw::Vector{Vector{Int32}},
                        ii::Int32,pr::UInt32,
                        arithm)
    pack=Int32(2^(floor(63-2*log(pr-1)/log(2))))
    d=Int32(length(i_xw[1])) 
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
        normalize_row!(gred[i-1],hom[i-1],pr,arithm)
    end
    #T^0 is stored at gred[0] virtually
    index[1]=i-1;
    #gred[dg+i] not used i=0..d-dg
    dg=Int32(i);
    #we are now going to compute T^2
    deg=Int32(2);
    w=Vector{UInt32}(undef,d)
    buf1=Vector{UInt64}(undef,d)
    buf2=Vector{UInt64}(undef,d)
    @inbounds while (continuer==1)
        v=mul_var_quo_UInt32(v,ii,t_v,i_xw,pr,arithm,buf1)
        # w=Vector{UInt32}(v)
        resize!(w, length(v))
        for ell in 1:length(v)
            w[ell] = v[ell]
        end
        #reduce with gred[0]
        if (w[1]!=0) w[1]=pr-w[1]; end
        #reduce with gred[i] i=1..dg-1
        #new_i = possible position to store the reduced vector in gred
        new_i=gauss_reduct(w,gred,dg,d,pack,pr,arithm,buf2);
        if (new_i<d)
            push!(free_set,v)
            gred[new_i]=copy(w);
            if (!(new_i<dg)) dg=Int32(new_i+1); end;
            index[deg]=Int32(new_i);
            hom[new_i]=invmod(gred[new_i][new_i+1],pr)%UInt32
            normalize_row!(gred[new_i],hom[new_i],pr,arithm)
            deg+=Int32(1);
        else
            continuer=0;
        end;
    end
    #set v[i] cofficient of T^(i-1) in the min poly
    v=Vector{UInt32}(undef,deg)
    v[1]=w[1];
#old    @inbounds for i in 2:deg v[i]=w[index[i-1]+1]; end;
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
    d=Int32(length(i_xw[1])) 
    new_free_set=copy(free_set);
    deg=length(free_set);
    new_generators=Vector{Vector{Int32}}()
    new_monomial_free_set=[[Int32(i-1),Int32(0)] for i=1:deg]
    new_monomial_basis=copy(new_monomial_free_set);
    new_leading_monomials=Vector{Vector{Int32}}()
    buf1=Vector{UInt64}()
    buf2=Vector{UInt64}()
    while (length(new_free_set)>0)
        tmp_set=Vector{Vector{UInt32}}()
        tmp_mon_set=Vector{Vector{Int32}}()
        @inbounds for j in eachindex(new_free_set)
            curr_mon=copy(new_monomial_free_set[j])
            curr_mon[2]+=Int32(1)
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
                normalize_row!(gred[new_i],hom[new_i],pr,arithm)
                deg+=Int32(1);
            else
                v=Vector{UInt32}(undef,deg)
                v[1]=w[1];
                @inbounds for i in 2:deg v[i]=Groebner.mod_p((w[index[i-1]+1]%UInt64)*(hom[index[i-1]]%UInt64), arithm)%UInt32 ; end;
                push!(new_generators,copy(v))
                push!(new_leading_monomials,copy(curr_mon));
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
            print("\n Error - variable is not separated ",Nemo.degree(f)," ** ",tmp);
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
    C, _Z = Nemo.polynomial_ring(Nemo.Native.GF(Int64(pr)))
    f=C(v)+_Z^(length(v));
    ifp=Nemo.derivative(f)
    f=f/Nemo.gcd(f,ifp)
    ifp=Nemo.derivative(f)
    push!(res,map(u->coeff_mod_p(u,pr),collect(Nemo.coefficients(f))));
    if (dd<0) flag=true else flag=(Int32(Nemo.degree(f))==dd) end
    if (flag)    
        @inbounds for j in 1:(ii-1)
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
                        if (!flag) print("\nError - System not separated"); end
                  end
                  #lc1=Nemo.mod(-d1*lc1,f1)
                  s1+=Nemo.mulmod(d1*lc1,pro,f)
                  #co0=Nemo.mod(co0,f1)
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
        print("Bad prime number for parameterization ",Int32(Nemo.degree(f)),"(",dd,")")
    end
    return(flag,res)
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

function general_param(sys_z, nn, dd, linform::Bool)::Vector{Vector{Rational{BigInt}}}
    qq_m=Vector{Vector{Rational{BigInt}}}()
    t_pr=Vector{UInt32}();
    t_param=Vector{Vector{Vector{UInt32}}}();
    pr=UInt32(Primes.prevprime(2^nn-1));
    arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
    graph,t_learn,t_v,q,i_xw,t_xw,pr,gb_expvecs=learn_zdim_quo(sys_z,pr,arithm,linform);
    backup=deepcopy(graph)
    expvecs,cfs_zz=extract_raw_data(sys_z)
    continuer=true
    bloc_p=Int32(2)
    lift_level=0
    kk=Int32(0)
    while(continuer)
        kk+=1
        pr=UInt32(Primes.prevprime(pr-1))
        arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
        redflag,cfs_zp=reduce_mod_p(cfs_zz,pr)
        if !redflag
            print("\n*** bad prime for lead detected ***\n")
            continue
        end
        # success,t_v=nothing,nothing
        # try
        success,t_v=apply_zdim_quo!(graph,t_learn,q,i_xw,t_xw,pr,arithm,gb_expvecs,cfs_zp,linform);
        # catch err
        #     println("***** Fatal error. Debug info will follow *****")
        #     println(err)
        #     println("nn=$nn, dd=$dd, linform=$linform, pr=$pr")
        #     println("sys_z=\n$sys_z")
        #     println("cfs_zp=\n$cfs_zp")
        #     rethrow(err)
        # end
        if !success
            print("\n*** bad prime for Gbasis detected ***\n")
            # The object may be corrupted after the failure. Revive it.
            graph=backup
            backup=deepcopy(backup)
            continue
        end
        flag,zp_param=zdim_parameterization(t_v,i_xw,pr,Int32(dd),Int32(0),arithm);
        if !flag
            print("\n*** bad prime for RUR detected ***\n")
            continue
        end
        length(t_param) > 0 && @assert map(length, t_param[end]) == map(length, zp_param)
        # From this point, we assume that the prime is OK
        push!(t_pr,UInt32(pr));
        push!(t_param,zp_param);
        kk != bloc_p && continue
        # Attempt reconstruction
        print("[",kk,",",length(t_pr),"]");
        Base.flush(Base.stdout)
        kk=0
        zz_p=BigInt(1);
        if (lift_level==0)
            print("(",0,")");
            zz_m=zz_mat_same_dims([t_param[1][1]]);
            qq_m=qq_mat_same_dims([t_param[1][1]]);
            tt=[[t_param[ij][1]] for ij=1:length(t_pr)]
            Groebner.crt_vec_full!(zz_m,zz_p,tt,t_pr);
            aa=Groebner.ratrec_vec_full!(qq_m,zz_m,zz_p)
            if (aa) lift_level=1 end
        end
        if (lift_level==1)
            print("(",1,")");
            zz_m=zz_mat_same_dims(t_param[1]);
            qq_m=qq_mat_same_dims(t_param[1]);
            Groebner.crt_vec_full!(zz_m,zz_p,t_param,t_pr);
            continuer=!Groebner.ratrec_vec_full!(qq_m,zz_m,zz_p);
        end
        #bloc_p=Int32(2^(Int32(floor(max(log(length(t_pr)/10)/log(2),2)))))
        bloc_p=Int32(max(floor(length(t_pr)/10),2))
    end;
    return qq_m;
end

function prepare_system(sys_z, nn,R,use_block)
    ls=Vector{Symbol}(AbstractAlgebra.symbols(R))
    C,ls2=polynomial_ring(AbstractAlgebra.ZZ,push!(ls,:_Z),ordering=:degrevlex);
    #C,ls2=polynomial_ring(AbstractAlgebra.ZZ,push!(ls,:_Z),ordering=:degrevlex,cached=false);
    lls=AbstractAlgebra.gens(C)
    sys=map(u->C(collect(AbstractAlgebra.coefficients(u)),map(u->push!(u,0),collect(AbstractAlgebra.exponent_vectors(u)))),sys_z);

    lf=lls[length(lls)]
    sep=Vector{BigInt}()
    push!(sep,BigInt(1))
    lf-=lls[1]
    for i in 2:length(ls)
        cc=-BigInt(rand(-50:50))
        lf+=C(cc)*lls[i]
        sep=push!(sep,cc)
    end
    push!(sys,lf)
  
    pr=UInt32(Primes.prevprime(2^nn-1));
    arithm=Groebner.ArithmeticZp(UInt64, UInt32, pr)
    sys_Int32=convert_to_mpol_UInt32(sys,pr)
    linform=false 
    if (use_block)
        linform=true    # if a linear form was appended
        gro=groebner_linform(sys_Int32,linform)
        quo,g=kbase_linform(gro,pr,linform)
    else 
        gro=Groebner.groebner(sys_Int32,ordering=Groebner.DegRevLex());
        quo=Groebner.kbase(gro,ordering=Groebner.DegRevLex());
        g=sys_mod_p(gro,pr);
    end
    ltg=map(u->u.exp[1],g);
    q=map(u->u.exp[1],sys_mod_p(map(u->AbstractAlgebra.leading_monomial(u),quo),pr));
    i_xw,t_xw=prepare_table_mxi(ltg,q);
    t_v=compute_fill_quo_gb!(t_xw,g,q,pr,arithm);
    t_learn=learn_compute_table!(t_v,t_xw,i_xw,q,pr,arithm); 
    ii=Int32(length(i_xw))
    
    flag,zp_param=zdim_parameterization(t_v,i_xw,pr,Int32(-1),Int32(1),arithm);

    if (!flag) 
        print("\n Error in the choice of the separating form \n") 
    end
    dd=length(zp_param[1])
    U, _Z = Nemo.polynomial_ring(Nemo.Native.GF(Int64(pr)))
    
    f=U(zp_param[1]);

    dd=Nemo.degree(f)
    ii-=1
    while (ii>0)
        v,gred,index,dg,hom=first_variable(t_v,i_xw,Int32(ii),pr,arithm)
        f=U(v)+_Z^(length(v));
        ifp=Nemo.derivative(f)
        f=f/Nemo.gcd(f,ifp)
        if (dd==Nemo.degree(f)) 
            print("\nvariable ",ii," is separating ",Nemo.degree(f),"(",length(q),")\n")
            break;
        end
        ii-=1
    end
    if (ii>0)
        ls3=Vector{Symbol}()
        for i in 1:(ii-1) push!(ls3,ls[i]) end
        for i in (ii+1):(length(ls)-1) push!(ls3,ls[i]) end
        push!(ls3,ls[ii])
        C,ls=polynomial_ring(AbstractAlgebra.ZZ,ls3,ordering=:degrevlex)
        sys=map(u->C(collect(AbstractAlgebra.coefficients(u)),collect(AbstractAlgebra.exponent_vectors(u))),sys_z);
        linform=false
    else 
        print("\nSeparating vector ",dd,"(",length(q),")\n")
    end
    return(dd,length(q),sys,AbstractAlgebra.symbols(C),linform)
end

function zdim_parameterization(sys,nn::Int32=Int32(28),use_block::Bool=false)
    @assert AbstractAlgebra.ordering(AbstractAlgebra.parent(sys[1])) == :degrevlex
    sys_z=convert_sys_to_sys_z(sys);
    dm,Dq,sys_T,_vars,linform=prepare_system(sys_z,nn,AbstractAlgebra.parent(sys[1]),use_block);
    qq_m=general_param(sys_T,nn,dm,linform);
    return(qq_m)
end
