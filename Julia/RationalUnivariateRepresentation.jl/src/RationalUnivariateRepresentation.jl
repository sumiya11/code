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

import Groebner
import Nemo
import Primes
import TimerOutputs
import AbstractAlgebra
import Base.Threads: nthreads, @threads, threadid

const to = TimerOutputs.TimerOutput()

# ************************************************************************
# Parameters of the packages
# ************************************************************************

const pr_max_bitsize = 28

const ModularCoeff = UInt32
const AccModularCoeff = UInt64
const CNCoeff{N,C} = Groebner.CompositeNumber{N,C} where {N,C}

const ModularArithmetic{Accum,Coeff} = Groebner.AbstractArithmeticZp{Accum,Coeff} where {Accum,Coeff}

const Deg = UInt32
const PP = Vector{Deg}

const ModularArithZp = Groebner.ArithmeticZp
const ModularArithMZp = Groebner.CompositeArithmeticZp
const ModularImage = Groebner.mod_p
const ModularPrime = Groebner.divisor
const ModularInverse = Groebner.inv_mod_p

mutable struct StackVect
    pos::Int32
    mon::PP
    prev::Int32
    var::Int32
end

# ************************************************************************
# Usefull functions
# ************************************************************************

const _verbose = Ref{Bool}(true)

rur_print(xs...) = rur_print(stdout, xs...)
function rur_print(io::IO, xs...)
    _verbose[] && print(io, xs...)
    Base.flush(io)
end

bitsize(a::Integer) = Int32(ceil(log(a) / log(2)))

#log[2](2^64/(2*(pr-1)^2))=log[2](2^63/((pr-1)^2))=64-1-2*log[2](pr-1)
function pack_value(arithm)
    Int32(2^(8 * sizeof(AccModularCoeff) - 1 - 2 * bitsize((2^pr_max_bitsize - 1) - 1)))
end

# ************************************************************************
# Interface with Primes
# ************************************************************************

const PrevPrime = Primes.prevprime
const PrevPrimes = Primes.prevprimes

# ************************************************************************
# Interface with AbstractAlgebra
# ************************************************************************

function coeff_mod_p(x)
    return (ModularCoeff(AbstractAlgebra.data(x)))
end

# ************************************************************************
# Groebner.jl interface
# ************************************************************************

TimerOutputs.@timeit to "Groebner Classic" function _gb_4_rur(de, co, arithm)
    ri = Groebner.PolyRing(length(de[1][1]), Groebner.DegRevLex(), Groebner.divisor(arithm))
    gro_exp, gro_coeffs = Groebner.groebner(ri, de, co, ordering = Groebner.DegRevLex(), threaded = :no)
    return (gro_exp, gro_coeffs)
end

TimerOutputs.@timeit to "Groebner learn" function _gb_4_rur_learn(de, co, arithm)
    ri = Groebner.PolyRing(length(de[1][1]), Groebner.DegRevLex(), Groebner.divisor(arithm))
    graph, gro_exp, gro_coeffs = Groebner.groebner_learn(ri, de, co, loglevel=:warn, ordering = Groebner.DegRevLex(), threaded = :no)
    return (graph, gro_exp, gro_coeffs)
end

TimerOutputs.@timeit to "Groebner Apply" function _gb_4_rur_apply!(graph, de, co, arithm)
    if arithm isa Vector  # apply with composite numbers
        ri = map(a -> Groebner.PolyRing(length(de[1][1]), Groebner.DegRevLex(), Groebner.divisor(a)), arithm)
        batch = tuple(map(ri_co -> (ri_co[1], de, ri_co[2]), zip(ri, co))...)
        success, gro_coeffs = Groebner.groebner_apply!(graph, batch, loglevel=:warn, ordering = Groebner.DegRevLex(), threaded = :no)
        return (success, gro_coeffs)
    else  # apply scalar
        ri = Groebner.PolyRing(length(de[1][1]), Groebner.DegRevLex(), Groebner.divisor(arithm))
        success, gro_coeffs = Groebner.groebner_apply!(graph, ri, de, co, loglevel=:warn, ordering = Groebner.DegRevLex(), threaded = :no)
        return (success, gro_coeffs)
    end
end

# ************************************************************************
# Extension for Power products
# ************************************************************************

total_degree(pp::PP) = sum(pp)

function pp_isless_drl(ea::PP, eb::PP)
    @assert length(ea) == length(eb)
    if sum(ea) < sum(eb)
        return true
    elseif sum(ea) != sum(eb)
        return false
    end
    i = length(ea)
    @inbounds while i > 1 && ea[i] == eb[i]
        i -= 1
    end
    @inbounds return ea[i] > eb[i]
end

function divides(a::PP, b::PP)
    var = 0
    @inbounds for j = 1:length(a)
        if a[j] > b[j]
            var = j
        end
        if a[j] < b[j]
            return false, var
        end
    end
    return true, var
end

function find_divisor(m::PP, lm::Vector{PP})
    pos = length(lm)
    @inbounds for j = pos:-1:1
        te, v = divides(m, lm[j])
        if te
            return (j)
        end
    end
    return (0)
end

function find_in_list(m::PP, lm::Vector{PP})
    pos = length(lm)
    @inbounds for j = pos:-1:1
        if m == lm[j]
            return (j)
        end
    end
    return (0)
end

function compute_quotient_basis(ltg::Vector{PP})
    if length(ltg) == 0
        return ([])
    end
    nbv = length(ltg[1])
    quo = Vector{PP}()
    todo = Vector{PP}([PP([Deg(0) for i = 1:nbv])])
    inspected = Vector{PP}([PP([Deg(0) for i = 1:nbv])])
    pos = 1
    while (length(todo) > 0)
        m = popfirst!(todo)
        pos = find_divisor(m, ltg)
        if pos == 0
            push!(quo, m)
            @inbounds for i = 1:nbv
                pp = PP(m)
                pp[i] += 1
                if (find_in_list(pp, inspected) == 0)
                    push!(todo, pp)
                    push!(inspected, pp)
                end
            end
        end
    end
    return (sort(quo, lt = (e1, e2) -> pp_isless_drl(e1, e2)))
end

function isuniv(v::PP)
    pos = 0
    for i in eachindex(v)
        if (v[i] != 0)
            if (pos > 0)
                return (0)
            else
                pos = i
            end
        end
    end
    return (1)
end

function count_univ(l::PP)
    nb = 0
    @inbounds for i in eachindex(l)
        nb += l[i]
    end
    return (nb)
end

function mul_pp_by_var!(m::PP, v_index)
    mm = PP(m)
    mm[v_index] += 1
    return (mm)
end

function find_in_border(m::PP, t_xw::Vector{StackVect})
    pos = length(t_xw)
    res_flag = 0
    res_dd = 0
    res_pos = 0
    # On y va par degrÃ©s dÃ©croissants pour limiter
    # le nombre de tests
    tm = total_degree(m)
    while (pos > 0)
        mm = t_xw[pos].mon
        tmm = total_degree(mm)
        if (tmm == tm)
            if (m == mm)
                return (1, m, pos)
            end
        elseif (tmm == tm - 1)
            #test if mm is not in the quotient basis
            if (!((t_xw[pos].prev > 0) && (t_xw[pos].var == 0)))
                flag, dd = divides(m, mm)
                if (flag)
                    res_flag = 2
                    res_dd = dd
                    res_pos = pos
                end
            end
        else
            #the total degree is to high to find a predecessor
            if (res_flag > 0)
                return (res_flag, res_dd, res_pos)
            else
                error("\nError find in border : ", tm, " ", tmm)
                return (0, 0, 0)
            end
        end
        pos = pos - 1
    end
    if (res_flag > 0)
        return (res_flag, res_dd, res_pos)
    else
        error("\n\nError in finding a predecessor")
        return (0, 0, 0)
    end
end

# ************************************************************************
# Generic coefficients
# ************************************************************************

reduce_mod!(vres, arithm) = reduce_mod!(vres, arithm, vres)

function reduce_mod!(v::Vector{Accum}, arithm, vres::Vector{Coeff}) where {Accum,Coeff}
    @fastmath @inbounds @simd for i in eachindex(vres)
        vres[i] = ModularImage(Accum(v[i]), arithm) % Coeff
    end
    return vres
end

@noinline function add_mul!(vres::Vector{Accum}, a::Coeff, v::Vector{Coeff}) where {Accum,Coeff}
    @inbounds for i in eachindex(vres)
        vres[i] += Accum(a) * Accum(v[i])
    end
end

function vectorize_pol_gro!(p_exps, p_coeffs, kb, arithm::ModularArithmetic{Accum,Coeff}, res) where {Accum,Coeff}
    pos = length(kb)
    @inbounds for i = 2:length(p_coeffs)
        m = p_exps[i]
        while ((pos > 0) && (kb[pos] != m))
            pos = pos - 1
        end
        if (pos < 1)
            error("\nVectorisation ", pos, " ", m)
        else
            res[pos] = (ModularPrime(arithm) % Coeff) - p_coeffs[i]
        end
    end
    return res
end

function compute_fill_quo_gb!(t_xw, gro_exps, gro_coeffs, quo, arithm::ModularArithmetic{Accum,Coeff}) where {Accum,Coeff}
    t_v = Vector{Vector{Coeff}}(undef, length(t_xw))
    @inbounds for i in eachindex(t_xw)
        if ((t_xw[i].var > 0) && (t_xw[i].prev == 0))
            #var>0 and prev=0 => groebner basis element at position var in the list gro
            t_v[i] = zeros(Coeff, length(quo))
            vectorize_pol_gro!(gro_exps[t_xw[i].var], gro_coeffs[t_xw[i].var], quo, arithm, t_v[i])
        elseif ((t_xw[i].var == 0) && (t_xw[i].prev > 0))
            #var=0 (and prev>0 )=> quotient element at position prev in the list quo
            t_v[i] = [t_xw[i].prev]
        else
            t_v[i] = Vector{Coeff}()
        end
    end
    return t_v
end

function prepare_table_mxi(ltg, kb)
    nbv = length(ltg[1])
    tablex = [[Int32(0) for i = 1:length(kb)] for j = 1:nbv]
    general_stack = Vector{StackVect}()
    nb_stack = 0
    for j in eachindex(kb)
        m = PP(kb[j])
        for ii = 1:nbv
            i = nbv - ii + 1
            nm = mul_pp_by_var!(m, i)
            pos = findfirst(item -> item == nm, kb)
            if isnothing(pos)
                pos = findfirst(item -> item == nm, ltg)
                if isnothing(pos)
                    flag, v, prec = find_in_border(nm, general_stack)
                    if (flag == 1)
                        tablex[i][j] = Int32(prec)
                    elseif (flag == 2)
                        nb_stack = nb_stack + 1
                        tablex[i][j] = Int32(nb_stack)
                        push!(general_stack, StackVect(nb_stack, PP(nm), Int32(prec), Int32(v)))
                    else
                        error("\n*** Error search table ***\n")
                    end
                else
                    #nm is a leading monomial of an element of the Gb
                    #we insert it with the flags prev=0 and var=pos in the gb
                    prev = findlast(item -> item.mon == nm, general_stack)
                    if isnothing(prev)
                        nb_stack = nb_stack + 1
                        tablex[i][j] = Int32(nb_stack)
                        push!(general_stack, StackVect(nb_stack, PP(ltg[pos]), Int32(0), Int32(pos)))
                    else
                        tablex[i][j] = Int32(prev)
                    end
                end
            else
                #nm is a element of the quotient 
                #we insert it with the flags prev=pos and var=0
                prev = findlast(item -> item.mon == nm, general_stack)
                if isnothing(prev)
                    nb_stack = nb_stack + 1
                    tablex[i][j] = Int32(nb_stack)
                    push!(general_stack, StackVect(nb_stack, PP(kb[pos]), Int32(pos), Int32(0)))
                else
                    tablex[i][j] = Int32(prev)
                end
            end
        end
    end
    return (tablex, general_stack)
end

function test_zdim(de, co)
    count_univ(map(w -> isuniv(w), map(u -> u[1], de))) == length(de[1][1])
end

TimerOutputs.@timeit to "Mul Var Quo" function _mul_var_quo!(
    v,
    ii,
    t_v,
    i_xw,
    arithm::ModularArithmetic{Accum,Coeff},
    vv,
    vres,
) where {Accum,Coeff}
    dim = length(v)
    @assert dim == length(vres)
    pack = pack_value(arithm)
    resize!(vv, dim)
    @inbounds for i = 1:dim
        vv[i] = zero(Accum)
    end
    continuer = true
    @inbounds for j = 1:dim
        iszero(v[j]) && continue
        if (length(t_v[i_xw[ii][j]]) > 1)
            add_mul!(vv, v[j], t_v[i_xw[ii][j]])
            if j % pack == 0
                reduce_mod!(vv, arithm)
            end
        elseif (length(t_v[i_xw[ii][j]]) == 1)
            kk = t_v[i_xw[ii][j]][1]
            vv[kk] = ModularImage(vv[kk] + Accum(v[j]), arithm)
        else
            continuer = false
            break
        end
    end
    return continuer, reduce_mod!(vv, arithm, vres)
end

function mul_var_quo(v, ii, t_v, i_xw, arithm::ModularArithmetic{Accum,Coeff}, vv) where {Accum,Coeff}
    vres = Vector{Coeff}(undef, length(v))
    f, r = _mul_var_quo!(v, ii, t_v, i_xw, arithm, vv, vres)
    return r
end

function learn_compute_table!(t_v, t_xw, i_xw, quo, arithm)
    nb = 1
    t_learn = Int32[]
    buf = Vector{AccModularCoeff}()
    while (nb > 0)
        nb = 0
        for i in eachindex(t_xw)
            #test if already computed
            if (length(t_v[i]) == 0)
                continuer = false
                #test if the ancestor is computed
                if (length(t_v[t_xw[i].prev]) > 0)
                    vv = Vector{ModularCoeff}(undef, length(t_v[t_xw[i].prev]))
                    continuer, vv = _mul_var_quo!(t_v[t_xw[i].prev], t_xw[i].var, t_v, i_xw, arithm, buf, vv)
                    if (continuer)
                        t_v[i] = vv
                        push!(t_learn, i)
                        nb = nb + 1
                    end
                end
            end
        end
    end
    return (t_learn)
end

TimerOutputs.@timeit to "Apply Table" function apply_compute_table!(
    t_v,
    t_learn,
    t_xw,
    i_xw,
    quo,
    arithm::ModularArithmetic{Accum,Coeff},
) where {Accum,Coeff}
    buf = Vector{Accum}()
    @inbounds for i in t_learn
        if (length(t_v[i]) == 0)
            t_v[i] = mul_var_quo(t_v[t_xw[i].prev], t_xw[i].var, t_v, i_xw, arithm, buf)
        end
    end
    return nothing
end

function normalize_row!(v::Vector{Coeff}, tmp::Coeff, arithm::ModularArithmetic{Accum,Coeff}) where {Accum,Coeff}
    @fastmath @inbounds @simd for i in eachindex(v)
        v[i] = ModularImage(Accum(v[i]) * Accum(tmp), arithm) % Coeff
    end
end

TimerOutputs.@timeit to "Gauss Reduct" function gauss_reduct(
    v::Vector{Coeff},
    gred::Vector{Vector{Coeff}},
    dg::Int32,
    dv::Int32,
    paquet::Int32,
    arithm::ModularArithmetic{Accum,Coeff},
    b::Vector{Accum},
) where {Coeff,Accum}
    j = Int32(0)
    last_nn = Int32(-1)
    resize!(b, dv)
    @inbounds for ell = 1:dv
        b[ell] = Accum(v[ell])
    end
    @inbounds for i = 1:(dg-1)
        b[i+1] = ModularImage(b[i+1], arithm)
        iszero(b[i+1]) && continue
        if (length(gred[i]) == 0)
            last_nn = i % Int32
            break
        end
        pivot::Coeff = ((ModularPrime(arithm) % Coeff) - (b[i+1] % Coeff)) % Coeff
        add_mul!(b, pivot, gred[i])
        b[i+1] = pivot
        if (j < paquet)
            j += Int32(1)
        else
            reduce_mod!(b, arithm)
            j = Int32(0)
        end
    end
    if (j > 0)
        @inbounds for k = 1:dv
            v[k] = ModularImage(b[k], arithm) % Coeff
        end
    else
        @inbounds for k = 1:dv
            v[k] = b[k] % Coeff
        end
    end
    if (last_nn == -1)
        @inbounds for ii = dg:(dv-1)
            if (!iszero(v[ii+1]))
                last_nn = ii % Int32
                break
            end
        end
    end
    if last_nn == -1
        return dv
    else
        return last_nn
    end
end

TimerOutputs.@timeit to "First Variable" function first_variable(
    t_v::Vector{Vector{Coeff}},
    i_xw::Vector{Vector{Int32}},
    ii::Int32,
    arithm::ModularArithmetic{Accum,Coeff},
) where {Accum,Coeff}
    pack = pack_value(arithm)
    d = Int32(length(i_xw[ii]))
    free_set = [append!([one(Coeff)], [zero(Coeff) for i = 2:d])]
    if (length(t_v[i_xw[ii][1]]) > 1)
        v = t_v[i_xw[ii][1]]
    else
        v = zeros(Coeff, d)
        v[t_v[i_xw[ii][1]][1]] = one(Coeff)
    end
    push!(free_set, v)
    # index[deg]=pos (pos=position in the gred with the convention that 0 is not stored)
    index = [Int32(i) for i = 1:d]
    # normalization values
    hom = [one(Coeff) for i = 1:d]
    gred = [Vector{Coeff}() for i = 1:d]

    i = Int32(2)
    continuer = 1
    dg = Int32(1)
    new_i = Int32(1)
    deg = Int32(0)
    while ((i < d) && (v[i] == 0))
        i += 1
    end
    gred[i-1] = Vector{Coeff}(v)
    if (gred[i-1][1] != 0)
        gred[i-1][1] = (ModularPrime(arithm) % Coeff) - gred[i-1][1]
    end
    if (gred[i-1][i] != 1)
        hom[i-1] = ModularInverse(Accum(gred[i-1][i]), arithm)
        normalize_row!(gred[i-1], hom[i-1], arithm)
    end
    #T^0 is stored at gred[0] virtually
    index[1] = i - 1
    #gred[dg+i] not used i=0..d-dg
    dg = Int32(i)
    #we are now going to compute T^2
    deg = Int32(2)
    w = Vector{Coeff}(undef, d)
    buf1 = Vector{Accum}(undef, d)
    buf2 = Vector{Accum}(undef, d)
    @inbounds while (continuer == 1)
        v = mul_var_quo(v, ii, t_v, i_xw, arithm, buf1)
        resize!(w, length(v))
        for ell = 1:length(v)
            w[ell] = v[ell]
        end
        #reduce with gred[0]
        if (w[1] != 0)
            w[1] = (ModularPrime(arithm) % Coeff) - w[1]
        end
        #reduce with gred[i] i=1..dg-1
        #new_i = possible position to store the reduced vector in gred
        new_i = gauss_reduct(w, gred, dg, d, pack, arithm, buf2)
        if (new_i < d)
            push!(free_set, v)
            gred[new_i] = copy(w)
            if (!(new_i < dg))
                dg = Int32(new_i + 1)
            end
            index[deg] = Int32(new_i)
            hom[new_i] = ModularInverse(Accum(gred[new_i][new_i+1]), arithm) % Coeff
            normalize_row!(gred[new_i], hom[new_i], arithm)
            deg += Int32(1)
        else
            continuer = 0
        end
    end
    #set v[i] cofficient of T^(i-1) in the min poly
    v = Vector{Coeff}(undef, deg)
    v[1] = w[1]
    @inbounds for i = 2:deg
        v[i] = Coeff(ModularImage(Accum(w[index[i-1]+1]) * Accum(hom[index[i-1]]), arithm))
    end
    return (v, gred, index, dg, hom, free_set)
end

TimerOutputs.@timeit to "Biv Lex" function biv_lex(t_v, i_xw, gred, index, dg, hom, free_set, ii, arithm)
    pack = pack_value(arithm)
    d = Int32(length(i_xw[length(i_xw)]))
    new_free_set = copy(free_set)
    deg = length(free_set)
    new_generators = Vector{Vector{Int32}}()
    new_monomial_free_set = [(Int32(i - 1), Int32(0)) for i = 1:deg]
    new_monomial_basis = copy(new_monomial_free_set)
    new_leading_monomials = Vector{Tuple{Int32,Int32}}()
    buf1 = Vector{AccModularCoeff}()
    buf2 = Vector{AccModularCoeff}()
    while (length(new_free_set) > 0)
        tmp_set = Vector{Vector{ModularCoeff}}()
        tmp_mon_set = Vector{Tuple{Int32,Int32}}()
        @inbounds for j in eachindex(new_free_set)
            curr_mon = new_monomial_free_set[j]
            curr_mon = (curr_mon[1], curr_mon[2] + Int32(1))
            v = mul_var_quo(new_free_set[j], ii, t_v, i_xw, arithm, buf1)
            w = Vector{ModularCoeff}(v)
            if (w[1] != 0)
                w[1] = ModularPrime(arithm) - w[1]
            end
            new_i = gauss_reduct(w, gred, dg, d, pack, arithm, buf2)
            if (new_i < d)
                push!(tmp_set, v)
                push!(tmp_mon_set, curr_mon)
                gred[new_i] = Vector{ModularCoeff}(w)
                if (!(new_i < dg))
                    dg = Int32(new_i + 1)
                end
                index[deg] = Int32(new_i)
                hom[new_i] = ModularCoeff(invmod(gred[new_i][new_i+1], ModularPrime(arithm)))
                normalize_row!(gred[new_i], hom[new_i], arithm)
                deg += Int32(1)
            else
                v = Vector{ModularCoeff}(undef, deg)
                v[1] = w[1]
                @inbounds for i = 2:deg
                    v[i] = ModularCoeff(
                        ModularImage((AccModularCoeff(w[index[i-1]+1])) * (AccModularCoeff(hom[index[i-1]])), arithm),
                    )
                end
                push!(new_generators, copy(v))
                push!(new_leading_monomials, curr_mon)
                break
            end
        end
        new_free_set = copy(tmp_set)
        new_monomial_free_set = copy(tmp_mon_set)
        append!(new_monomial_basis, copy(tmp_mon_set))
    end
    #    rur_print("{",length(new_monomial_basis),"}")
    return (new_monomial_basis, new_leading_monomials, new_generators)
end

function convert_biv_lex_2_biv_pol(n_g, m_b, lt_g)
    l_base = Vector{Vector{Vector{Int32}}}()
    for kk in eachindex(n_g)
        pp = n_g[kk]
        p_mat = Vector{Vector{ModularCoeff}}()
        ldeg2 = Int32(0)
        lco = [ModularCoeff(0) for i = 1:length(m_b)]
        for i in eachindex(pp)
            deg2 = m_b[i][2]
            if deg2 > ldeg2
                push!(p_mat, lco)
                for j = (ldeg2+1):(deg2-1)
                    push!(p_math, Vector{ModularCoeff}())
                end
                lco = [ModularCoeff(0) for j = 1:length(m_b)]
                ldeg2 = deg2
            end
            lco[m_b[i][1]+1] = pp[i]
        end
        deg2 = lt_g[kk][2]
        if (deg2 == ldeg2)
            lco[lt_g[kk][1]+1] = ModularCoeff(1)
        else
            push!(p_mat, lco)
            for i = (ldeg2+1):(deg2-1)
                push!(p_math, Vector{ModularCoeff}())
            end
            lco = [ModularCoeff(0) for i = 1:length(m_b)]
            lco[lt_g[kk][1]+1] = ModularCoeff(1)
        end
        push!(p_mat, lco)
        push!(l_base, p_mat)
    end
    return (l_base)
end

function check_separation_biv(bli, f, C)
    k = length(bli) - 1
    akk = C(bli[k+1])
    invakk = Nemo.invmod(k * akk, f)
    b = C(bli[k]) * invakk
    for i = (k-1):-1:1
        tmp = Nemo.mod((k - i + 1) * C(bli[i]) - (i) * b * C(bli[i+1]), f)
        if (tmp != 0)
            return (false)
        end
    end
    return (true)
end

#Nemo functions : derivative, gcd, degree mulmod invmod 
function zdim_parameterization(t_v, i_xw, dd, check, arithm)
    res = Vector{Vector{ModularCoeff}}()
    ii = Int32(length(i_xw))
    v, gred, index, dg, hom, free_set = first_variable(t_v, i_xw, ii, arithm)
    # here we can process with a loop over several primes to get f and ifp
    # as vectors of Composite numbers : just convert the vector v to a list of vectors

    if ModularCoeff isa CNCoeff
        l_pr = ModularPrime(arithm).data
    else
        l_pr = [ModularPrime(arithm)]
    end
    flag = true
    l_res = []
    for pr in l_pr
        res = []
        C, _Z = Nemo.polynomial_ring(Nemo.Native.GF(Int64(pr), cached = false), cached = false)
        f = C(v) + _Z^(length(v))
        ifp = Nemo.derivative(f)
        f = f / Nemo.gcd(f, ifp)
        ifp = Nemo.derivative(f)
        push!(res, map(u -> coeff_mod_p(u), collect(Nemo.coefficients(f))))
        push!(l_res, res)
        if (dd < 0)
            flag = true
        else
            flag = (flag && (Int32(Nemo.degree(f)) == dd))
        end
    end
    # end of the first loop over several primes
    if (flag)
        @inbounds for j = (ii-1):-1:1
            m_b, lt_b, n_g = biv_lex(t_v, i_xw, gred, index, dg, hom, free_set, Int32(j), arithm)
            bl = [convert_biv_lex_2_biv_pol(n_g, m_b, lt_b)]
            nbp = 0
            # as above we can process with a loop over several primes to get
            for pr in l_pr
                nbp += 1
                C, _Z = Nemo.polynomial_ring(Nemo.Native.GF(Int64(pr), cached = false), cached = false)
                s1 = C([Int32(0)])
                s0 = C([Int32(0)])
                pro = C([Int32(1)])
                #                ft=C(f)
                f = C(l_res[nbp][1])
                ft = C(l_res[nbp][1])
                ifp = Nemo.derivative(ft)
                @inbounds for i = 1:length(bl[nbp])
                    d1 = length(bl[nbp][i]) - 1
                    lc1 = C(bl[nbp][i][d1+1])
                    co0 = C(bl[nbp][i][d1])
                    f2 = Nemo.gcd(ft, lc1)
                    f1 = ft / f2
                    # in some non radical case
                    if (Nemo.degree(f1) > 0)
                        if (check > 0)
                            flag = check_separation_biv(bl[nbp][i], f1, C)
                            if (!flag)
                                return (false, [], j)
                            end
                        end
                        s1 += Nemo.mulmod(d1 * lc1, pro, f)
                        s0 += Nemo.mulmod(co0, pro, f)
                        pro = pro * f1
                    end
                    ft = C(f2)
                end
                is1 = Nemo.invmod(s1, f)
                s0 = Nemo.mulmod(s0, is1, f)
                s0 = -Nemo.mulmod(s0, ifp, f)
                #                push!(res,map(u->coeff_mod_p(u),collect(Nemo.coefficients(s0))))
                push!(l_res[nbp], map(u -> coeff_mod_p(u), collect(Nemo.coefficients(s0))))
            end
        end
    else
        rur_print("Bad prime number for parameterization (", dd, ",", length(v), ")")
        #error("Bad prime number for parameterization ")
    end
    if (check>0) return(flag,l_res,ii)
    else return (flag, l_res) end
end


function _zdim_modular_RUR_LV(de, cco, arithm)
    rur_print("Gb - ")
    ex, co = _gb_4_rur(de, cco, arithm)
    ltg = map(u -> u[1], ex)
    rur_print("Quo - ")
    q = compute_quotient_basis(ltg)
    rur_print("Prep table - ")
    i_xw, t_xw = prepare_table_mxi(ltg, q)
    rur_print("Fill Gb - ")
    t_v = compute_fill_quo_gb!(t_xw, ex, co, q, arithm)
    rur_print("Learn comp Table - ")
    t_learn = learn_compute_table!(t_v, t_xw, i_xw, q, arithm)
    rur_print("Learn Parameterization \n")
    flag, zp_param,uu = zdim_parameterization(t_v, i_xw, Int32(-1), Int32(1), arithm)
    return (flag, zp_param, ltg, q, i_xw, t_xw, t_learn,uu)
end

function _zdim_modular_RUR_LV_apply!(de, co, arithm, dd, ltg, q, i_xw, t_xw, t_learn, graph)
    t_v = compute_fill_quo_gb!(t_xw, graph.gb_support, co, q, arithm)
    apply_compute_table!(t_v, t_learn, t_xw, i_xw, q, arithm)
    success, zp_param = zdim_parameterization(t_v, i_xw, Int32(dd), Int32(0), arithm)
    return (success, zp_param)
end

function swap_vars(m::PP, ii, jj)
    res = PP(m)
    u = res[ii]
    res[ii] = res[jj]
    res[jj] = u
    return (res)
end

swap_vars(p::Vector{PP}, ii, jj) = map(u -> swap_vars(u, ii, jj), p)
swap_vars(de::Vector{Vector{PP}}, ii, jj) = map(u -> swap_vars(u, ii, jj), de)

extend_vars(m::PP, ii) = (vcat(m, [Deg(0) for i = 1:ii]))
extend_vars(p::Vector{PP}, ii) = map(u -> extend_vars(u, ii), p)
extend_vars(de::Vector{Vector{PP}}, ii) = map(u -> extend_vars(u, ii), de)

function modular_coeffs_vect(v, pr)
    map(u -> (u < 0 ? ModularCoeff(pr - ((-u) % pr)) : ModularCoeff(u % pr)), v)
end

#should be independent from the used arithmetic

function extend_system_with_sepform(de, co, lt, minus_sep_lin)
    nbv = length(de[1][1])
    dde = extend_vars(de, 1)
    pp = [Deg(0) for j = 1:(nbv+1)]
    pp[nbv+1] = 1
    co_sep_pol = [lt]
    exp_sep_pol = [pp]
    for i = 1:(nbv)
        if !iszero(minus_sep_lin[i])
            pp = [Deg(0) for j = 1:(nbv+1)]
            pp[i] = 1
            push!(exp_sep_pol, pp)
            push!(co_sep_pol, minus_sep_lin[i])
        end
    end
    push!(dde, exp_sep_pol)
    cco = copy(co)
    push!(cco, co_sep_pol)
    return (dde, cco)
end

function _zdim_modular_RUR(de, co, arithm, learn = false)
    nbv = length(de[1][1])
    ii = nbv
    flag, zp_param, ltg, q, i_xw, t_xw, t_learn = _zdim_modular_RUR_LV(de, co, arithm)
    while ((!flag) && (ii > 1))
        ii = ii - 1
        dde = swap_vars(de, ii, nbv)
        flag, zp_param, ltg, q, i_xw, t_xw, t_learn = _zdim_modular_RUR_LV(dde, co, arithm)
    end
    sep_lin = [0 for i = 1:nbv]
    if (flag)
        sep_lin[ii] = 1
    else
        sep_lin=[ 0 for i=1:nbv]
        sep_lin[nbv]=1
        sep_lin[nbv-1]=-1
        vv=nbv-2
        # naive strategy  sep_lin = [rand(-100:100) for i = 1:nbv]
        while (!flag)
            dde, cco = extend_system_with_sepform(de,co,ModularCoeff(1),modular_coeffs_vect(map(u -> -u, sep_lin), ModularPrime(arithm)))
            flag, zp_param, ltg, q, i_xw, t_xw, t_learn, uu = _zdim_modular_RUR_LV(dde, cco, arithm)           
            if (!flag) 
                #  naive strategy          sep_lin = [rand(-100:100) for i = 1:nbv]
                rur_print(" (",uu,",",vv,")")
                if (sep_lin[uu]<0) sep_lin[uu]=-sep_lin[uu]
                else sep_lin[uu]=-sep_lin[uu]-1 end
                vv=uu
                if (vv<1)
                    rur_print("\n Error in the choice of the separating form \n") 
                    return(-1,-1,sys,AbstractAlgebra.symbols(C),linform)
                end
            else
                dd=length(zp_param[1])-1
            end  
        end
    end
    if (learn)
        return (flag, sep_lin, zp_param, ltg, q, i_xw, t_xw, t_learn)
    else
        return (flag, sep_lin, zp_param)
    end
end

function zz_mat_same_dims(zpm)::Vector{Vector{BigInt}}
    zzm = Vector{Vector{BigInt}}(undef, length(zpm))
    @inbounds for i in eachindex(zpm)
        zzm[i] = [BigInt(0) for j = 1:length(zpm[i])]
    end
    return (zzm)
end

function qq_mat_same_dims(zpm)::Vector{Vector{Rational{BigInt}}}
    zzm = Vector{Vector{Rational{BigInt}}}(undef, length(zpm))
    @inbounds for i in eachindex(zpm)
        zzm[i] = [Rational{BigInt}(0) for j = 1:length(zpm[i])]
    end
    return (zzm)
end

function learn_compute_table_cyclic!(t_v, t_xw, i_xw, quo, arithm)
    nb = 1
    t_learn = Int32[]
    buf = Vector{AccModularCoeff}()
    #xnwi, i=1..D
    to_compute = copy(i_xw[end])
    #xiw1 i=1..n-1
    for i = 1:(length(i_xw)-1)
        push!(to_compute, i_xw[i][1])
    end
    while (nb > 0)
        nb = 0
        for i in eachindex(to_compute)
            #test if already computed
            if (length(t_v[i]) == 0)
                continuer = false
                #test if the ancestor is computed
                if (length(t_v[t_xw[i].prev]) > 0)
                    vv = Vector{ModularCoeff}(undef, length(t_v[t_xw[i].prev]))
                    continuer, vv = _mul_var_quo!(t_v[t_xw[i].prev], t_xw[i].var, t_v, i_xw, arithm, buf, vv)
                    if (continuer)
                        t_v[i] = vv
                        push!(t_learn, i)
                        nb = nb + 1
                    end
                end
            end
        end
    end
    return (t_learn)
end

function _list_zdim_modular_RUR_LV_apply_serial!(de, cco, lp, dd, ltg, q, i_xw, t_xw, t_learn, graph, composite)
    success = Vector{Bool}(undef, length(lp))
    res = Vector{Vector{Vector{ModularCoeff}}}(undef, length(lp))

    @assert length(lp) % composite == 0
    for i in 1:composite:length(lp)
        arithm_Nx = map(pr -> ModularArithZp(AccModularCoeff, ModularCoeff, ModularCoeff(pr)), lp[i:i+composite-1])
        co_Nx = map(pr -> map(v -> modular_coeffs_vect(v, pr), cco), lp[i:i+composite-1])
        flag1, gb_Nx = _gb_4_rur_apply!(graph, de, co_Nx, arithm_Nx)
        j = 1
        while j <= composite
            !flag1 && (success[i+j-1] = false; continue)
            pr = lp[i+j-1]
            arithm = ModularArithZp(AccModularCoeff, ModularCoeff, ModularCoeff(pr))
            flag2, l_zp_param = _zdim_modular_RUR_LV_apply!(de, gb_Nx[j], arithm, dd, ltg, q, i_xw, t_xw, t_learn, graph)
            !flag2 && (success[i+j-1] = false; continue)
            success[i+j-1] = true
            res[i+j-1] = l_zp_param[1]
            j += 1
        end
    end

    return (success, res)
end

function _list_zdim_modular_RUR_LV_apply_parallel!(de, cco, lp, dd, ltg, q, i_xw, t_xw, t_learn, graph, composite, threads)
    success = Vector{Bool}(undef, length(lp))
    res = Vector{Vector{Vector{ModularCoeff}}}(undef, length(lp))

    @assert length(lp) % (threads * composite) == 0
    tstep = div(length(lp), threads)
    @threads :static for T in 1:threads
        for i in (1+(T-1)*tstep):composite:((T)*tstep)
            t_graph = graph[threadid()]
            arithm_Nx = map(pr -> ModularArithZp(AccModularCoeff, ModularCoeff, ModularCoeff(pr)), lp[i:i+composite-1])
            co_Nx = map(pr -> map(v -> modular_coeffs_vect(v, pr), cco), lp[i:i+composite-1])
            flag1, gb_Nx = _gb_4_rur_apply!(t_graph, de, co_Nx, arithm_Nx)
            j = 1
            while j <= composite
                !flag1 && (success[i+j-1] = false; continue)
                pr = lp[i+j-1]
                arithm = ModularArithZp(AccModularCoeff, ModularCoeff, ModularCoeff(pr))
                flag2, l_zp_param = _zdim_modular_RUR_LV_apply!(de, gb_Nx[j], arithm, dd, ltg, q, i_xw, t_xw, t_learn, t_graph)
                !flag2 && (success[i+j-1] = false; continue)
                success[i+j-1] = true
                res[i+j-1] = l_zp_param[1]
                j += 1
            end
        end
    end

    return (success, res)
end

function rur_check(de, cco, pr, qq_m)
    co = map(_c -> map(__c -> mod(numerator(__c) * invmod(denominator(__c), pr), pr) % UInt32, _c), cco)
    arithm = ModularArithZp(UInt64, UInt32, UInt32(pr))
    flag, rur1, _, _, _, _, _, _ = _zdim_modular_RUR_LV(de, co, arithm)
    @assert flag
    rur2 = map(_c -> map(__c -> mod(numerator(__c) * invmod(denominator(__c), pr), pr) % UInt32, _c), qq_m)
    rur1[1] == rur2
end

TimerOutputs.@timeit to "MM loop" function _zdim_multi_modular_RUR!(de, cco, bit_pr=pr_max_bitsize,parallelism=:serial,composite=4,threads=1)
    nbv_ori = length(de[1][1])
    pr = ModularCoeff(PrevPrime(2^bit_pr - 1))
    rur_print("primes of bitsize ", bit_pr, "\n")
    arithm = ModularArithZp(UInt64, UInt32, pr)
    co = map(v -> modular_coeffs_vect(v, pr), cco)
    flag, sep_lin, l_zp_param, ltg, q, i_xw, t_xw, t_learn = _zdim_modular_RUR(de, co, arithm, true)
    dd = length(l_zp_param[1][1]) - 1

    cyclic = (dd == length(q))

    if (length(ltg[1]) == nbv_ori)
        if (iszero(sep_lin[nbv_ori]))
            i = findfirst(item -> !iszero(item), de)
            if (isnothing(i))
                error("Misseformed separating element")
            else
                de = swap_vars(de, i, nbv_ori)
                rur_print("Use variable", i, " as separating element\n")
            end
        else
            rur_print("Use last variable as separating element\n")
        end
    else
        rur_print("Use non trivial separating element\n")
        de, cco = extend_system_with_sepform(de, cco, BigInt(1), map(u -> -u, sep_lin))
    end

    rur_print("Run Groebner learn\n")
    co = map(v -> modular_coeffs_vect(v, pr), cco)
    graph, gex, gco = _gb_4_rur_learn(de, co, arithm)

    if (cyclic)
        rur_print("Test cyclic optimization")
        t_v2 = compute_fill_quo_gb!(t_xw, gex, gco, q, arithm)
        t_learn2 = learn_compute_table_cyclic!(t_v2, t_xw, i_xw, q, arithm)
        success, l_zp_param2 = zdim_parameterization(t_v2, i_xw, Int32(dd), Int32(1), arithm)
        if (success)
            rur_print("\nApply cyclic optimization \n")
            t_v = t_v2
            t_learn = t_learn2
            l_zp_param = l_zp_param2
        else
            rur_print("\nSwitch off cyclic optimization \n")
        end
    end

    rur_print("Multi-modular computation ($threads threads): ")
    t_pr = Vector{ModularCoeff}()
    t_param = Vector{Vector{Vector{ModularCoeff}}}()
    continuer = true

    zz_m = zz_mat_same_dims(l_zp_param[1])
    qq_m = qq_mat_same_dims(l_zp_param[1])

    witness = reduce(vcat, [[(poly,rand(1:length(l_zp_param[1][poly]))) for _ in 1:ceil(log2(length(l_zp_param[1][poly])))] for poly in 1:length(l_zp_param[1])])
    crt_mask = [falses(length(_x)) for _x in zz_m]
    ratrec_mask = [falses(length(_x)) for _x in zz_m]

    append!(t_pr, pr)
    append!(t_param, l_zp_param)
   
    bloc_p = 2
    ALIGN_BLOC_TO = 1 * composite
   
    if parallelism!=:serial
        TimerOutputs.disable_timer!(to)
        graph = map(_ -> deepcopy(graph), 1:nthreads())
        ALIGN_BLOC_TO *= threads
    end

    align_to(x,n) = (x + (n - 1)) & (~(n - 1))
    bloc_p = align_to(bloc_p, ALIGN_BLOC_TO)

    while (continuer)
        l_pr = PrevPrimes(pr - 1, bloc_p)
        pr = l_pr[end]

        if parallelism==:serial
             success, l_zp_param = _list_zdim_modular_RUR_LV_apply_serial!(de, cco, l_pr, dd, ltg, q, i_xw, t_xw, t_learn, graph, composite)
        elseif parallelism==:multithreading
            success, l_zp_param = _list_zdim_modular_RUR_LV_apply_parallel!(de, cco, l_pr, dd, ltg, q, i_xw, t_xw, t_learn, graph, composite, threads)
        else
            error("Unknown parallelism strategy ",parallelism)
        end

        for i in eachindex(success)
            !success[i] && (rur_print("bad prime $(l_pr[i])-"); continue)
            push!(t_pr, l_pr[i])
            push!(t_param, l_zp_param[i])
        end
        
        rur_print(length(t_pr), "-")
        zz_p = BigInt(1)
        TimerOutputs.@timeit to "crt partial" Groebner.crt_vec_partial!(zz_m, zz_p, t_param, t_pr, witness, crt_mask)
        TimerOutputs.@timeit to "rat rec partial" continuer = !Groebner.ratrec_vec_partial!(qq_m, zz_m, zz_p, witness, ratrec_mask)
        
        if !continuer
            TimerOutputs.@timeit to "crt" Groebner.crt_vec_full!(zz_m, zz_p, t_param, t_pr, crt_mask)
            TimerOutputs.@timeit to "rat rec" continuer = !Groebner.ratrec_vec_full!(qq_m, zz_m, zz_p, ratrec_mask)
            if !continuer
                continuer = !rur_check(de, cco, PrevPrime(pr - 1), qq_m)
            end
        end

        bloc_p = max(floor(Int, length(t_pr) / 10), 2)
        bloc_p = align_to(bloc_p, ALIGN_BLOC_TO)
    end
    rur_print("\n")
    return qq_m,sep_lin
end

# ************************************************************************
# User Interfaces
# ************************************************************************


# ************************************************************************
# Abstract Algebra Interface
# ************************************************************************

function zdim_parameterization(
        sys_ori;
        nn::Int32=Int32(28),
        verbose::Bool=true,
        parallelism=:serial,
        threads=(parallelism == :multithreading ? nthreads() : 1),
	    get_separating_element::Bool=false,
        composite=4)
    @assert 1 <= nn <= 32
    @assert parallelism in (:serial, :multithreading) && 1 <= threads <= nthreads()
    parallelism == :serial && threads > 1 && rur_print("WARN: threads=$threads was ignored\n")
    @assert 1 <= composite && ispow2(composite)
    @assert AbstractAlgebra.base_ring(AbstractAlgebra.parent(sys_ori[1])) == AbstractAlgebra.QQ
    _verbose[]=verbose
    TimerOutputs.enable_timer!(to)
    sys=sys_ori.*(map(v->lcm(map(w->denominator(w),collect(AbstractAlgebra.coefficients(v)))),sys_ori))
    de = map(p -> collect(AbstractAlgebra.exponent_vectors(p)), sys)
    de = map(u -> map(v -> map(w -> map(h -> Deg(h), w), v), u), de)
    co = map(p -> collect(AbstractAlgebra.coefficients(p)), sys)
    res,sep = _zdim_multi_modular_RUR!(de, co, nn , parallelism, composite, threads)
    if (get_separating_element) 
       return(res,sep)
    else
       return(res)
    end
end

using PrecompileTools
include("precompile.jl")

end # module

