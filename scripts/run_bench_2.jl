# julia -t 4 -p 4 run_bench_2.jl 
using Base.Threads, Distributed, SharedArrays

function f_threads()
    k = Array{Float64, 3}(undef, 100000, 100, 100)
    @threads for i = 1:10000
        k[i, :, :] = rand(Float64, (100, 100))
        @inbounds for j = 1:100
            k[i, :, :] .= sin.(rand(Float64, (100, 100))).^2
        end
    end
    return k
end

function f_procs()
    k = SharedArray{Float64}(100000, 100, 100)
    @sync @distributed for i = 1:10000
        k[i, :, :] = rand(Float64, (100, 100))
        @inbounds for j = 1:100
            k[i, :, :] .= sin.(rand(Float64, (100, 100))).^2
        end
    end
    return k
end

@info "Using $(nthreads()) threads"
@info "Using $(nworkers()) processes"

@time f_threads();
@time f_procs();

@time f_threads();
@time f_procs();

#=
N = 1:
 43.657884 seconds (2.76 M allocations: 82.792 GiB, 0.63% gc time, 0.49% compilation time)
 47.551134 seconds (889.50 k allocations: 60.183

N = 2:
 36.572905 seconds (2.68 M allocations: 82.787 GiB, 1.03% gc time, 1.00% compilation time)
 27.566359 seconds (890.43 k allocations: 60.244 MiB, 0.03% gc time, 1.65% compilation time)

N = 4:
 28.541435 seconds (2.68 M allocations: 82.787 GiB, 1.56% gc time, 2.73% compilation time)
 14.646976 seconds (891.22 k allocations: 60.300 MiB, 0.06% gc time, 3.12% compilation time)

N = 8:
 21.123327 seconds (2.03 M allocations: 82.747 GiB, 1.60% gc time, 0.14% compilation time)
 11.319487 seconds (1.84 k allocations: 99.297 KiB)
=#

#=
On another PC:

N = 1:
 244.747607 seconds (2.71 M allocations: 82.790 GiB, 2.75% gc time, 0.24% compilation time)
 281.859264 seconds (1.06 M allocations: 70.604 MiB, 0.06% gc time, 0.66% compilation time)

N = 2:
 120.446679 seconds (2.71 M allocations: 82.790 GiB, 4.88% gc time, 1.04% compilation time)
 135.672287 seconds (1.06 M allocations: 70.658 MiB, 0.06% gc time, 1.31% compilation time)

N = 4:
 54.183032 seconds (2.05 M allocations: 82.749 GiB, 7.41% gc time, 0.67% compilation time)
 59.482889 seconds (931 allocations: 40.305 KiB)

N = 8:
 24.010636 seconds (2.05 M allocations: 82.749 GiB, 15.24% gc time, 2.17% compilation time)
 22.953416 seconds (1.86 k allocations: 82.539 KiB)
=#
