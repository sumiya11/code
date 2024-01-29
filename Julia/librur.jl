module LibRUR
using Libdl

# *****************************************************
# Charge des fonctions de la bibliothèque RS 
# et initialise RS
# *****************************************************

function local_install()
   if Sys.iswindows()
      return("Windows_x86_64")
   else 
      os=readchomp(`uname -s`)
      mach=readchomp(`uname -m`)
      return(os*"_"*mach)
   end
end

function __init__()
#    rur_install_path=dirname(pwd())
#    lib_so_path=joinpath(rur_install_path,joinpath("bin",local_install()))
#    lib_so_rur=joinpath(lib_so_path,"librur")
    #    librur=Libdl.dlopen(lib_so_rur)
    librur=Libdl.dlopen("/Users/rouillie/Travail/Git/rur/C/librur.dylib")
    global _rur_test_lib = Libdl.dlsym(librur, :rur_test_lib)
    global _new_table_zp =  Libdl.dlsym(librur, :new_table_zp)
    global _table_zp_init_set_row =  Libdl.dlsym(librur, :table_zp_init_set_row)
    global _table_zp_get_row =  Libdl.dlsym(librur, :table_zp_get_row)
    global _min_poly_zp =  Libdl.dlsym(librur, :min_poly_zp)
end

function rur_test_lib(id::Vector{Int32},v::Int32)
    @ccall $_rur_test_lib(Base.unsafe_convert(Ptr{Int32}, id)::Ptr{Cvoid},v::Int32)::Cvoid
end

function new_table_zp(d::Int32)::Ptr{Cvoid}
    @ccall $_new_table_zp(d::Int32)::Ptr{Cvoid}
end

function table_zp_init_set_row(t::Ptr{Cvoid},index::Int32,v::Vector{Int32},d::Int32)::Cvoid
    @ccall $_table_zp_init_set_row(t::Ptr{Cvoid},index::Int32,Base.unsafe_convert(Ptr{Int32}, v)::Ptr{Cvoid},d::Int32)::Cvoid
end

function table_zp_get_row(v::Vector{Int32},t::Ptr{Cvoid},index::Int32,d::Int32)::Cvoid
    @ccall $_table_zp_get_row(Base.unsafe_convert(Ptr{Int32}, v)::Ptr{Cvoid},t::Ptr{Cvoid},index::Int32,d::Int32)::Cvoid
end
function min_poly_zp(v::Vector{Int32},t::Ptr{Cvoid},d::Int32,pr::Int32)::Cvoid
    @ccall $_min_poly_zp(Base.unsafe_convert(Ptr{Int32}, v)::Ptr{Cvoid},t::Ptr{Cvoid},d::Int32,pr::Int32)::Cvoid
end

end
