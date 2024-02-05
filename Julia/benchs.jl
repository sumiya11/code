using Dates
import AbstractAlgebra,Primes,CPUTime,HostCPUFeatures
include("./Groebner.jl/src/Groebner.jl")
include("./code/Julia/rur.jl")

path_ex="./code/Data/Systems/";
l_ex=["caprasse","chandran9","reimer6","eco11","henrion6","chandran10","cp_d_3_n_6_p_6","cp_d_3_n_6_p_2","katsura11","eco12","fab_4","chandran11","eco13","katsura12","Noon7","Reimer
7","phuoc1"];
QQ=AbstractAlgebra.QQ;
polynomial_ring=AbstractAlgebra.polynomial_ring;

nn=Int32(28)
print("\nUsing ",nn," bit primes")
for i in 1:length(path_ex)
   print("\n\nProcess ",l_ex[i]," ",Dates.Time(Dates.now()))
   ex_name=string(path_ex,l_ex[i],".jl")
   include(ex_name);
   sys_z=convert_sys_to_sys_z(sys);
   CPUTime.CPUtic();
   dm,Dq,sys_T,_vars=prepare_system(sys_z,nn,R);
   qq=test_param(sys_T,nn);
   t1=CPUTime.CPUtoc();
   print("Time ",l_ex[i]," : ",t1)
   Base.flush(Base.stdout)
end
