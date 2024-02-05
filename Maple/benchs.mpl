restart:
kernelopts(numcpus = 1):
kernelopts(printbytes=false):
rep:="/Users/rouillie/Travail/Git/code/Data/Systems/":
read("/Users/rouillie/Travail/Git/code/Maple/zds.mpl"):
tab:=["caprasse","chandran9","reimer6","eco11","henrion6","chandran10","cp_d_3_n_6_p_6","cp_d_3_n_6_p_2","katsura11","eco12","fab_4","chandran11","eco13","katsura12","Noon7","Reimer7","phuoc1"];
for i from 1 to nops(tab) do
      fname:=cat(rep,tab[i],".mpl");
      read(fname):
      t11:=time():
      hh:=zds:-rur(numer(sys),vars):
      tt:=time()-t11:
      printf("%s : %g",tab[i],tt):
end:
