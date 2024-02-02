R,(x1,x2,x3,x4)=polynomial_ring(QQ,["x1","x2","x3","x4"],ordering=:degrevlex)

sys=[-13872*x1^6-25721*x1^5*x2+16473*x1^5*x3-11560*x1^5*x4+17918*x1^4*x2^2+22542*x1^4*x2*x3-18207*x1^4*x2*x4-23409*x1^4*x3^2-15317*x1^4*x3*x4+27166*x1^4*x4^2+46400*x1^2*x2^4-372006*x1^2*x3^4-29095*x1^2*x4^4+100800*x1*x2^5+59200*x1*x2^4*x3-57600*x1*x2^4*x4+241115*x1*x2*x3^4+51313*x1*x2*x4^4+585565*x1*x3^5+13778*x1*x3^4*x4-37030*x1*x3*x4^4-49197*x1*x4^5+72000*x2^6-70400*x2^5*x3+73600*x2^4*x3^2+6400*x2^4*x3*x4+97600*x2^4*x4^2-337561*x2^2*x3^4+44436*x2^2*x4^4+89557*x2*x3^5-103335*x2*x3^4*x4+3174*x2*x3*x4^4+49726*x2*x4^5+158447*x3^6-454674*x3^5*x4+668233*x3^4*x4^2+46023*x3^2*x4^4-13225*x3*x4^5-32798*x4^6-22712*x1^5-79152*x1^4*x2+65263*x1^4*x3-29427*x1^4*x4+44268*x1^3*x2^2+55692*x1^3*x2*x3-44982*x1^3*x2*x4-57834*x1^3*x3^2-37842*x1^3*x3*x4+67116*x1^3*x4^2-211120*x1^2*x2^3-528876*x1^2*x3^3-15180*x1^2*x4^3-508240*x1*x2^4-269360*x1*x2^3*x3+262080*x1*x2^3*x4+342790*x1*x2*x3^3+26772*x1*x2*x4^3+1328498*x1*x3^4+19588*x1*x3^3*x4-19320*x1*x3*x4^3+19297*x1*x4^4-330800*x2^5+400320*x2^4*x3-158400*x2^4*x4-334880*x2^3*x3^2-29120*x2^3*x3*x4-444080*x2^3*x4^2-479906*x2^2*x3^3+23184*x2^2*x4^3+44654*x2*x3^4-146910*x2*x3^3*x4+1656*x2*x3*x4^3+30176*x2*x4^4+872828*x3^5-667071*x3^4*x4+950018*x3^3*x4^2+24012*x3^2*x4^3+44413*x3*x4^4-49381*x4^5-120193*x1^4-365275*x1^3*x2+269937*x1^3*x3-148982*x1^3*x4+254311*x1^2*x2^2+286338*x1^2*x2*x3-231273*x1^2*x2*x4-861813*x1^2*x3^2-194563*x1^2*x3*x4+292494*x1^2*x4^2+283703*x1*x2^3+34077*x1*x2^2*x3-33156*x1*x2^2*x4+365855*x1*x2*x3^2+92732*x1*x2*x4^2+1593673*x1*x3^3+20906*x1*x3^2*x4-66920*x1*x3*x4^2-65448*x1*x4^3+22405*x2^4-404524*x2^3*x3+720720*x2^3*x4-469831*x2^2*x3^2+3684*x2^2*x3*x4+136485*x2^2*x4^2+18361*x2*x3^3-156795*x2*x3^2*x4+5736*x2*x3*x4^2+92072*x2*x4^3+1016386*x3^4-719280*x3^3*x4+1097113*x3^2*x4^2+2872*x3*x4^3-53890*x4^4+23150*x1^3-67768*x1^2*x2+271841*x1^2*x3-183813*x1^2*x4+1273701*x1*x2^2+1104208*x1*x2*x3-830874*x1*x2*x4+850686*x1*x3^2-218358*x1*x3*x4+434000*x1*x4^2+904518*x2^3-933530*x2^2*x3-71019*x2^2*x4+709216*x2*x3^2-5924*x2*x3*x4+1051592*x2*x4^2+890896*x3^3-337575*x3^2*x4+567464*x3*x4^2-61604*x4^3+203345*x1^2-903977*x1*x2+1645515*x1*x3-690946*x1*x4+834765*x2^2+1134594*x2*x3-2213171*x2*x4-19958*x3^2-562481*x3*x4+1536474*x4^2+638674*x1-873870*x2+1290865*x3-884623*x4+659387, -10404*x1^6+18207*x1^5*x2+27455*x1^5*x3+23987*x1^5*x4+19363*x1^4*x2^2+21675*x1^4*x2*x3+4624*x1^4*x2*x4-4046*x1^4*x3^2-21964*x1^4*x3*x4+2312*x1^4*x4^2-134400*x1^2*x2^4-27556*x1^2*x3^4-19573*x1^2*x4^4+132800*x1*x2^5+116800*x1*x2^4*x3-11200*x1*x2^4*x4+68890*x1*x2*x3^4-26450*x1*x2*x4^4+654455*x1*x3^5+172225*x1*x3^4*x4-50255*x1*x3*x4^4+43378*x1*x4^5-56000*x2^6+94400*x2^5*x3-88000*x2^5*x4+148800*x2^4*x3^2+142400*x2^4*x3*x4+140800*x2^4*x4^2+516675*x2^2*x3^4+24334*x2^2*x4^4+620010*x2*x3^5+82668*x2*x3^4*x4-51842*x2*x3*x4^4-5819*x2*x4^5+661344*x3^6+151558*x3^5*x4-344450*x3^4*x4^2-39146*x3^2*x4^4+35972*x3*x4^5-43907*x4^6-5763*x1^5+18105*x1^4*x2+85459*x1^4*x3+74290*x1^4*x4+47838*x1^3*x2^2+53550*x1^3*x2*x3+11424*x1^3*x2*x4-9996*x1^3*x3^2-54264*x1^3*x3*x4+5712*x1^3*x4^2+611520*x1^2*x2^3-39176*x1^2*x3^3-10212*x1^2*x4^3-588240*x1*x2^4-531440*x1*x2^3*x3+50960*x1*x2^3*x4+97940*x1*x2*x3^3-13800*x1*x2*x4^3+413755*x1*x3^4+244850*x1*x3^3*x4-26220*x1*x3*x4^3+35857*x1*x4^4+266000*x2^5-323920*x2^4*x3+310800*x2^4*x4-677040*x2^3*x3^2-647920*x2^3*x3*x4-640640*x2^3*x4^2+734550*x2^2*x3^3+12696*x2^2*x4^3+302784*x2*x3^4+117528*x2*x3^3*x4-27048*x2*x3*x4^3-5152*x2*x4^4+1312230*x3^5-232317*x3^4*x4-489700*x3^3*x4^2-20424*x3^2*x4^3+53682*x3*x4^4-32959*x4^5-106588*x1^4+164871*x1^3*x2+392299*x1^3*x3+341821*x1^3*x4+168593*x1^2*x2^2+275325*x1^2*x2*x3+58736*x1^2*x2*x4-93206*x1^2*x3^2-278996*x1^2*x3*x4-6004*x1^2*x4^2+3643*x1*x2^3+67233*x1*x2^2*x3-6447*x1*x2^2*x4+104530*x1*x2*x3^2-47800*x1*x2*x4^2+258485*x1*x3^3+261325*x1*x3^2*x4-90820*x1*x3*x4^2+85292*x1*x4^3-57595*x2^4-426141*x2^3*x3+357025*x2^3*x4+869628*x2^2*x3^2+81969*x2^2*x3*x4+125024*x2^2*x4^2+118074*x2*x3^3+125436*x2*x3^2*x4-93688*x2*x3*x4^2-11620*x2*x4^3+1105246*x3^4-406644*x3^3*x4-593394*x3^2*x4^2+83224*x3*x4^3-89353*x4^4+51111*x1^3-1496529*x1^2*x2+583157*x1^2*x3+513182*x1^2*x4+1666292*x1*x2^2+1571122*x1*x2*x3-65368*x1*x2*x4-369015*x1*x3^2-202140*x1*x3*x4+75500*x1*x4^2-696073*x2^3+1420382*x2^2*x3-961456*x2^2*x4+1125180*x2*x3^2+1526168*x2*x3*x4+1467008*x2*x4^2+433010*x3^3-588173*x3^2*x4-168384*x3*x4^2-40568*x4^3-1083444*x1^2+1065097*x1*x2+1476517*x1*x3+980207*x1*x4+591079*x2^2+1984611*x2*x3-1242976*x2*x4+420084*x3^2-172896*x3*x4+682468*x4^2+257885*x1-661949*x2+923533*x3-129104*x4-717594, -14450*x1^6-3757*x1^5*x2-1445*x1^5*x3+11849*x1^5*x4+867*x1^4*x2^2+578*x1^4*x2*x3-289*x1^4*x2*x4+26877*x1^4*x3^2-25143*x1^4*x3*x4+16473*x1^4*x4^2+32000*x1^2*x2^4+413340*x1^2*x3^4-11638*x1^2*x4^4-40000*x1*x2^5+140800*x1*x2^4*x3+112000*x1*x2^4*x4+558009*x1*x2*x3^4-50784*x1*x2*x4^4-558009*x1*x3^5-303116*x1*x3^4*x4+4761*x1*x3*x4^4+49197*x1*x4^5+137600*x2^6+3200*x2^5*x3-105600*x2^5*x4+123200*x2^4*x3^2-158400*x2^4*x3*x4+20800*x2^4*x4^2-633788*x2^2*x3^4+8464*x2^2*x4^4+461563*x2*x3^5-606232*x2*x3^4*x4+45494*x2*x3*x4^4-40733*x2*x4^5-27556*x3^6-399562*x3^5*x4+585565*x3^4*x4^2+5290*x3^2*x4^4-17986*x3*x4^5+24863*x4^6-33099*x1^5-24021*x1^4*x2+11169*x1^4*x3+31586*x1^4*x4+2142*x1^3*x2^2+1428*x1^3*x2*x3-714*x1^3*x2*x4+66402*x1^3*x3^2-62118*x1^3*x3*x4+40698*x1^3*x4^2-145600*x1^2*x2^3+587640*x1^2*x3^3-6072*x1^2*x4^3+106800*x1*x2^4-640640*x1*x2^3*x3-509600*x1*x2^3*x4+793314*x1*x2*x3^3-26496*x1*x2*x4^3-641756*x1*x3^4-430936*x1*x3^3*x4+2484*x1*x3*x4^3+68517*x1*x4^4-734880*x2^5+86240*x2^4*x3+464480*x2^4*x4-560560*x2^3*x3^2+720720*x2^3*x3*x4-94640*x2^3*x4^2-901048*x2^2*x3^3+4416*x2^2*x4^3+1303764*x2*x3^4-861872*x2*x3^3*x4+23736*x2*x3*x4^3+22655*x2*x4^4-376737*x3^5-85822*x3^4*x4+832490*x3^3*x4^2+2760*x3^2*x4^3-39537*x3*x4^4-31464*x4^5-200533*x1^4-84137*x1^3*x2+18059*x1^3*x3+156223*x1^3*x4+29433*x1^2*x2^2+7342*x1^2*x2*x3-3671*x1^2*x2*x4+968583*x1^2*x3^2-319377*x1^2*x3*x4+188215*x1^2*x4^2+319135*x1*x2^3+81048*x1*x2^2*x3+64470*x1*x2^2*x4+846693*x1*x2*x3^2-91776*x1*x2*x4^2-631225*x1*x3^3-459932*x1*x3^2*x4+8604*x1*x3*x4^2+111264*x1*x4^3+471846*x2^4-456798*x2^3*x3+12014*x2^3*x4-890759*x2^2*x3^2-91179*x2^2*x3*x4+27269*x2^2*x4^2+1620987*x2*x3^3-919864*x2*x3^2*x4+82216*x2*x3*x4^2-50704*x2*x4^3-879946*x3^4+79306*x3^3*x4+898065*x3^2*x4^2-48236*x3*x4^3+32328*x4^4-224295*x1^3+95789*x1^2*x2+464631*x1^2*x3+187678*x1^2*x4-449917*x1*x2^2+1882888*x1*x2*x3+1145050*x1*x2*x4+199600*x1*x3^2-563034*x1*x3*x4+327186*x1*x4^2+1843276*x2^3-364441*x2^2*x3-1110474*x2^2*x4+2603922*x2*x3^2-2073144*x2*x3*x4+278540*x2*x4^2-1041309*x3^3+446662*x3^2*x4+358608*x3*x4^2-63504*x4^3-446371*x1^2-1214899*x1*x2+872945*x1*x3+973449*x1*x4-598445*x2^2+1708302*x2*x3-901201*x2*x4+701597*x3^2-1405783*x3*x4+792157*x4^2-568565*x1-1908427*x2+626559*x3+82240*x4-1356449, -8381*x1^6+15895*x1^5*x2-25721*x1^5*x3-26010*x1^5*x4-17918*x1^4*x2^2-21386*x1^4*x2*x3+8670*x1^4*x2*x4+12716*x1^4*x3^2-24854*x1^4*x3*x4-867*x1^4*x4^2-81600*x1^2*x2^4+137780*x1^2*x3^4-39146*x1^2*x4^4-134400*x1*x2^5+28800*x1*x2^4*x3+67200*x1*x2^4*x4+124002*x1*x2*x3^4-40204*x1*x2*x4^4+179114*x1*x3^5-447785*x1*x3^4*x4+36501*x1*x3*x4^4+20102*x1*x4^5-113600*x2^6-123200*x2^5*x3-38400*x2^5*x4+64000*x2^4*x3^2+131200*x2^4*x3*x4-59200*x2^4*x4^2+413340*x2^2*x3^4-24334*x2^2*x4^4-530453*x2*x3^5-289338*x2*x3^4*x4+16399*x2*x3*x4^4-28037*x2*x4^5-447785*x3^6-20667*x3^5*x4+489119*x3^4*x4^2+19573*x3^2*x4^4+11638*x3*x4^5+35972*x4^6-28509*x1^5+30600*x1^4*x2-77707*x1^4*x3-89403*x1^4*x4-44268*x1^3*x2^2-52836*x1^3*x2*x3+21420*x1^3*x2*x4+31416*x1^3*x3^2-61404*x1^3*x3*x4-2142*x1^3*x4^2+371280*x1^2*x2^3+195880*x1^2*x3^3-20424*x1^2*x4^3+768320*x1*x2^4-131040*x1*x2^3*x3-305760*x1*x2^3*x4+176292*x1*x2*x3^3-20976*x1*x2*x4^3+619761*x1*x3^4-636610*x1*x3^3*x4+19044*x1*x3*x4^3+49105*x1*x4^4+361680*x2^5+664560*x2^4*x3+72320*x2^4*x4-291200*x2^3*x3^2-596960*x2^3*x3*x4+269360*x2^3*x4^2+587640*x2^2*x3^3-12696*x2^2*x4^3-402799*x2*x3^4-411348*x2*x3^3*x4+8556*x2*x3*x4^3+29808*x2*x4^4-1029283*x3^5-608058*x3^4*x4+695374*x3^3*x4^2+10212*x3^2*x4^3+37812*x3*x4^4+31993*x4^5-97126*x1^4+180485*x1^3*x2-361705*x1^3*x3-392508*x1^3*x4-274573*x1^2*x2^2-271654*x1^2*x2*x3+110130*x1^2*x2*x4+370584*x1^2*x3^2-315706*x1^2*x3*x4-81757*x1^2*x4^2-790804*x1*x2^3+16578*x1*x2^2*x3+38682*x1*x2^2*x4+188154*x1*x2*x3^2-72656*x1*x2*x4^2+790860*x1*x3^3-679445*x1*x3^2*x4+65964*x1*x3*x4^2+56476*x1*x4^3+531969*x2^4-544117*x2^3*x3+443816*x2^3*x4+664020*x2^2*x3^2+75522*x2^2*x3*x4-78053*x2^2*x4^2-305387*x2*x3^3-439026*x2*x3^2*x4+29636*x2*x3*x4^2-27484*x2*x4^3-1409928*x3^4-854055*x3^3*x4+777535*x3^2*x4^2+37592*x3*x4^3+64502*x4^4-144141*x1^3-744624*x1^2*x2-435869*x1^2*x3-696237*x1^2*x4-1563618*x1*x2^2+95340*x1*x2*x3+804708*x1*x2*x4+858425*x1*x3^2-648720*x1*x3*x4+66938*x1*x4^2-783121*x2^3-932063*x2^2*x3-471840*x2^2*x4+821251*x2*x3^2+1172296*x2*x3*x4-551944*x2*x4^2-1162811*x3^3-884040*x3^2*x4+414516*x3*x4^2+36356*x4^3-432010*x1^2+1307963*x1*x2-510251*x1*x3-885832*x1*x4-2759850*x2^2-101890*x2*x3-1079130*x2*x4+91983*x3^2-480498*x3*x4-195183*x4^2+1103499*x1-2106786*x2-92513*x3-1468407*x4+268223]
