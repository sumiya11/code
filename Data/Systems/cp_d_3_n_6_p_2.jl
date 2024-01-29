R,(x1, x2, x3, x4, x5, x6)=polynomial_ring(QQ, ["x1","x2","x3","x4","x5","x6"],ordering=:degrevlex)

sys=[-7*x1^3+22*x1^2*x2-62*x1*x2^2+6*x2^3-55*x1^2*x3+97*x1*x2*x3+74*x2^2*x3+62*x1*x3^2+44*x2*x3^2+68*x3^3-94*x1^2*x4-73*x1*x2*x4+72*x2^2*x4-82*x1*x3*x4+29*x2*x3*x4-10*x3^2*x4-17*x1*x4^2-61*x2*x4^2+95*x3*x4^2-96*x4^3+87*x1^2*x5-4*x1*x2*x5+37*x2^2*x5+80*x1*x3*x5+98*x2*x3*x5+31*x3^2*x5-75*x1*x4*x5-8*x2*x4*x5+x3*x4*x5+72*x4^2*x5-40*x1*x5^2+11*x2*x5^2-28*x3*x5^2-90*x4*x5^2+53*x5^3-56*x1^2*x6-83*x1*x2*x6-23*x2^2*x6-44*x1*x3*x6-23*x2*x3*x6-51*x3^2*x6-10*x1*x4*x6-29*x2*x4*x6+x3*x4*x6-87*x4^2*x6+42*x1*x5*x6-49*x2*x5*x6+16*x3*x5*x6+43*x4*x5*x6-28*x5^2*x6+23*x1*x6^2+40*x2*x6^2-27*x3*x6^2-91*x4*x6^2+13*x5*x6^2+71*x6^3-10*x1*x2+87*x2^2+71*x1*x3+10*x2*x3+77*x3^2-7*x1*x4+95*x2*x4+55*x3*x4+47*x4^2-50*x1*x5-47*x2*x5+30*x3*x5+92*x4*x5+5*x5^2+75*x1*x6-81*x2*x6-15*x3*x6-88*x4*x6-10*x5*x6+16*x6^2-92*x1+91*x2-59*x3-48*x4-82*x5+83*x6+9,
-60*x1^3-83*x1^2*x2+5*x1*x2^2+33*x2^3+98*x1^2*x3+96*x1*x2*x3+10*x2^2*x3+98*x1*x3^2-25*x2*x3^2-89*x3^3-48*x1^2*x4-17*x1*x2*x4+7*x2^2*x4-64*x1*x3*x4-96*x2*x3*x4-77*x3^2*x4-34*x1*x4^2+7*x2*x4^2-42*x3*x4^2-69*x4^3-19*x1^2*x5+25*x1*x2*x5-89*x2^2*x5+64*x1*x3*x5+50*x2*x3*x5+69*x3^2*x5-13*x1*x4*x5-89*x2*x4*x5-33*x3*x4*x5-46*x4^2*x5+71*x1*x5^2-68*x2*x5^2+97*x3*x5^2-34*x4*x5^2-85*x5^3+62*x1^2*x6+91*x1*x2*x6+65*x2^2*x6-90*x1*x3*x6-60*x2*x3*x6+80*x3^2*x6+44*x1*x4*x6-70*x2*x4*x6+21*x3*x4*x6-33*x4^2*x6-47*x1*x5*x6-60*x2*x5*x6+30*x3*x5*x6+40*x4*x5*x6+54*x5^2*x6-53*x1*x6^2+52*x2*x6^2+89*x3*x6^2+x4*x6^2+18*x5*x6^2+91*x6^3+37*x1^2+12*x2^2-60*x1*x3-42*x2*x3+28*x3^2-2*x1*x4+34*x2*x4-35*x3*x4+87*x4^2-39*x1*x5+16*x2*x5-64*x3*x5+77*x4*x5-72*x1*x6-20*x2*x6-16*x3*x6-10*x4*x6+52*x5*x6-22*x6^2-97*x1-4*x2+59*x3-65*x4+36*x5+51*x6-27,
-2409*x1^4-1439*x1^3*x2+697*x1^2*x2^2-9855*x1*x2^3-7146*x2^4+29390*x1^3*x3-3736*x1^2*x2*x3+2706*x1*x2^2*x3-9612*x2^3*x3+21103*x1^2*x3^2+50947*x1*x2*x3^2-31872*x2^2*x3^2-33759*x1*x3^3-43596*x2*x3^3-6648*x3^4-16303*x1^3*x4+18574*x1^2*x2*x4+31574*x1*x2^2*x4-4195*x2^3*x4-18022*x1^2*x3*x4+43098*x1*x2*x3*x4-16618*x2^2*x3*x4+22663*x1*x3^2*x4-60197*x2*x3^2*x4+4565*x3^3*x4+4646*x1^2*x4^2-2844*x1*x2*x4^2-25519*x2^2*x4^2-24844*x1*x3*x4^2-27578*x2*x3*x4^2+9000*x3^2*x4^2+9159*x1*x4^3-1725*x2*x4^3+17436*x3*x4^3+1897*x4^4+9031*x1^3*x5-4849*x1^2*x2*x5+8532*x1*x2^2*x5+5110*x2^3*x5+15180*x1^2*x3*x5+24716*x1*x2*x3*x5+11030*x2^2*x3*x5+21226*x1*x3^2*x5+31088*x2*x3^2*x5-28744*x3^3*x5-7328*x1^2*x4*x5-1614*x1*x2*x4*x5+9699*x2^2*x4*x5+5631*x1*x3*x4*x5+7376*x2*x3*x4*x5+14727*x3^2*x4*x5-11031*x1*x4^2*x5+8649*x2*x4^2*x5-19127*x3*x4^2*x5+10797*x4^3*x5-5108*x1^2*x5^2+12230*x1*x2*x5^2+30804*x2^2*x5^2+22855*x1*x3*x5^2+41598*x2*x3*x5^2+24927*x3^2*x5^2-7122*x1*x4*x5^2+21334*x2*x4*x5^2-1799*x3*x4*x5^2+630*x4^2*x5^2+6456*x1*x5^3+9408*x2*x5^3+16640*x3*x5^3-3563*x4*x5^3-837*x5^4-8761*x1^3*x6-5782*x1^2*x2*x6-14184*x1*x2^2*x6-8883*x2^3*x6-42558*x1^2*x3*x6-60254*x1*x2*x3*x6-372*x2^2*x3*x6+26781*x1*x3^2*x6+14697*x2*x3^2*x6+22871*x3^3*x6+12449*x1^2*x4*x6+11970*x1*x2*x4*x6-2503*x2^2*x4*x6-4157*x1*x3*x4*x6+44850*x2*x3*x4*x6+20162*x3^2*x4*x6-4761*x1*x4^2*x6+1227*x2*x4^2*x6+1391*x3*x4^2*x6+6580*x4^3*x6-17246*x1^2*x5*x6-42853*x1*x2*x5*x6-20668*x2^2*x5*x6-19268*x1*x3*x5*x6-8874*x2*x3*x5*x6+48369*x3^2*x5*x6-615*x1*x4*x5*x6+17503*x2*x4*x5*x6+1559*x3*x4*x5*x6+6764*x4^2*x5*x6-8341*x1*x5^2*x6+5452*x2*x5^2*x6-9989*x3*x5^2*x6-1613*x4*x5^2*x6-5015*x5^3*x6+17971*x1^2*x6^2+4963*x1*x2*x6^2+6577*x2^2*x6^2+8049*x1*x3*x6^2+13036*x2*x3*x6^2-27847*x3^2*x6^2-7556*x1*x4*x6^2+6880*x2*x4*x6^2-17334*x3*x4*x6^2-12399*x4^2*x6^2-3457*x1*x5*x6^2-3216*x2*x5*x6^2-1322*x3*x5*x6^2-5206*x4*x5*x6^2+3969*x5^2*x6^2-6242*x1*x6^3-1788*x2*x6^3+8037*x3*x6^3-3683*x4*x6^3-5613*x5*x6^3+4964*x6^4+3593*x1^3+24048*x1^2*x2+11275*x1*x2^2-1782*x2^3-1912*x1^2*x3+12844*x1*x2*x3-28258*x2^2*x3-379*x1*x3^2-44548*x2*x3^2+12212*x3^3+21202*x1^2*x4-4994*x1*x2*x4-31229*x2^2*x4+6791*x1*x3*x4-34430*x2*x3*x4-18438*x3^2*x4+3781*x1*x4^2-23012*x2*x4^2-10609*x3*x4^2-5470*x4^3-4819*x1^2*x5+18404*x1*x2*x5-756*x2^2*x5-31168*x1*x3*x5+44950*x2*x3*x5+8991*x3^2*x5+18756*x1*x4*x5-1510*x2*x4*x5+26374*x3*x4*x5+6154*x4^2*x5-1584*x1*x5^2+14454*x2*x5^2-4368*x3*x5^2+18239*x4*x5^2-2775*x5^3-11216*x1^2*x6-21600*x1*x2*x6-7389*x2^2*x6-30878*x1*x3*x6+11434*x2*x3*x6+29896*x3^2*x6+1803*x1*x4*x6+1352*x2*x4*x6+42351*x3*x4*x6+13209*x4^2*x6+11271*x1*x5*x6+3782*x2*x5*x6-5904*x3*x5*x6+11760*x4*x5*x6-6343*x5^2*x6+2121*x1*x6^2+21020*x2*x6^2-21544*x3*x6^2+2846*x4*x6^2-10097*x5*x6^2-7609*x6^3+5699*x1^2-9326*x1*x2+561*x2^2+31541*x1*x3+11762*x2*x3-15332*x3^2-19226*x1*x4-11038*x2*x4-16003*x3*x4-11823*x4^2+9707*x1*x5-11236*x2*x5+17062*x3*x5-15057*x4*x5+7880*x5^2-1454*x1*x6+582*x2*x6+7009*x3*x6-1001*x4*x6+3139*x5*x6+14415*x6^2-5766*x1+7900*x2+3824*x3+4646*x4-7533*x5-7475*x6+5133,
-8858*x1^4+459*x1^3*x2+17410*x1^2*x2^2+5333*x1*x2^3-7002*x2^4-3846*x1^3*x3+9186*x1^2*x2*x3+11291*x1*x2^2*x3-5003*x2^3*x3-5322*x1^2*x3^2-12853*x1*x2*x3^2-13076*x2^2*x3^2-11375*x1*x3^3-14695*x2*x3^3-3638*x3^4-2412*x1^3*x4-6642*x1^2*x2*x4+913*x1*x2^2*x4+12330*x2^3*x4+3452*x1^2*x3*x4+4134*x1*x2*x3*x4-22925*x2^2*x3*x4-34507*x1*x3^2*x4-29614*x2*x3^2*x4-2139*x3^3*x4-20486*x1^2*x4^2+17684*x1*x2*x4^2+27579*x2^2*x4^2+16173*x1*x3*x4^2-45285*x2*x3*x4^2+4263*x3^2*x4^2+14601*x1*x4^3-25776*x2*x4^3-29857*x3*x4^3+14643*x4^4-3969*x1^3*x5-18651*x1^2*x2*x5+2227*x1*x2^2*x5+12524*x2^3*x5+7598*x1^2*x3*x5-21896*x1*x2*x3*x5-18561*x2^2*x3*x5-7358*x1*x3^2*x5-27356*x2*x3^2*x5-8473*x3^3*x5+2742*x1^2*x4*x5+1056*x1*x2*x4*x5-42944*x2^2*x4*x5-44063*x1*x3*x4*x5+13814*x2*x3*x4*x5-19315*x3^2*x4*x5+16028*x1*x4^2*x5-87331*x2*x4^2*x5+24858*x3*x4^2*x5-19372*x4^3*x5-13211*x1^2*x5^2-13791*x1*x2*x5^2+5261*x2^2*x5^2+1645*x1*x3*x5^2-12902*x2*x3*x5^2-8557*x3^2*x5^2-11911*x1*x4*x5^2+7046*x2*x4*x5^2-13493*x3*x4*x5^2-5605*x4^2*x5^2-2857*x1*x5^3-20059*x2*x5^3+873*x3*x5^3+1042*x4*x5^3-6494*x5^4+12676*x1^3*x6+13179*x1^2*x2*x6+13891*x1*x2^2*x6-8071*x2^3*x6+14011*x1^2*x3*x6+15036*x1*x2*x3*x6-4696*x2^2*x3*x6+7480*x1*x3^2*x6+8073*x2*x3^2*x6+2120*x3^3*x6-15726*x1^2*x4*x6+33888*x1*x2*x4*x6+26417*x2^2*x4*x6-5576*x1*x3*x4*x6-31756*x2*x3*x4*x6+8316*x3^2*x4*x6+42227*x1*x4^2*x6+35421*x2*x4^2*x6-16689*x3*x4^2*x6-8913*x4^3*x6+9139*x1^2*x5*x6+12007*x1*x2*x5*x6-4768*x2^2*x5*x6+1118*x1*x3*x5*x6+10781*x2*x3*x5*x6+8835*x3^2*x5*x6-7866*x1*x4*x5*x6-53192*x2*x4*x5*x6+34740*x3*x4*x5*x6-12088*x4^2*x5*x6+5718*x1*x5^2*x6+25017*x2*x5^2*x6-872*x3*x5^2*x6-1217*x4*x5^2*x6-370*x5^3*x6-7305*x1^2*x6^2+11627*x1*x2*x6^2+12553*x2^2*x6^2+7091*x1*x3*x6^2-4606*x2*x3*x6^2-5214*x3^2*x6^2+16764*x1*x4*x6^2+33978*x2*x4*x6^2-31408*x3*x4*x6^2-2994*x4^2*x6^2-4338*x1*x5*x6^2-25008*x2*x5*x6^2+3967*x3*x5*x6^2-24631*x4*x5*x6^2-2237*x5^2*x6^2+10478*x1*x6^3+10492*x2*x6^3-4695*x3*x6^3+9*x4*x6^3-6145*x5*x6^3+4772*x6^4-145*x1^3+3025*x1^2*x2-5785*x1*x2^2-9303*x2^3+485*x1^2*x3-13632*x1*x2*x3-17249*x2^2*x3-12252*x1*x3^2-15309*x2*x3^2-1355*x3^3+10973*x1^2*x4-33862*x1*x2*x4+973*x2^2*x4+6054*x1*x3*x4-6098*x2*x3*x4+14436*x3^2*x4-14167*x1*x4^2-2627*x2*x4^2-12786*x3*x4^2-21145*x4^3+13403*x1^2*x5-11978*x1*x2*x5-5071*x2^2*x5-1074*x1*x3*x5+12394*x2*x3*x5+2999*x3^2*x5-1735*x1*x4*x5+21238*x2*x4*x5+31459*x3*x4*x5+2334*x4^2*x5-955*x1*x5^2+29463*x2*x5^2+3716*x3*x5^2+17728*x4*x5^2+10141*x5^3-5153*x1^2*x6+168*x1*x2*x6-15993*x2^2*x6+8358*x1*x3*x6+9726*x2*x3*x6+7754*x3^2*x6-14186*x1*x4*x6-31314*x2*x4*x6+2362*x3*x4*x6+13413*x4^2*x6-15600*x1*x5*x6-11141*x2*x5*x6+16080*x3*x5*x6+12125*x4*x5*x6-5961*x5^2*x6+5348*x1*x6^2+15768*x2*x6^2-14803*x3*x6^2+1257*x4*x6^2-7465*x5*x6^2+2275*x6^3-10138*x1^2+6181*x1*x2+8143*x2^2-7813*x1*x3-20360*x2*x3-9147*x3^2-4087*x1*x4+20118*x2*x4-12884*x3*x4-2354*x4^2-587*x1*x5-13413*x2*x5-1570*x3*x5-17043*x4*x5-12524*x5^2+13849*x1*x6+2262*x2*x6+669*x3*x6-18349*x4*x6+1598*x5*x6-1327*x6^2+440*x1-6684*x2-5631*x3+11667*x4+11198*x5+3043*x6-6107,
6803*x1^4+1704*x1^3*x2-10902*x1^2*x2^2+11512*x1*x2^3-5265*x2^4-2147*x1^3*x3-1245*x1^2*x2*x3-22353*x1*x2^2*x3-22714*x2^3*x3+3958*x1^2*x3^2-4552*x1*x2*x3^2+622*x2^2*x3^2+8533*x1*x3^3+14242*x2*x3^3+3811*x3^4-3645*x1^3*x4-6107*x1^2*x2*x4+29089*x1*x2^2*x4-14144*x2^3*x4+9785*x1^2*x3*x4+2246*x1*x2*x3*x4-6906*x2^2*x3*x4-718*x1*x3^2*x4+11340*x2*x3^2*x4+3550*x3^3*x4+5188*x1^2*x4^2+9026*x1*x2*x4^2-15490*x2^2*x4^2-20989*x1*x3*x4^2-20099*x2*x3*x4^2-5511*x3^2*x4^2+5900*x1*x4^3-2147*x2*x4^3+7584*x3*x4^3+2302*x4^4-5615*x1^3*x5-3894*x1^2*x2*x5+27689*x1*x2^2*x5-4626*x2^3*x5+12606*x1^2*x3*x5+2792*x1*x2*x3*x5-960*x2^2*x3*x5+29663*x1*x3^2*x5+35022*x2*x3^2*x5+12348*x3^3*x5-18340*x1^2*x4*x5+27790*x1*x2*x4*x5-7301*x2^2*x4*x5+1199*x1*x3*x4*x5+18160*x2*x3*x4*x5-8319*x3^2*x4*x5-14385*x1*x4^2*x5+10282*x2*x4^2*x5-38449*x3*x4^2*x5+12184*x4^3*x5+14726*x1^2*x5^2+26295*x1*x2*x5^2-24942*x2^2*x5^2-15315*x1*x3*x5^2-43746*x2*x3*x5^2+17434*x3^2*x5^2+12591*x1*x4*x5^2-74495*x2*x4*x5^2+3374*x3*x4*x5^2+3356*x4^2*x5^2-6833*x1*x5^3+9432*x2*x5^3-34614*x3*x5^3+3203*x4*x5^3+8007*x5^4-3888*x1^3*x6-12126*x1^2*x2*x6+6316*x1*x2^2*x6+3055*x2^3*x6-13538*x1^2*x3*x6-40019*x1*x2*x3*x6-19717*x2^2*x3*x6-4864*x1*x3^2*x6-2219*x2*x3^2*x6+1993*x3^3*x6+23139*x1^2*x4*x6+8259*x1*x2*x4*x6-1186*x2^2*x4*x6+2694*x1*x3*x4*x6+11077*x2*x3*x4*x6+6229*x3^2*x4*x6-7223*x1*x4^2*x6+3938*x2*x4^2*x6+9751*x3*x4^2*x6+3633*x4^3*x6-1489*x1^2*x5*x6+4320*x1*x2*x5*x6+4303*x2^2*x5*x6-8782*x1*x3*x5*x6+24976*x2*x3*x5*x6-3851*x3^2*x5*x6+2486*x1*x4*x5*x6+58962*x2*x4*x5*x6-16395*x3*x4*x5*x6-6743*x4^2*x5*x6-1755*x1*x5^2*x6-6924*x2*x5^2*x6+17341*x3*x5^2*x6+8573*x4*x5^2*x6+19415*x5^3*x6-3730*x1^2*x6^2+4987*x1*x2*x6^2+2683*x2^2*x6^2-1447*x1*x3*x6^2-5712*x2*x3*x6^2+2535*x3^2*x6^2-643*x1*x4*x6^2-9854*x2*x4*x6^2+2308*x3*x4*x6^2-4923*x4^2*x6^2+10398*x1*x5*x6^2-626*x2*x5*x6^2+5432*x3*x5*x6^2+1221*x4*x5*x6^2-26038*x5^2*x6^2-6741*x1*x6^3-2370*x2*x6^3+734*x3*x6^3-248*x4*x6^3+7130*x5*x6^3+44*x6^4-4818*x1^3-3857*x1^2*x2+8070*x1*x2^2-11433*x2^3+4923*x1^2*x3+17726*x1*x2*x3+6198*x2^2*x3-8744*x1*x3^2+5335*x2*x3^2-74*x3^3+6694*x1^2*x4-14412*x1*x2*x4-29767*x2^2*x4-1702*x1*x3*x4-9702*x2*x3*x4+11925*x3^2*x4+447*x1*x4^2-8762*x2*x4^2+14154*x3*x4^2-12159*x4^3+317*x1^2*x5+16102*x1*x2*x5-28773*x2^2*x5-10944*x1*x3*x5+35068*x2*x3*x5-11673*x3^2*x5+12125*x1*x4*x5+754*x2*x4*x5+20717*x3*x4*x5+8172*x4^2*x5-6923*x1*x5^2-43386*x2*x5^2-3258*x3*x5^2-15562*x4*x5^2+10121*x5^3+11850*x1^2*x6-5796*x1*x2*x6+5985*x2^2*x6+3893*x1*x3*x6+3840*x2*x3*x6+1313*x3^2*x6-27838*x1*x4*x6-2685*x2*x4*x6+15458*x3*x4*x6+8609*x4^2*x6-14602*x1*x5*x6+32360*x2*x5*x6-8640*x3*x5*x6+12645*x4*x5*x6+20147*x5^2*x6+2101*x1*x6^2+8692*x2*x6^2-7300*x3*x6^2-5024*x4*x6^2-14590*x5*x6^2+1402*x6^3-7005*x1^2-7131*x1*x2+4727*x2^2+15658*x1*x3-1760*x2*x3+6557*x3^2-8280*x1*x4+12509*x2*x4-12293*x3*x4-1333*x4^2+17141*x1*x5-24460*x2*x5+28006*x3*x5-19925*x4*x5-27909*x5^2+2004*x1*x6+10400*x2*x6+2930*x3*x6-2089*x4*x6+836*x5*x6+2982*x6^2-4109*x1+9500*x2-8788*x3+13583*x4-340*x5+136*x6+2948,
-3284*x1^4-12015*x1^3*x2-4273*x1^2*x2^2+2025*x1*x2^3+3447*x2^4+5758*x1^3*x3+35462*x1^2*x2*x3+34047*x1*x2^2*x3+11277*x2^3*x3-5651*x1^2*x3^2-23533*x1*x2*x3^2+354*x2^2*x3^2+7596*x1*x3^3+9645*x2*x3^3+2245*x3^4-5340*x1^3*x4-7645*x1^2*x2*x4+19882*x1*x2^2*x4+11293*x2^3*x4+8017*x1^2*x3*x4-13592*x1*x2*x3*x4-18142*x2^2*x3*x4-9914*x1*x3^2*x4+7569*x2*x3^2*x4-1627*x3^3*x4-14719*x1^2*x4^2+11085*x1*x2*x4^2-5459*x2^2*x4^2+9749*x1*x3*x4^2-1127*x2*x3*x4^2-7445*x3^2*x4^2-1684*x1*x4^3+939*x2*x4^3-10597*x3*x4^3+2622*x4^4+3604*x1^3*x5-3648*x1^2*x2*x5-4799*x1*x2^2*x5+4487*x2^3*x5+3733*x1^2*x3*x5-13401*x1*x2*x3*x5-9958*x2^2*x3*x5-5309*x1*x3^2*x5-7633*x2*x3^2*x5+12110*x3^3*x5+3188*x1^2*x4*x5-14833*x1*x2*x4*x5-24400*x2^2*x4*x5+2010*x1*x3*x4*x5-2533*x2*x3*x4*x5+2070*x3^2*x4*x5+1449*x1*x4^2*x5-10788*x2*x4^2*x5+4383*x3*x4^2*x5-10220*x4^3*x5-5124*x1^2*x5^2-5596*x1*x2*x5^2-10267*x2^2*x5^2-3282*x1*x3*x5^2+7966*x2*x3*x5^2+1228*x3^2*x5^2-1735*x1*x4*x5^2+12159*x2*x4*x5^2+2131*x3*x4*x5^2-5870*x4^2*x5^2+2823*x1*x5^3-4980*x2*x5^3+8110*x3*x5^3+440*x4*x5^3-1310*x5^4+1436*x1^3*x6+26040*x1^2*x2*x6-16856*x1*x2^2*x6-6048*x2^3*x6-8576*x1^2*x3*x6-25872*x1*x2*x3*x6+25217*x2^2*x3*x6+16367*x1*x3^2*x6+36950*x2*x3^2*x6+1582*x3^3*x6-15002*x1^2*x4*x6-21292*x1*x2*x4*x6+35405*x2^2*x4*x6+1180*x1*x3*x4*x6+39924*x2*x3*x4*x6-10797*x3^2*x4*x6+11584*x1*x4^2*x6+8760*x2*x4^2*x6-32894*x3*x4^2*x6-3981*x4^3*x6-4095*x1^2*x5*x6-10094*x1*x2*x5*x6+24575*x2^2*x5*x6-9629*x1*x3*x5*x6+11100*x2*x3*x5*x6+15668*x3^2*x5*x6+1612*x1*x4*x5*x6-28570*x2*x4*x5*x6+7547*x3*x4*x5*x6-20345*x4^2*x5*x6+4057*x1*x5^2*x6+15032*x2*x5^2*x6-2918*x3*x5^2*x6-13234*x4*x5^2*x6-2162*x5^3*x6+33689*x1^2*x6^2-45042*x1*x2*x6^2-27561*x2^2*x6^2+59*x1*x3*x6^2+36180*x2*x3*x6^2+15855*x3^2*x6^2+8662*x1*x4*x6^2+61190*x2*x4*x6^2+9245*x3*x4*x6^2-27738*x4^2*x6^2-7881*x1*x5*x6^2+52932*x2*x5*x6^2+5242*x3*x5*x6^2+5895*x4*x5*x6^2+20899*x5^2*x6^2-48674*x1*x6^3-40248*x2*x6^3+16429*x3*x6^3+16537*x4*x6^3-509*x5*x6^3-156*x6^4+4021*x1^3+12237*x1^2*x2+11745*x1*x2^2+19521*x2^3-16613*x1^2*x3-20766*x1*x2*x3-10347*x2^2*x3-4953*x1*x3^2+8605*x2*x3^2-2421*x3^3+6361*x1^2*x4+10848*x1*x2*x4+2259*x2^2*x4+6259*x1*x3*x4-13376*x2*x3*x4+5042*x3^2*x4+7951*x1*x4^2-7739*x2*x4^2-9680*x3*x4^2+1049*x4^3-2821*x1^2*x5-7756*x1*x2*x5-25923*x2^2*x5+1565*x1*x3*x5+11498*x2*x3*x5-752*x3^2*x5-2676*x1*x4*x5-10177*x2*x4*x5+6270*x3*x4*x5-5573*x4^2*x5+5347*x1*x5^2+10012*x2*x5^2+2098*x3*x5^2-3886*x4*x5^2-2198*x5^3-4243*x1^2*x6-21890*x1*x2*x6+17941*x2^2*x6+7011*x1*x3*x6+32202*x2*x3*x6-9656*x3^2*x6+4514*x1*x4*x6+16992*x2*x4*x6+5601*x3*x4*x6+3901*x4^2*x6+12531*x1*x5*x6+820*x2*x5*x6-15984*x3*x5*x6-824*x4*x5*x6-8498*x5^2*x6+736*x1*x6^2+36842*x2*x6^2-750*x3*x6^2+22583*x4*x6^2-11959*x5*x6^2-21277*x6^3+14149*x1^2-13333*x1*x2-3012*x2^2-8797*x1*x3-5690*x2*x3+10605*x3^2-7638*x1*x4+922*x2*x4+6556*x3*x4-5001*x4^2-4724*x1*x5+24416*x2*x5+4734*x3*x5+17949*x4*x5+8723*x5^2-13476*x1*x6-11776*x2*x6+21689*x3*x6-2433*x4*x6+3005*x5*x6+27623*x6^2-6762*x1+4738*x2+2480*x3+761*x4+967*x5-6347*x6+4973,
11852*x1^4+12457*x1^3*x2-4694*x1^2*x2^2-6761*x1*x2^3-202*x2^4+24028*x1^3*x3+7378*x1^2*x2*x3-34402*x1*x2^2*x3-3178*x2^3*x3-21539*x1^2*x3^2-56788*x1*x2*x3^2+8056*x2^2*x3^2-42538*x1*x3^3-19117*x2*x3^3-18378*x3^4+4992*x1^3*x4-5840*x1^2*x2*x4+4911*x1*x2^2*x4+9371*x2^3*x4-29284*x1^2*x3*x4-19326*x1*x2*x3*x4-2036*x2^2*x3*x4-76280*x1*x3^2*x4-24417*x2*x3^2*x4+33594*x3^3*x4+34501*x1^2*x4^2-11304*x1*x2*x4^2-20055*x2^2*x4^2+36428*x1*x3*x4^2-43782*x2*x3*x4^2-95919*x3^2*x4^2-9346*x1*x4^3-37445*x2*x4^3-40212*x3*x4^3-31761*x4^4+10241*x1^3*x5+14926*x1^2*x2*x5-10141*x1*x2^2*x5-9420*x2^3*x5+24929*x1^2*x3*x5-20347*x1*x2*x3*x5-31044*x2^2*x3*x5-25137*x1*x3^2*x5-40146*x2*x3^2*x5-9859*x3^3*x5-19200*x1^2*x4*x5-22001*x1*x2*x4*x5-1742*x2^2*x4*x5-69326*x1*x3*x4*x5-376*x2*x3*x4*x5-11341*x3^2*x4*x5+13057*x1*x4^2*x5-7533*x2*x4^2*x5+54019*x3*x4^2*x5-12403*x4^3*x5+24912*x1^2*x5^2+8767*x1*x2*x5^2-17118*x2^2*x5^2+30010*x1*x3*x5^2-15315*x2*x3*x5^2-30024*x3^2*x5^2-16834*x1*x4*x5^2-14753*x2*x4*x5^2-54834*x3*x4*x5^2+31382*x4^2*x5^2+10679*x1*x5^3+4436*x2*x5^3+11139*x3*x5^3-14396*x4*x5^3+9682*x5^4-7788*x1^3*x6+1562*x1^2*x2*x6+1533*x1*x2^2*x6-731*x2^3*x6+21535*x1^2*x3*x6+19117*x1*x2*x3*x6-14352*x2^2*x3*x6+30850*x1*x3^2*x6-13802*x2*x3^2*x6+14005*x3^3*x6+18292*x1^2*x4*x6+3946*x1*x2*x4*x6-17105*x2^2*x4*x6+56672*x1*x3*x4*x6+13850*x2*x3*x4*x6-81887*x3^2*x4*x6-18130*x1*x4^2*x6-36429*x2*x4^2*x6+39681*x3*x4^2*x6-7737*x4^3*x6-6380*x1^2*x5*x6+1733*x1*x2*x5*x6-3361*x2^2*x5*x6+18560*x1*x3*x5*x6+16111*x2*x3*x5*x6+23079*x3^2*x5*x6+23544*x1*x4*x5*x6+21212*x2*x4*x5*x6+5042*x3*x4*x5*x6+2010*x4^2*x5*x6-4376*x1*x5^2*x6+741*x2*x5^2*x6+13171*x3*x5^2*x6+16249*x4*x5^2*x6-3135*x5^3*x6+15689*x1^2*x6^2+14647*x1*x2*x6^2-5743*x2^2*x6^2+23264*x1*x3*x6^2+6906*x2*x3*x6^2-23426*x3^2*x6^2-13546*x1*x4*x6^2-6610*x2*x4*x6^2+5896*x3*x4*x6^2+31082*x4^2*x6^2+16044*x1*x5*x6^2+9173*x2*x5*x6^2+2768*x3*x5*x6^2-10033*x4*x5*x6^2+17077*x5^2*x6^2-8532*x1*x6^3-1012*x2*x6^3+13802*x3*x6^3+19180*x4*x6^3-2161*x5*x6^3+8072*x6^4-8252*x1^3-20717*x1^2*x2-4319*x1*x2^2+4660*x2^3-13933*x1^2*x3-31295*x1*x2*x3+1906*x2^2*x3-28728*x1*x3^2+14373*x2*x3^2-3753*x3^3-31864*x1^2*x4-1508*x1*x2*x4+19963*x2^2*x4-1018*x1*x3*x4+50090*x2*x3*x4+41603*x3^2*x4-45643*x1*x4^2+3624*x2*x4^2-3739*x3*x4^2-987*x4^3-25842*x1^2*x5-19830*x1*x2*x5+6952*x2^2*x5-21215*x1*x3*x5-20376*x2*x3*x5+22424*x3^2*x5+1797*x1*x4*x5+19151*x2*x4*x5-2148*x3*x4*x5-10207*x4^2*x5-11997*x1*x5^2-14523*x2*x5^2-13399*x3*x5^2-9441*x4*x5^2-17860*x5^3+10372*x1^2*x6+6555*x1*x2*x6+4187*x2^2*x6+32357*x1*x3*x6-26843*x2*x3*x6-26593*x3^2*x6-16443*x1*x4*x6-24056*x2*x4*x6-42260*x3*x4*x6-17669*x4^2*x6+14997*x1*x5*x6+4320*x2*x5*x6-8303*x3*x5*x6-10532*x4*x5*x6+10310*x5^2*x6-13012*x1*x6^2-17649*x2*x6^2+16101*x3*x6^2-16150*x4*x6^2-13493*x5*x6^2+6631*x6^3+16095*x1^2+11413*x1*x2-4661*x2^2+10861*x1*x3-2291*x2*x3-29413*x3^2+25915*x1*x4+10762*x2*x4+9186*x3*x4+33874*x4^2+13543*x1*x5+13487*x2*x5+13658*x3*x5+13974*x4*x5+21990*x5^2-9538*x1*x6+1670*x2*x6+17805*x3*x6+10367*x4*x6-10112*x5*x6+10079*x6^2-6964*x1-10277*x2-8502*x3-21067*x4-14993*x5+5989*x6+6667,
-7481*x1^4-11178*x1^3*x2+1802*x1^2*x2^2-10295*x1*x2^3-6956*x2^4-30768*x1^3*x3-7264*x1^2*x2*x3-16910*x1*x2^2*x3-3262*x2^3*x3+4776*x1^2*x3^2+4373*x1*x2*x3^2+5819*x2^2*x3^2+36896*x1*x3^3+43988*x2*x3^3+22353*x3^4+15191*x1^3*x4+17113*x1^2*x2*x4+1930*x1*x2^2*x4-5535*x2^3*x4+28455*x1^2*x3*x4-7067*x1*x2*x3*x4+7652*x2^2*x3*x4-19599*x1*x3^2*x4-4077*x2*x3^2*x4-3071*x3^3*x4-6411*x1^2*x4^2-9958*x1*x2*x4^2-14374*x2^2*x4^2-18896*x1*x3*x4^2+8105*x2*x3*x4^2+18511*x3^2*x4^2+3995*x1*x4^3-3213*x2*x4^3+8915*x3*x4^3-1346*x4^4-7058*x1^3*x5+22822*x1^2*x2*x5-10834*x1*x2^2*x5-20856*x2^3*x5+14922*x1^2*x3*x5+13616*x1*x2*x3*x5-6576*x2^2*x3*x5+39104*x1*x3^2*x5-12810*x2*x3^2*x5+24624*x3^3*x5+11228*x1^2*x4*x5+15443*x1*x2*x4*x5-12254*x2^2*x4*x5+1740*x1*x3*x4*x5-13040*x2*x3*x4*x5-75528*x3^2*x4*x5-6590*x1*x4^2*x5-39709*x2*x4^2*x5-23070*x3*x4^2*x5-11690*x4^3*x5+7016*x1^2*x5^2-34683*x1*x2*x5^2-35985*x2^2*x5^2-33388*x1*x3*x5^2-15052*x2*x3*x5^2+5250*x3^2*x5^2+42307*x1*x4*x5^2+14063*x2*x4*x5^2+49383*x3*x4*x5^2-29251*x4^2*x5^2-26792*x1*x5^3-31266*x2*x5^3-37752*x3*x5^3+24356*x4*x5^3-8283*x5^4+7135*x1^3*x6+3708*x1^2*x2*x6+1417*x1*x2^2*x6+317*x2^3*x6-24876*x1^2*x3*x6+9540*x1*x2*x3*x6+2218*x2^2*x3*x6-17364*x1*x3^2*x6-42390*x2*x3^2*x6-1606*x3^3*x6-7896*x1^2*x4*x6+3309*x1*x2*x4*x6-3213*x2^2*x4*x6+15756*x1*x3*x4*x6+10899*x2*x3*x4*x6+24129*x3^2*x4*x6+6837*x1*x4^2*x6+2987*x2*x4^2*x6+2462*x3*x4^2*x6+4048*x4^3*x6-23262*x1^2*x5*x6+7900*x1*x2*x5*x6+7036*x2^2*x5*x6+398*x1*x3*x5*x6+12276*x2*x3*x5*x6-3922*x3^2*x5*x6-19997*x1*x4*x5*x6-13169*x2*x4*x5*x6+22868*x3*x4*x5*x6+10183*x4^2*x5*x6+39668*x1*x5^2*x6+32386*x2*x5^2*x6+17386*x3*x5^2*x6-6313*x4*x5^2*x6-6442*x5^3*x6-3646*x1^2*x6^2+2010*x1*x2*x6^2-1248*x2^2*x6^2-10970*x1*x3*x6^2+6392*x2*x3*x6^2-3099*x3^2*x6^2+7563*x1*x4*x6^2+7514*x2*x4*x6^2-8822*x3*x4*x6^2-3773*x4^2*x6^2-7910*x1*x5*x6^2-2506*x2*x5*x6^2-2988*x3*x5*x6^2+18937*x4*x5*x6^2-5623*x5^2*x6^2-2091*x1*x6^3+6347*x2*x6^3-6150*x3*x6^3-5162*x4*x6^3+1966*x5*x6^3-1643*x6^4+10916*x1^3+9742*x1^2*x2-339*x1*x2^2+2318*x2^3+7090*x1^2*x3+13680*x1*x2*x3-17140*x2^2*x3-22987*x1*x3^2-9213*x2*x3^2+3844*x3^3-16676*x1^2*x4-20458*x1*x2*x4-4096*x2^2*x4-2001*x1*x3*x4-2524*x2*x3*x4+45914*x3^2*x4-8517*x1*x4^2+8000*x2*x4^2-5088*x3*x4^2+11169*x4^3+9380*x1^2*x5-7776*x1*x2*x5+3080*x2^2*x5+39284*x1*x3*x5-14610*x2*x3*x5+31628*x3^2*x5-16585*x1*x4*x5-15761*x2*x4*x5+1936*x3*x4*x5-3279*x4^2*x5-4123*x1*x5^2+5067*x2*x5^2-48436*x3*x5^2-32770*x4*x5^2+1556*x5^3-3804*x1^2*x6-7378*x1*x2*x6+29*x2^2*x6+19722*x1*x3*x6+8530*x2*x3*x6+12851*x3^2*x6+4299*x1*x4*x6+3698*x2*x4*x6-19375*x3*x4*x6+8212*x4^2*x6+8752*x1*x5*x6+4046*x2*x5*x6+18966*x3*x5*x6+4056*x4*x5*x6+5239*x5^2*x6+5750*x1*x6^2+2797*x2*x6^2-2796*x3*x6^2-8472*x4*x6^2-902*x5*x6^2-576*x6^3-3725*x1^2+5951*x1*x2+4738*x2^2+6090*x1*x3-3948*x2*x3-31986*x3^2+4084*x1*x4+2764*x2*x4-7224*x3*x4+5897*x4^2+700*x1*x5+12246*x2*x5+4846*x3*x5+20510*x4*x5+13250*x5^2-5792*x1*x6-209*x2*x6+16742*x3*x6-312*x4*x6+1048*x5*x6+3557*x6^2+2887*x1-1255*x2+12142*x3-10861*x4-4758*x5-4330*x6+2714,
2078*x1^4+14519*x1^3*x2+20622*x1^2*x2^2+16077*x1*x2^3+5040*x2^4+27926*x1^3*x3+31256*x1^2*x2*x3+6594*x1*x2^2*x3+360*x2^3*x3-4242*x1^2*x3^2-3993*x1*x2*x3^2+7119*x2^2*x3^2-10192*x1*x3^3-13891*x2*x3^3+2703*x3^4-10108*x1^3*x4-4432*x1^2*x2*x4-12781*x1*x2^2*x4-5213*x2^3*x4+863*x1^2*x3*x4-15671*x1*x2*x3*x4-14856*x2^2*x3*x4-6086*x1*x3^2*x4-25043*x2*x3^2*x4-4903*x3^3*x4+9631*x1^2*x4^2+14510*x1*x2*x4^2-1177*x2^2*x4^2-1516*x1*x3*x4^2-16281*x2*x3*x4^2-24769*x3^2*x4^2+898*x1*x4^3-17177*x2*x4^3-10701*x3*x4^3-6789*x4^4+7013*x1^3*x5+20979*x1^2*x2*x5+14726*x1*x2^2*x5+3570*x2^3*x5-10090*x1^2*x3*x5+8650*x1*x2*x3*x5-3196*x2^2*x3*x5+12366*x1*x3^2*x5-12039*x2*x3^2*x5+22390*x3^3*x5+2502*x1^2*x4*x5+2625*x1*x2*x4*x5-10018*x2^2*x4*x5+6686*x1*x3*x4*x5+2581*x2*x3*x4*x5+21066*x3^2*x4*x5-587*x1*x4^2*x5-2381*x2*x4^2*x5+19358*x3*x4^2*x5+2702*x4^3*x5-2978*x1^2*x5^2+5059*x1*x2*x5^2+1257*x2^2*x5^2+11638*x1*x3*x5^2+12445*x2*x3*x5^2+5899*x3^2*x5^2-4695*x1*x4*x5^2+3744*x2*x4*x5^2-8973*x3*x4*x5^2+14776*x4^2*x5^2+3354*x1*x5^3+13125*x2*x5^3+4820*x3*x5^3-6161*x4*x5^3+1204*x5^4-6446*x1^3*x6-44518*x1^2*x2*x6-17899*x1*x2^2*x6+4021*x2^3*x6-24022*x1^2*x3*x6+16636*x1*x2*x3*x6+23914*x2^2*x3*x6+31424*x1*x3^2*x6+60440*x2*x3^2*x6+21894*x3^3*x6+27764*x1^2*x4*x6+16292*x1*x2*x4*x6+13082*x2^2*x4*x6+24184*x1*x3*x4*x6+13074*x2*x3*x4*x6-61213*x3^2*x4*x6-26074*x1*x4^2*x6-8096*x2*x4^2*x6+3860*x3*x4^2*x6-5660*x4^3*x6-7432*x1^2*x5*x6-2715*x1*x2*x5*x6+8766*x2^2*x5*x6+2218*x1*x3*x5*x6+33980*x2*x3*x5*x6+29964*x3^2*x5*x6+14117*x1*x4*x5*x6+17959*x2*x4*x5*x6+15960*x3*x4*x5*x6-273*x4^2*x5*x6-7186*x1*x5^2*x6-10856*x2*x5^2*x6-2130*x3*x5^2*x6+18484*x4*x5^2*x6-1826*x5^3*x6-23775*x1^2*x6^2+18785*x1*x2*x6^2+20772*x2^2*x6^2-10790*x1*x3*x6^2+7599*x2*x3*x6^2+105426*x3^2*x6^2-26592*x1*x4*x6^2+20294*x2*x4*x6^2+56914*x3*x4*x6^2+47339*x4^2*x6^2+3419*x1*x5*x6^2+22081*x2*x5*x6^2-18066*x3*x5*x6^2+7377*x4*x5*x6^2-27475*x5^2*x6^2+5926*x1*x6^3-3427*x2*x6^3-61926*x3*x6^3+11944*x4*x6^3-5308*x5*x6^3-26328*x6^4-2348*x1^3-5397*x1^2*x2+417*x1*x2^2-986*x2^3-17624*x1^2*x3+20972*x1*x2*x3+2888*x2^2*x3-2483*x1*x3^2-37159*x2*x3^2+7907*x3^3+23852*x1^2*x4+4456*x1*x2*x4-7364*x2^2*x4+31747*x1*x3*x4-32267*x2*x3*x4-21733*x3^2*x4-13995*x1*x4^2-22889*x2*x4^2-14522*x3*x4^2-9506*x4^3-14981*x1^2*x5-9706*x1*x2*x5+3858*x2^2*x5-18742*x1*x3*x5+4400*x2*x3*x5+11876*x3^2*x5+7316*x1*x4*x5-5490*x2*x4*x5+15089*x3*x4*x5-1247*x4^2*x5+2973*x1*x5^2+8441*x2*x5^2+18315*x3*x5^2+14480*x4*x5^2-658*x5^3+2610*x1^2*x6+2769*x1*x2*x6-4919*x2^2*x6-25886*x1*x3*x6+23776*x2*x3*x6+32020*x3^2*x6-19289*x1*x4*x6-91*x2*x4*x6+39308*x3*x4*x6-8155*x4^2*x6-6901*x1*x5*x6+4210*x2*x5*x6-4876*x3*x5*x6-4856*x4*x5*x6+746*x5^2*x6+34574*x1*x6^2+22077*x2*x6^2+27715*x3*x6^2+26914*x4*x6^2+19520*x5*x6^2-2347*x6^3-11905*x1^2-7343*x1*x2-3136*x2^2-19362*x1*x3+14201*x2*x3+29230*x3^2-8201*x1*x4+7557*x2*x4+12447*x3*x4+11781*x4^2+4795*x1*x5+1595*x2*x5-3878*x3*x5-5529*x4*x5-10093*x5^2+9842*x1*x6-7141*x2*x6-34366*x3*x6+6370*x4*x6-5544*x5*x6-36266*x6^2+8424*x1+9955*x2+5035*x3+11492*x4+4364*x5+1271*x6-7906,
5962*x1^4+324*x1^3*x2+6272*x1^2*x2^2+8954*x1*x2^3-6667*x2^4+4950*x1^3*x3+2187*x1^2*x2*x3+12071*x1*x2^2*x3+3885*x2^3*x3+1763*x1^2*x3^2+6640*x1*x2*x3^2+19348*x2^2*x3^2+1846*x1*x3^3+12023*x2*x3^3+1697*x3^4+4184*x1^3*x4+7634*x1^2*x2*x4+8498*x1*x2^2*x4+3988*x2^3*x4+6378*x1^2*x3*x4+2085*x1*x2*x3*x4-27006*x2^2*x3*x4+15767*x1*x3^2*x4+8293*x2*x3^2*x4+16121*x3^3*x4+26603*x1^2*x4^2+1672*x1*x2*x4^2+40445*x2^2*x4^2-1072*x1*x3*x4^2-2106*x2*x3*x4^2-13637*x3^2*x4^2-5321*x1*x4^3+28580*x2*x4^3+7019*x3*x4^3+28152*x4^4-14632*x1^3*x5+8082*x1^2*x2*x5+27686*x1*x2^2*x5-5941*x2^3*x5-38596*x1^2*x3*x5-2091*x1*x2*x3*x5+21982*x2^2*x3*x5-29140*x1*x3^2*x5+8795*x2*x3^2*x5-5160*x3^3*x5-7248*x1^2*x4*x5-4948*x1*x2*x4*x5+3236*x2^2*x4*x5+20488*x1*x3*x4*x5-49912*x2*x3*x4*x5+31764*x3^2*x4*x5-71770*x1*x4^2*x5+47762*x2*x4^2*x5-97834*x3*x4^2*x5-17676*x4^3*x5+24580*x1^2*x5^2+21162*x1*x2*x5^2-7159*x2^2*x5^2+10270*x1*x3*x5^2+755*x2*x3*x5^2+7983*x3^2*x5^2+33950*x1*x4*x5^2+3586*x2*x4*x5^2-15314*x3*x4*x5^2+86589*x4^2*x5^2+5692*x1*x5^3+29179*x2*x5^3-14372*x3*x5^3-22092*x4*x5^3+28356*x5^4+2796*x1^3*x6+14000*x1^2*x2*x6-1874*x1*x2^2*x6+1194*x2^3*x6-1516*x1^2*x3*x6+1576*x1*x2*x3*x6+148*x2^2*x3*x6-1402*x1*x3^2*x6-2606*x2*x3^2*x6+350*x3^3*x6+15236*x1^2*x4*x6-6580*x1*x2*x4*x6+30534*x2^2*x4*x6-9870*x1*x3*x4*x6-15238*x2*x3*x4*x6-59*x3^2*x4*x6+18398*x1*x4^2*x6+22987*x2*x4^2*x6+10134*x3*x4^2*x6+10137*x4^3*x6-10966*x1^2*x5*x6-8950*x1*x2*x5*x6+4464*x2^2*x5*x6-11245*x1*x3*x5*x6-14651*x2*x3*x5*x6-1737*x3^2*x5*x6-31772*x1*x4*x5*x6-11052*x2*x4*x5*x6-12132*x3*x4*x5*x6-37886*x4^2*x5*x6+1690*x1*x5^2*x6+9683*x2*x5^2*x6+3092*x3*x5^2*x6+67402*x4*x5^2*x6-28949*x5^3*x6-804*x1^2*x6^2+3695*x1*x2*x6^2+7577*x2^2*x6^2-8481*x1*x3*x6^2-1659*x2*x3*x6^2-5795*x3^2*x6^2+10188*x1*x4*x6^2+14785*x2*x4*x6^2+2487*x3*x4*x6^2-2501*x4^2*x6^2-16340*x1*x5*x6^2+5695*x2*x5*x6^2-15217*x3*x5*x6^2-12332*x4*x5*x6^2+28752*x5^2*x6^2+3483*x1*x6^3+5897*x2*x6^3-3001*x3*x6^3-5957*x4*x6^3-9518*x5*x6^3-1651*x6^4+1573*x1^3-6709*x1^2*x2-1217*x1*x2^2-8232*x2^3+9166*x1^2*x3-628*x1*x2*x3-11048*x2^2*x3+9607*x1*x3^2+8926*x2*x3^2+7830*x3^3-21879*x1^2*x4+2629*x1*x2*x4-19381*x2^2*x4-18831*x1*x3*x4-3151*x2*x3*x4-4014*x3^2*x4+16814*x1*x4^2-38811*x2*x4^2+43714*x3*x4^2-19984*x4^3-6846*x1^2*x5+14035*x1*x2*x5-29086*x2^2*x5+7409*x1*x3*x5+13285*x2*x3*x5+14367*x3^2*x5+17848*x1*x4*x5-19112*x2*x4*x5+12324*x3*x4*x5+36774*x4^2*x5+23267*x1*x5^2-45985*x2*x5^2+20810*x3*x5^2-46914*x4*x5^2-35363*x5^3+177*x1^2*x6-12238*x1*x2*x6+4228*x2^2*x6-10698*x1*x3*x6-622*x2*x3*x6-5536*x3^2*x6-15306*x1*x4*x6+5093*x2*x4*x6+26956*x3*x4*x6-23326*x4^2*x6-25869*x1*x5*x6+24567*x2*x5*x6-16564*x3*x5*x6+31308*x4*x5*x6+32858*x5^2*x6+7975*x1*x6^2+2441*x2*x6^2+5021*x3*x6^2-20467*x4*x6^2-6783*x5*x6^2-6176*x6^3-580*x1^2-7693*x1*x2+12961*x2^2-9459*x1*x3-10723*x2*x3-10441*x3^2-6372*x1*x4+14377*x2*x4-2180*x3*x4-29224*x4^2-15500*x1*x5+5123*x2*x5-23470*x3*x5-12536*x4*x5+15777*x5^2+10782*x1*x6-3687*x2*x6+9800*x3*x6-10029*x4*x6+1658*x5*x6-7889*x6^2-1794*x1+2385*x2+4132*x3+19936*x4+10276*x5-7134*x6-7058,
-8516*x1^4-18016*x1^3*x2-10412*x1^2*x2^2+1997*x1*x2^3+4841*x2^4-2320*x1^3*x3-5994*x1^2*x2*x3-14314*x1*x2^2*x3-4482*x2^3*x3-10336*x1^2*x3^2-17394*x1*x2*x3^2-252*x2^2*x3^2-12312*x1*x3^3-3747*x2*x3^3-4727*x3^4-10532*x1^3*x4-13712*x1^2*x2*x4-5859*x1*x2^2*x4-12445*x2^3*x4+970*x1^2*x3*x4+25074*x1*x2*x3*x4+14751*x2^2*x3*x4-29852*x1*x3^2*x4-23206*x2*x3^2*x4+10783*x3^3*x4-32698*x1^2*x4^2-47279*x1*x2*x4^2-16302*x2^2*x4^2+20824*x1*x3*x4^2-15102*x2*x3*x4^2-35892*x3^2*x4^2-19536*x1*x4^3+19401*x2*x4^3-19419*x3*x4^3-8505*x4^4+1056*x1^3*x5-5951*x1^2*x2*x5-13420*x1*x2^2*x5-6544*x2^3*x5+8882*x1^2*x3*x5+892*x1*x2*x3*x5-6657*x2^2*x3*x5-5937*x1*x3^2*x5-6765*x2*x3^2*x5-671*x3^3*x5+3104*x1^2*x4*x5+10764*x1*x2*x4*x5+15808*x2^2*x4*x5-24718*x1*x3*x4*x5-25974*x2*x3*x4*x5+16837*x3^2*x4*x5+30554*x1*x4^2*x5-18572*x2*x4^2*x5+6096*x3*x4^2*x5-15375*x4^3*x5-9833*x1^2*x5^2-7453*x1*x2*x5^2-6429*x2^2*x5^2-319*x1*x3*x5^2+3003*x2*x3*x5^2-11072*x3^2*x5^2-13385*x1*x4*x5^2-10523*x2*x4*x5^2+13303*x3*x4*x5^2-11620*x4^2*x5^2+1244*x1*x5^3+810*x2*x5^3-3026*x3*x5^3+3062*x4*x5^3-5812*x5^4+14016*x1^3*x6-392*x1^2*x2*x6-22273*x1*x2^2*x6+3433*x2^3*x6-3614*x1^2*x3*x6-13002*x1*x2*x3*x6+24568*x2^2*x3*x6-11172*x1*x3^2*x6-369*x2*x3^2*x6-4787*x3^3*x6-16676*x1^2*x4*x6-10360*x1*x2*x4*x6-25218*x2^2*x4*x6-25056*x1*x3*x4*x6-2972*x2*x3*x4*x6-2036*x3^2*x4*x6+23448*x1*x4^2*x6-5955*x2*x4^2*x6-79144*x3*x4^2*x6-38250*x4^3*x6+9940*x1^2*x5*x6+10288*x1*x2*x5*x6+10723*x2^2*x5*x6-17971*x1*x3*x5*x6-387*x2*x3*x5*x6+5212*x3^2*x5*x6-2054*x1*x4*x5*x6+8572*x2*x4*x5*x6+19520*x3*x4*x5*x6-23503*x4^2*x5*x6+5733*x1*x5^2*x6-8760*x2*x5^2*x6-15670*x3*x5^2*x6-10036*x4*x5^2*x6+1086*x5^3*x6-21988*x1^2*x6^2-22772*x1*x2*x6^2+14857*x2^2*x6^2-996*x1*x3*x6^2+23330*x2*x3*x6^2+7754*x3^2*x6^2+30676*x1*x4*x6^2-55503*x2*x4*x6^2+37138*x3*x4*x6^2-43803*x4^2*x6^2-21373*x1*x5*x6^2+24330*x2*x5*x6^2+13860*x3*x5*x6^2+58043*x4*x5*x6^2-21706*x5^2*x6^2-2502*x1*x6^3-2551*x2*x6^3-20344*x3*x6^3-33444*x4*x6^3-83*x5*x6^3-25056*x6^4+9822*x1^3+11514*x1^2*x2+5340*x1*x2^2+6084*x2^3+13480*x1^2*x3-2099*x1*x2*x3-11551*x2^2*x3-305*x1*x3^2-2712*x2*x3^2+1620*x3^3+19508*x1^2*x4+30088*x1*x2*x4+7918*x2^2*x4-11547*x1*x3*x4-21722*x2*x3*x4+6608*x3^2*x4+36550*x1*x4^2-10266*x2*x4^2-10849*x3*x4^2-3300*x4^3+11436*x1^2*x5+12924*x1*x2*x5+482*x2^2*x5-7681*x1*x3*x5-9810*x2*x3*x5+11696*x3^2*x5-13412*x1*x4*x5-23271*x2*x4*x5+8038*x3*x4*x5-26641*x4^2*x5-2992*x1*x5^2+2075*x2*x5^2+4170*x3*x5^2+14793*x4*x5^2+2104*x5^3-2090*x1^2*x6-6182*x1*x2*x6-7272*x2^2*x6+6735*x1*x3*x6+35207*x2*x3*x6+3553*x3^2*x6+2804*x1*x4*x6+5288*x2*x4*x6+15990*x3*x4*x6+49118*x4^2*x6-11398*x1*x5*x6+14106*x2*x5*x6+22130*x3*x5*x6-42*x4*x5*x6+3962*x5^2*x6+13812*x1*x6^2+15758*x2*x6^2+7021*x3*x6^2-2630*x4*x6^2-2087*x5*x6^2-17922*x6^3-6772*x1^2-21487*x1*x2-670*x2^2+1337*x1*x3+5887*x2*x3-2679*x3^2-18776*x1*x4+8347*x2*x4+13195*x3*x4+12794*x4^2-10543*x1*x5+16351*x2*x5+4583*x3*x5+28339*x4*x5-626*x5^2+11374*x1*x6+221*x2*x6-13788*x3*x6-25026*x4*x6-12353*x5*x6+209*x6^2+8140*x1-2282*x2+5503*x3-14888*x4-4845*x5+534*x6+2947,
4330*x1^4+7492*x1^3*x2+4239*x1^2*x2^2-3705*x1*x2^3+358*x2^4-122*x1^3*x3+17271*x1^2*x2*x3+13309*x1*x2^2*x3+3253*x2^3*x3+7393*x1^2*x3^2-445*x1*x2*x3^2-2707*x2^2*x3^2+9910*x1*x3^3+10117*x2*x3^3+5999*x3^4-1740*x1^3*x4-19951*x1^2*x2*x4-11546*x1*x2^2*x4-7738*x2^3*x4+10398*x1^2*x3*x4-684*x1*x2*x3*x4-6805*x2^2*x3*x4-4535*x1*x3^2*x4-4159*x2*x3^2*x4-1021*x3^3*x4-6066*x1^2*x4^2+8672*x1*x2*x4^2-7363*x2^2*x4^2-7424*x1*x3*x4^2-5368*x2*x3*x4^2+8448*x3^2*x4^2+4052*x1*x4^3-13853*x2*x4^3-1346*x3*x4^3-6378*x4^4-299*x1^3*x5-8759*x1^2*x2*x5-7756*x1*x2^2*x5-8279*x2^3*x5+17306*x1^2*x3*x5+2218*x1*x2*x3*x5-4522*x2^2*x3*x5+11439*x1*x3^2*x5+6307*x2*x3^2*x5+5240*x3^3*x5-8700*x1^2*x4*x5-5956*x1*x2*x4*x5-17322*x2^2*x4*x5+10607*x1*x3*x4*x5+20657*x2*x3*x4*x5-20407*x3^2*x4*x5+2501*x1*x4^2*x5-4997*x2*x4^2*x5+19369*x3*x4^2*x5+4882*x4^3*x5-2460*x1^2*x5^2+10224*x1*x2*x5^2-4008*x2^2*x5^2-29606*x1*x3*x5^2+6989*x2*x3*x5^2-1463*x3^2*x5^2+2042*x1*x4*x5^2-7253*x2*x4*x5^2-12170*x3*x4*x5^2-29108*x4^2*x5^2+2893*x1*x5^3-24655*x2*x5^3+11258*x3*x5^3+5701*x4*x5^3+1446*x5^4-8376*x1^3*x6+3365*x1^2*x2*x6-10034*x1*x2^2*x6+6403*x2^3*x6-140*x1^2*x3*x6-6735*x1*x2*x3*x6+11262*x2^2*x3*x6+12079*x1*x3^2*x6+10598*x2*x3^2*x6+12054*x3^3*x6+11548*x1^2*x4*x6+2898*x1*x2*x4*x6-4431*x2^2*x4*x6-1009*x1*x3*x4*x6+2931*x2*x3*x4*x6+16802*x3^2*x4*x6-11215*x1*x4^2*x6-10499*x2*x4^2*x6+7273*x3*x4^2*x6-6167*x4^3*x6+8150*x1^2*x5*x6-12682*x1*x2*x5*x6+15658*x2^2*x5*x6+1448*x1*x3*x5*x6-16700*x2*x3*x5*x6+858*x3^2*x5*x6+44146*x1*x4*x5*x6-28910*x2*x4*x5*x6-906*x3*x4*x5*x6+2296*x4^2*x5*x6-12648*x1*x5^2*x6+45590*x2*x5^2*x6+5768*x3*x5^2*x6-54246*x4*x5^2*x6+12354*x5^3*x6+27322*x1^2*x6^2+12342*x1*x2*x6^2+30021*x2^2*x6^2+9692*x1*x3*x6^2+3040*x2*x3*x6^2+192*x3^2*x6^2-31822*x1*x4*x6^2+6639*x2*x4*x6^2+22863*x3*x4*x6^2+37957*x4^2*x6^2-49751*x1*x5*x6^2+20408*x2*x5*x6^2-60848*x3*x5*x6^2-14858*x4*x5*x6^2+94104*x5^2*x6^2+19271*x1*x6^3-685*x2*x6^3+1264*x3*x6^3+6521*x4*x6^3-38292*x5*x6^3-285*x6^4-10123*x1^3-14671*x1^2*x2-980*x1*x2^2-10636*x2^3-11177*x1^2*x3-2328*x1*x2*x3+3829*x2^2*x3-19232*x1*x3^2-2065*x2*x3^2-325*x3^3+11259*x1^2*x4+25162*x1*x2*x4-3746*x2^2*x4+3017*x1*x3*x4-9188*x2*x3*x4+17232*x3^2*x4+947*x1*x4^2-14182*x2*x4^2-9459*x3*x4^2-1105*x4^3+4052*x1^2*x5+25296*x1*x2*x5-6168*x2^2*x5-2806*x1*x3*x5+12836*x2*x3*x5+8832*x3^2*x5+19885*x1*x4*x5-20417*x2*x4*x5+24904*x3*x4*x5-861*x4^2*x5+675*x1*x5^2-26741*x2*x5^2-7213*x3*x5^2-26546*x4*x5^2+6258*x5^3+6667*x1^2*x6+8329*x1*x2*x6-8282*x2^2*x6-18401*x1*x3*x6-1384*x2*x3*x6+358*x3^2*x6-27348*x1*x4*x6+9224*x2*x4*x6+12088*x3*x4*x6+20446*x4^2*x6-6938*x1*x5*x6+12184*x2*x5*x6+3296*x3*x5*x6+21490*x4*x5*x6+608*x5^2*x6-9305*x1*x6^2-16165*x2*x6^2+21248*x3*x6^2+16441*x4*x6^2+882*x5*x6^2-14954*x6^3+9471*x1^2-6728*x1*x2+7008*x2^2+10587*x1*x3+1804*x2*x3-10310*x3^2-21325*x1*x4+20038*x2*x4-5217*x3*x4+19184*x4^2-17234*x1*x5+16610*x2*x5-21234*x3*x5-2910*x4*x5+26374*x5^2+13347*x1*x6-2959*x2*x6-12658*x3*x6+3425*x4*x6-16148*x5*x6-32109*x6^2+3891*x1+831*x2+8694*x3+2289*x4-3394*x5-2370*x6-7170]
