static void n1fv_2(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(4, is), MAKE_VOLATILE_STRIDE(4, os)) {
	       Var T1, T2;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	       ST(&(xo[WS(os, 1)]), VSUB(T1, T2), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[0]), VADD(T1, T2), ovs, &(xo[0]));
	  }
}
static void n1fv_3(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(6, is), MAKE_VOLATILE_STRIDE(6, os)) {
	       Var T1, T4, T6, T2, T3, T5;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
	       T3 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
	       T4 = VADD(T2, T3);
	       T6 = VBYI(VMUL(LDK(KP866025403), VSUB(T3, T2)));
	       ST(&(xo[0]), VADD(T1, T4), ovs, &(xo[0]));
	       T5 = VFNMS(LDK(KP500000000), T4, T1);
	       ST(&(xo[WS(os, 2)]), VSUB(T5, T6), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 1)]), VADD(T5, T6), ovs, &(xo[WS(os, 1)]));
	  }
}
static void n1fv_4(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(8, is), MAKE_VOLATILE_STRIDE(8, os)) {
	       Var T3, T7, T6, T8;
	       {
		    Var T1, T2, T4, T5;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T3 = VSUB(T1, T2);
		    T7 = VADD(T1, T2);
		    T4 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T6 = VBYI(VSUB(T4, T5));
		    T8 = VADD(T4, T5);
	       }
	       ST(&(xo[WS(os, 1)]), VSUB(T3, T6), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[0]), VADD(T7, T8), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 3)]), VADD(T3, T6), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 2)]), VSUB(T7, T8), ovs, &(xo[0]));
	  }
}
static void n1fv_5(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(10, is), MAKE_VOLATILE_STRIDE(10, os)) {
	       Var T8, T7, Td, T9, Tc;
	       T8 = LD(&(xi[0]), ivs, &(xi[0]));
	       {
		    Var T1, T2, T3, T4, T5, T6;
		    T1 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    T3 = VADD(T1, T2);
		    T4 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T6 = VADD(T4, T5);
		    T7 = VMUL(LDK(KP559016994), VSUB(T3, T6));
		    Td = VSUB(T4, T5);
		    T9 = VADD(T3, T6);
		    Tc = VSUB(T1, T2);
	       }
	       ST(&(xo[0]), VADD(T8, T9), ovs, &(xo[0]));
	       {
		    Var Te, Tf, Tb, Tg, Ta;
		    Te = VBYI(VFMA(LDK(KP951056516), Tc, VMUL(LDK(KP587785252), Td)));
		    Tf = VBYI(VFNMS(LDK(KP587785252), Tc, VMUL(LDK(KP951056516), Td)));
		    Ta = VFNMS(LDK(KP250000000), T9, T8);
		    Tb = VADD(T7, Ta);
		    Tg = VSUB(Ta, T7);
		    ST(&(xo[WS(os, 1)]), VSUB(Tb, Te), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VSUB(Tg, Tf), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 4)]), VADD(Te, Tb), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(Tf, Tg), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_6(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(12, is), MAKE_VOLATILE_STRIDE(12, os)) {
	       Var T3, Td, T6, Te, T9, Tf, Ta, Tg, T1, T2;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       T3 = VSUB(T1, T2);
	       Td = VADD(T1, T2);
	       {
		    Var T4, T5, T7, T8;
		    T4 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T5 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T6 = VSUB(T4, T5);
		    Te = VADD(T4, T5);
		    T7 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    T8 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T9 = VSUB(T7, T8);
		    Tf = VADD(T7, T8);
	       }
	       Ta = VADD(T6, T9);
	       Tg = VADD(Te, Tf);
	       ST(&(xo[WS(os, 3)]), VADD(T3, Ta), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[0]), VADD(Td, Tg), ovs, &(xo[0]));
	       {
		    Var Tb, Tc, Th, Ti;
		    Tb = VFNMS(LDK(KP500000000), Ta, T3);
		    Tc = VBYI(VMUL(LDK(KP866025403), VSUB(T9, T6)));
		    ST(&(xo[WS(os, 5)]), VSUB(Tb, Tc), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(Tb, Tc), ovs, &(xo[WS(os, 1)]));
		    Th = VFNMS(LDK(KP500000000), Tg, Td);
		    Ti = VBYI(VMUL(LDK(KP866025403), VSUB(Tf, Te)));
		    ST(&(xo[WS(os, 2)]), VSUB(Th, Ti), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(Th, Ti), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_7(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DVK(KP222520933, +0.222520933956314404288902564496794759466355569);
     DVK(KP623489801, +0.623489801858733530525004884004239810632274731);
     DVK(KP781831482, +0.781831482468029808708444526674057750232334519);
     DVK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DVK(KP433883739, +0.433883739117558120475768332848358754609990728);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(14, is), MAKE_VOLATILE_STRIDE(14, os)) {
	       Var T1, Ta, Td, T4, Tc, T7, Te, T8, T9, Tj, Ti;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T8 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
	       T9 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       Ta = VADD(T8, T9);
	       Td = VSUB(T9, T8);
	       {
		    Var T2, T3, T5, T6;
		    T2 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T3 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    Tc = VSUB(T3, T2);
		    T5 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T6 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T7 = VADD(T5, T6);
		    Te = VSUB(T6, T5);
	       }
	       ST(&(xo[0]), VADD(T1, VADD(T4, VADD(T7, Ta))), ovs, &(xo[0]));
	       Tj = VBYI(VFMA(LDK(KP433883739), Tc, VFNMS(LDK(KP781831482), Te, VMUL(LDK(KP974927912), Td))));
	       Ti = VFMA(LDK(KP623489801), T7, VFNMS(LDK(KP222520933), Ta, VFNMS(LDK(KP900968867), T4, T1)));
	       ST(&(xo[WS(os, 4)]), VSUB(Ti, Tj), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 3)]), VADD(Ti, Tj), ovs, &(xo[WS(os, 1)]));
	       {
		    Var Tf, Tb, Th, Tg;
		    Tf = VBYI(VFNMS(LDK(KP781831482), Td, VFNMS(LDK(KP433883739), Te, VMUL(LDK(KP974927912), Tc))));
		    Tb = VFMA(LDK(KP623489801), Ta, VFNMS(LDK(KP900968867), T7, VFNMS(LDK(KP222520933), T4, T1)));
		    ST(&(xo[WS(os, 5)]), VSUB(Tb, Tf), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 2)]), VADD(Tb, Tf), ovs, &(xo[0]));
		    Th = VBYI(VFMA(LDK(KP781831482), Tc, VFMA(LDK(KP974927912), Te, VMUL(LDK(KP433883739), Td))));
		    Tg = VFMA(LDK(KP623489801), T4, VFNMS(LDK(KP900968867), Ta, VFNMS(LDK(KP222520933), T7, T1)));
		    ST(&(xo[WS(os, 6)]), VSUB(Tg, Th), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 1)]), VADD(Tg, Th), ovs, &(xo[WS(os, 1)]));
	       }
	  }
}
static void n1fv_8(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(16, is), MAKE_VOLATILE_STRIDE(16, os)) {
	       Var T3, Tj, Tf, Tk, Ta, Tn, Tc, Tm;
	       {
		    Var T1, T2, Td, Te;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    T3 = VSUB(T1, T2);
		    Tj = VADD(T1, T2);
		    Td = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    Te = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    Tf = VSUB(Td, Te);
		    Tk = VADD(Td, Te);
		    {
			 Var T4, T5, T6, T7, T8, T9;
			 T4 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 T5 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 T6 = VSUB(T4, T5);
			 T7 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 T8 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T9 = VSUB(T7, T8);
			 Ta = VMUL(LDK(KP707106781), VADD(T6, T9));
			 Tn = VADD(T7, T8);
			 Tc = VMUL(LDK(KP707106781), VSUB(T9, T6));
			 Tm = VADD(T4, T5);
		    }
	       }
	       {
		    Var Tb, Tg, Tp, Tq;
		    Tb = VADD(T3, Ta);
		    Tg = VBYI(VSUB(Tc, Tf));
		    ST(&(xo[WS(os, 7)]), VSUB(Tb, Tg), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(Tb, Tg), ovs, &(xo[WS(os, 1)]));
		    Tp = VSUB(Tj, Tk);
		    Tq = VBYI(VSUB(Tn, Tm));
		    ST(&(xo[WS(os, 6)]), VSUB(Tp, Tq), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(Tp, Tq), ovs, &(xo[0]));
	       }
	       {
		    Var Th, Ti, Tl, To;
		    Th = VSUB(T3, Ta);
		    Ti = VBYI(VADD(Tf, Tc));
		    ST(&(xo[WS(os, 5)]), VSUB(Th, Ti), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VADD(Th, Ti), ovs, &(xo[WS(os, 1)]));
		    Tl = VADD(Tj, Tk);
		    To = VADD(Tm, Tn);
		    ST(&(xo[WS(os, 4)]), VSUB(Tl, To), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(Tl, To), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_9(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP342020143, +0.342020143325668733044099614682259580763083368);//0
     DVK(KP939692620, +0.939692620785908384054109277324731469936208134);//1
     DVK(KP642787609, +0.642787609686539326322643409907263432907559884);//2
     DVK(KP766044443, +0.766044443118978035202392650555416673935832457);//3
     DVK(KP984807753, +0.984807753012208059366743024589523013670643252);//4
     DVK(KP173648177, +0.173648177666930348851716626769314796000375677);//5
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);//6
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);//7
     DVK(KP852868531, +0.852868531952443209628250963940074071936020296);//d
     DVK(KP150383733, +0.150383733180435296639271897612501926072238258);//8
     DVK(KP296198132, +0.296198132726023843175338011893050938967728390);//a
     DVK(KP813797681, +0.813797681349373692844693217248393223289101568);//9
     DVK(KP556670399, +0.556670399226419366452912952047023132968291906);//c
     DVK(KP663413948, +0.663413948168938396205421319635891297216863310);//b
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(18, is), MAKE_VOLATILE_STRIDE(18, os)) {
	       Var T5, Ts, Tj, To, Tf, Tn, Tp, Tu, Tl, Ta, Tk, Tm, Tt;
	       {
		    Var T1, T2, T3, T4;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T3 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    T5 = VADD(T1, T4);
		    Ts = VMUL(LDK(KP866025403), VSUB(T3, T2));
		    Tj = VFNMS(LDK(KP500000000), T4, T1);
	       }
	       {
		    Var Tb, Te, Tc, Td;
		    Tb = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    Tc = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    Td = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    Te = VADD(Tc, Td);
		    To = VSUB(Td, Tc);
		    Tf = VADD(Tb, Te);
		    Tn = VFNMS(LDK(KP500000000), Te, Tb);
		    Tp = VFMA(LDK(KP173648177), Tn, VMUL(LDK(KP852868531), To));
		    Tu = VFNMS(LDK(KP984807753), Tn, VMUL(LDK(KP150383733), To));
	       }
	       {
		    Var T6, T9, T7, T8;
		    T6 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T7 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    T8 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
		    T9 = VADD(T7, T8);
		    Tl = VSUB(T8, T7);
		    Ta = VADD(T6, T9);
		    Tk = VFNMS(LDK(KP500000000), T9, T6);
		    Tm = VFMA(LDK(KP766044443), Tk, VMUL(LDK(KP556670399), Tl));
		    Tt = VFNMS(LDK(KP642787609), Tk, VMUL(LDK(KP663413948), Tl));
	       }
	       {
		    Var Ti, Tg, Th, Tz, TA;
		    Ti = VBYI(VMUL(LDK(KP866025403), VSUB(Tf, Ta)));
		    Tg = VADD(Ta, Tf);
		    Th = VFNMS(LDK(KP500000000), Tg, T5);
		    ST(&(xo[0]), VADD(T5, Tg), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 3)]), VADD(Th, Ti), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 6)]), VSUB(Th, Ti), ovs, &(xo[0]));
		    Tz = VFMA(LDK(KP173648177), Tk, VFNMS(LDK(KP296198132), To, VFNMS(LDK(KP939692620), Tn, VFNMS(LDK(KP852868531), Tl, Tj))));
		    TA = VBYI(VSUB(VFNMS(LDK(KP342020143), Tn, VFNMS(LDK(KP150383733), Tl, VFNMS(LDK(KP984807753), Tk, VMUL(LDK(KP813797681), To)))), Ts));
		    ST(&(xo[WS(os, 7)]), VSUB(Tz, TA), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 2)]), VADD(Tz, TA), ovs, &(xo[0]));
		    {
			 Var Tr, Tx, Tw, Ty, Tq, Tv;
			 Tq = VADD(Tm, Tp);
			 Tr = VADD(Tj, Tq);
			 Tx = VFMA(LDK(KP866025403), VSUB(Tt, Tu), VFNMS(LDK(KP500000000), Tq, Tj));
			 Tv = VADD(Tt, Tu);
			 Tw = VBYI(VADD(Ts, Tv));
			 Ty = VBYI(VADD(Ts, VFNMS(LDK(KP500000000), Tv, VMUL(LDK(KP866025403), VSUB(Tp, Tm)))));
			 ST(&(xo[WS(os, 8)]), VSUB(Tr, Tw), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 4)]), VADD(Tx, Ty), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 1)]), VADD(Tw, Tr), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 5)]), VSUB(Tx, Ty), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
static void n1fv_10(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(20, is), MAKE_VOLATILE_STRIDE(20, os)) {
	       Var Ti, Ty, Tm, Tn, Tw, Tt, Tz, TA, TB, T7, Te, Tj, Tg, Th;
	       Tg = LD(&(xi[0]), ivs, &(xi[0]));
	       Th = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
	       Ti = VSUB(Tg, Th);
	       Ty = VADD(Tg, Th);
	       {
		    Var T3, Tu, Td, Ts, T6, Tv, Ta, Tr;
		    {
			 Var T1, T2, Tb, Tc;
			 T1 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T2 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 T3 = VSUB(T1, T2);
			 Tu = VADD(T1, T2);
			 Tb = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Tc = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Td = VSUB(Tb, Tc);
			 Ts = VADD(Tb, Tc);
		    }
		    {
			 Var T4, T5, T8, T9;
			 T4 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T6 = VSUB(T4, T5);
			 Tv = VADD(T4, T5);
			 T8 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T9 = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 Ta = VSUB(T8, T9);
			 Tr = VADD(T8, T9);
		    }
		    Tm = VSUB(T3, T6);
		    Tn = VSUB(Ta, Td);
		    Tw = VSUB(Tu, Tv);
		    Tt = VSUB(Tr, Ts);
		    Tz = VADD(Tu, Tv);
		    TA = VADD(Tr, Ts);
		    TB = VADD(Tz, TA);
		    T7 = VADD(T3, T6);
		    Te = VADD(Ta, Td);
		    Tj = VADD(T7, Te);
	       }
	       ST(&(xo[WS(os, 5)]), VADD(Ti, Tj), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[0]), VADD(Ty, TB), ovs, &(xo[0]));
	       {
		    Var To, Tq, Tl, Tp, Tf, Tk;
		    To = VBYI(VFMA(LDK(KP951056516), Tm, VMUL(LDK(KP587785252), Tn)));
		    Tq = VBYI(VFNMS(LDK(KP587785252), Tm, VMUL(LDK(KP951056516), Tn)));
		    Tf = VMUL(LDK(KP559016994), VSUB(T7, Te));
		    Tk = VFNMS(LDK(KP250000000), Tj, Ti);
		    Tl = VADD(Tf, Tk);
		    Tp = VSUB(Tk, Tf);
		    ST(&(xo[WS(os, 1)]), VSUB(Tl, To), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 7)]), VADD(Tq, Tp), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 9)]), VADD(To, Tl), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VSUB(Tp, Tq), ovs, &(xo[WS(os, 1)]));
	       }
	       {
		    Var Tx, TF, TE, TG, TC, TD;
		    Tx = VBYI(VFNMS(LDK(KP587785252), Tw, VMUL(LDK(KP951056516), Tt)));
		    TF = VBYI(VFMA(LDK(KP951056516), Tw, VMUL(LDK(KP587785252), Tt)));
		    TC = VFNMS(LDK(KP250000000), TB, Ty);
		    TD = VMUL(LDK(KP559016994), VSUB(Tz, TA));
		    TE = VSUB(TC, TD);
		    TG = VADD(TD, TC);
		    ST(&(xo[WS(os, 2)]), VADD(Tx, TE), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 6)]), VSUB(TG, TF), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VSUB(TE, Tx), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(TF, TG), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_11(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP959492973, +0.959492973614497389890368057066327699062454848);//2
     DVK(KP654860733, +0.654860733945285064056925072466293553183791199);//0
     DVK(KP415415013, +0.415415013001886425529274149229623203524004910);//3
     DVK(KP142314838, +0.142314838273285140443792668616369668791051361);//1
     DVK(KP841253532, +0.841253532831181168861811648919367717513292498);//4
     DVK(KP281732556, +0.281732556841429697711417915346616899035777899);//7
     DVK(KP755749574, +0.755749574354258283774035843972344420179717445);//9
     DVK(KP909631995, +0.909631995354518371411715383079028460060241051);//6
     DVK(KP989821441, +0.989821441880932732376092037776718787376519372);//5
     DVK(KP540640817, +0.540640817455597582107635954318691695431770608);//8
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(22, is), MAKE_VOLATILE_STRIDE(22, os)) {
	       Var T1, T4, Ti, Tg, Tl, Td, Tk, Ta, Tj, T7, Tm, Tb, Tc, Tt, Ts;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       {
		    Var T2, T3, Te, Tf;
		    T2 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    T3 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    Ti = VSUB(T3, T2);
		    Te = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    Tf = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    Tg = VADD(Te, Tf);
		    Tl = VSUB(Tf, Te);
	       }
	       Tb = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
	       Tc = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
	       Td = VADD(Tb, Tc);
	       Tk = VSUB(Tc, Tb);
	       {
		    Var T8, T9, T5, T6;
		    T8 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    T9 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    Ta = VADD(T8, T9);
		    Tj = VSUB(T9, T8);
		    T5 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T6 = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    T7 = VADD(T5, T6);
		    Tm = VSUB(T6, T5);
	       }
	       ST(&(xo[0]), VADD(T1, VADD(T4, VADD(T7, VADD(Ta, VADD(Td, Tg))))), ovs, &(xo[0]));
	       {
		    Var Tn, Th, Tv, Tu;
		    Tn = VBYI(VFMA(LDK(KP755749574), Ti, VFMA(LDK(KP540640817), Tj, VFNMS(LDK(KP909631995), Tl, VFNMS(LDK(KP989821441), Tm, VMUL(LDK(KP281732556), Tk))))));
		    Th = VFMA(LDK(KP841253532), Ta, VFMA(LDK(KP415415013), Tg, VFNMS(LDK(KP959492973), Td, VFNMS(LDK(KP142314838), T7, VFNMS(LDK(KP654860733), T4, T1)))));
		    ST(&(xo[WS(os, 7)]), VSUB(Th, Tn), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 4)]), VADD(Th, Tn), ovs, &(xo[0]));
		    Tv = VBYI(VFMA(LDK(KP281732556), Ti, VFMA(LDK(KP755749574), Tj, VFNMS(LDK(KP909631995), Tk, VFNMS(LDK(KP540640817), Tm, VMUL(LDK(KP989821441), Tl))))));
		    Tu = VFMA(LDK(KP841253532), T7, VFMA(LDK(KP415415013), Td, VFNMS(LDK(KP142314838), Tg, VFNMS(LDK(KP654860733), Ta, VFNMS(LDK(KP959492973), T4, T1)))));
		    ST(&(xo[WS(os, 6)]), VSUB(Tu, Tv), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 5)]), VADD(Tu, Tv), ovs, &(xo[WS(os, 1)]));
	       }
	       Tt = VBYI(VFMA(LDK(KP989821441), Ti, VFMA(LDK(KP540640817), Tk, VFNMS(LDK(KP909631995), Tj, VFNMS(LDK(KP281732556), Tm, VMUL(LDK(KP755749574), Tl))))));
	       Ts = VFMA(LDK(KP415415013), Ta, VFMA(LDK(KP841253532), Td, VFNMS(LDK(KP654860733), Tg, VFNMS(LDK(KP959492973), T7, VFNMS(LDK(KP142314838), T4, T1)))));
	       ST(&(xo[WS(os, 8)]), VSUB(Ts, Tt), ovs, &(xo[0]));
	       ST(&(xo[WS(os, 3)]), VADD(Ts, Tt), ovs, &(xo[WS(os, 1)]));
	       {
		    Var Tr, Tq, Tp, To;
		    Tr = VBYI(VFMA(LDK(KP540640817), Ti, VFMA(LDK(KP909631995), Tm, VFMA(LDK(KP989821441), Tj, VFMA(LDK(KP755749574), Tk, VMUL(LDK(KP281732556), Tl))))));
		    Tq = VFMA(LDK(KP841253532), T4, VFMA(LDK(KP415415013), T7, VFNMS(LDK(KP959492973), Tg, VFNMS(LDK(KP654860733), Td, VFNMS(LDK(KP142314838), Ta, T1)))));
		    ST(&(xo[WS(os, 10)]), VSUB(Tq, Tr), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 1)]), VADD(Tq, Tr), ovs, &(xo[WS(os, 1)]));
		    Tp = VBYI(VFMA(LDK(KP909631995), Ti, VFNMS(LDK(KP540640817), Tl, VFNMS(LDK(KP989821441), Tk, VFNMS(LDK(KP281732556), Tj, VMUL(LDK(KP755749574), Tm))))));
		    To = VFMA(LDK(KP415415013), T4, VFMA(LDK(KP841253532), Tg, VFNMS(LDK(KP142314838), Td, VFNMS(LDK(KP959492973), Ta, VFNMS(LDK(KP654860733), T7, T1)))));
		    ST(&(xo[WS(os, 9)]), VSUB(To, Tp), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 2)]), VADD(To, Tp), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_12(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(24, is), MAKE_VOLATILE_STRIDE(24, os)) {
	       Var T5, Ta, TJ, Ty, Tq, Tp, Tg, Tl, TI, TA, Tz, Tu;
	       {
		    Var T1, T6, T4, Tw, T9, Tx;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T6 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
		    {
			 Var T2, T3, T7, T8;
			 T2 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T3 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T4 = VADD(T2, T3);
			 Tw = VSUB(T3, T2);
			 T7 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T9 = VADD(T7, T8);
			 Tx = VSUB(T8, T7);
		    }
		    T5 = VADD(T1, T4);
		    Ta = VADD(T6, T9);
		    TJ = VADD(Tw, Tx);
		    Ty = VMUL(LDK(KP866025403), VSUB(Tw, Tx));
		    Tq = VFNMS(LDK(KP500000000), T9, T6);
		    Tp = VFNMS(LDK(KP500000000), T4, T1);
	       }
	       {
		    Var Tc, Th, Tf, Ts, Tk, Tt;
		    Tc = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    Th = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
		    {
			 Var Td, Te, Ti, Tj;
			 Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 Te = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Tf = VADD(Td, Te);
			 Ts = VSUB(Te, Td);
			 Ti = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tj = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 Tk = VADD(Ti, Tj);
			 Tt = VSUB(Tj, Ti);
		    }
		    Tg = VADD(Tc, Tf);
		    Tl = VADD(Th, Tk);
		    TI = VADD(Ts, Tt);
		    TA = VFNMS(LDK(KP500000000), Tk, Th);
		    Tz = VFNMS(LDK(KP500000000), Tf, Tc);
		    Tu = VMUL(LDK(KP866025403), VSUB(Ts, Tt));
	       }
	       {
		    Var Tb, Tm, Tn, To;
		    Tb = VSUB(T5, Ta);
		    Tm = VBYI(VSUB(Tg, Tl));
		    ST(&(xo[WS(os, 9)]), VSUB(Tb, Tm), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VADD(Tb, Tm), ovs, &(xo[WS(os, 1)]));
		    Tn = VADD(T5, Ta);
		    To = VADD(Tg, Tl);
		    ST(&(xo[WS(os, 6)]), VSUB(Tn, To), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(Tn, To), ovs, &(xo[0]));
	       }
	       {
		    Var Tv, TE, TC, TD, Tr, TB;
		    Tr = VSUB(Tp, Tq);
		    Tv = VSUB(Tr, Tu);
		    TE = VADD(Tr, Tu);
		    TB = VSUB(Tz, TA);
		    TC = VBYI(VADD(Ty, TB));
		    TD = VBYI(VSUB(Ty, TB));
		    ST(&(xo[WS(os, 5)]), VSUB(Tv, TC), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 11)]), VSUB(TE, TD), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 7)]), VADD(TC, Tv), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(TD, TE), ovs, &(xo[WS(os, 1)]));
	       }
	       {
		    Var TK, TM, TH, TL, TF, TG;
		    TK = VBYI(VMUL(LDK(KP866025403), VSUB(TI, TJ)));
		    TM = VBYI(VMUL(LDK(KP866025403), VADD(TJ, TI)));
		    TF = VADD(Tp, Tq);
		    TG = VADD(Tz, TA);
		    TH = VSUB(TF, TG);
		    TL = VADD(TF, TG);
		    ST(&(xo[WS(os, 10)]), VSUB(TH, TK), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(TL, TM), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(TH, TK), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VSUB(TL, TM), ovs, &(xo[0]));
	       }
	  }
}
static void n1fv_13(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP174138601, +0.174138601152135905005660794929264742616964676);//0
     DVK(KP2_000000000, +2.000000000000000000000000000000000000000000000); //1
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);//2
     DVK(KP575140729, +0.575140729474003121368385547455453388461001608);//3
     DVK(KP1_732050807, +1.732050807568877293527446341505872366942805254);//4
     DVK(KP083333333, +0.083333333333333333333333333333333333333333333);//5
     DVK(KP251768516, +0.251768516431883313623436926934233488546674281);//6
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);//7
     DVK(KP503537032, +0.503537032863766627246873853868466977093348562);//8
     DVK(KP300462606, +0.300462606288665774426601772289207995520941381);//9
     DVK(KP258260390, +0.258260390311744861420450644284508567852516811);//a
     DVK(KP387390585, +0.387390585467617292130675966426762851778775217);//b
     DVK(KP075902986, +0.075902986037193865983102897245103540356428373);//c
     DVK(KP132983124, +0.132983124607418643793760531921092974399165133);//d
     DVK(KP300238635, +0.300238635966332641462884626667381504676006424);//e
     DVK(KP011599105, +0.011599105605768290721655456654083252189827041);//f
     DVK(KP156891391, +0.156891391051584611046832726756003269660212636);//g
     DVK(KP256247671, +0.256247671582936600958684654061725059144125175);//h
     DVK(KP113854479, +0.113854479055790798974654345867655310534642560);//i
     DVK(KP265966249, +0.265966249214837287587521063842185948798330267);//j

	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(26, is), MAKE_VOLATILE_STRIDE(26, os)) {
	       Var TW, Tb, Tm, Tu, TC, TR, TX, TK, TU, Tz, TB, TN, TT;
	       TW = LD(&(xi[0]), ivs, &(xi[0]));
	       {
		    Var T3, TH, Tl, Tw, Tp, Tg, Tv, To, T6, Tr, T9, Ts, Ta, TI, T1;
		    Var T2, Tq, Tt;
		    T1 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T3 = VSUB(T1, T2);
		    TH = VADD(T1, T2);
		    {
			 Var Th, Ti, Tj, Tk;
			 Th = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 Ti = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 Tj = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 Tk = VADD(Ti, Tj);
			 Tl = VADD(Th, Tk);
			 Tw = VSUB(Ti, Tj);
			 Tp = VFNMS(LDK(KP500000000), Tk, Th);
		    }
		    {
			 Var Tc, Td, Te, Tf;
			 Tc = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Td = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 Te = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 Tf = VADD(Td, Te);
			 Tg = VADD(Tc, Tf);
			 Tv = VSUB(Td, Te);
			 To = VFNMS(LDK(KP500000000), Tf, Tc);
		    }
		    {
			 Var T4, T5, T7, T8;
			 T4 = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 T5 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 T6 = VSUB(T4, T5);
			 Tr = VADD(T4, T5);
			 T7 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 T8 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T9 = VSUB(T7, T8);
			 Ts = VADD(T7, T8);
		    }
		    Ta = VADD(T6, T9);
		    TI = VADD(Tr, Ts);
		    Tb = VADD(T3, Ta);
		    Tm = VSUB(Tg, Tl);
		    Tq = VSUB(To, Tp);
		    Tt = VMUL(LDK(KP866025403), VSUB(Tr, Ts));
		    Tu = VADD(Tq, Tt);
		    TC = VSUB(Tq, Tt);
		    {
			 Var TP, TQ, TG, TJ;
			 TP = VADD(Tg, Tl);
			 TQ = VADD(TH, TI);
			 TR = VMUL(LDK(KP300462606), VSUB(TP, TQ));
			 TX = VADD(TP, TQ);
			 TG = VADD(To, Tp);
			 TJ = VFNMS(LDK(KP500000000), TI, TH);
			 TK = VSUB(TG, TJ);
			 TU = VADD(TG, TJ);
		    }
		    {
			 Var Tx, Ty, TL, TM;
			 Tx = VMUL(LDK(KP866025403), VSUB(Tv, Tw));
			 Ty = VFNMS(LDK(KP500000000), Ta, T3);
			 Tz = VSUB(Tx, Ty);
			 TB = VADD(Tx, Ty);
			 TL = VADD(Tv, Tw);
			 TM = VSUB(T6, T9);
			 TN = VSUB(TL, TM);
			 TT = VADD(TL, TM);
		    }
	       }
	       ST(&(xo[0]), VADD(TW, TX), ovs, &(xo[0]));
	       {
		    Var T19, T1n, T14, T13, T1f, T1k, Tn, TE, T1e, T1j, TS, T1m, TZ, T1c, TA;
		    Var TD;
		    {
			 Var T17, T18, T11, T12;
			 T17 = VFMA(LDK(KP387390585), TN, VMUL(LDK(KP265966249), TK));
			 T18 = VFNMS(LDK(KP503537032), TU, VMUL(LDK(KP113854479), TT));
			 T19 = VSUB(T17, T18);
			 T1n = VADD(T17, T18);
			 T14 = VFMA(LDK(KP575140729), Tm, VMUL(LDK(KP174138601), Tb));
			 T11 = VFNMS(LDK(KP156891391), TB, VMUL(LDK(KP256247671), TC));
			 T12 = VFMA(LDK(KP011599105), Tz, VMUL(LDK(KP300238635), Tu));
			 T13 = VSUB(T11, T12);
			 T1f = VADD(T14, T13);
			 T1k = VMUL(LDK(KP1_732050807), VADD(T11, T12));
		    }
		    Tn = VFNMS(LDK(KP174138601), Tm, VMUL(LDK(KP575140729), Tb));
		    TA = VFNMS(LDK(KP300238635), Tz, VMUL(LDK(KP011599105), Tu));
		    TD = VFMA(LDK(KP256247671), TB, VMUL(LDK(KP156891391), TC));
		    TE = VSUB(TA, TD);
		    T1e = VMUL(LDK(KP1_732050807), VADD(TD, TA));
		    T1j = VSUB(Tn, TE);
		    {
			 Var TO, T1b, TV, TY, T1a;
			 TO = VFNMS(LDK(KP132983124), TN, VMUL(LDK(KP258260390), TK));
			 T1b = VSUB(TR, TO);
			 TV = VFMA(LDK(KP251768516), TT, VMUL(LDK(KP075902986), TU));
			 TY = VFNMS(LDK(KP083333333), TX, TW);
			 T1a = VSUB(TY, TV);
			 TS = VFMA(LDK(KP2_000000000), TO, TR);
			 T1m = VADD(T1b, T1a);
			 TZ = VFMA(LDK(KP2_000000000), TV, TY);
			 T1c = VSUB(T1a, T1b);
		    }
		    {
			 Var TF, T10, T1l, T1o;
			 TF = VBYI(VFMA(LDK(KP2_000000000), TE, Tn));
			 T10 = VADD(TS, TZ);
			 ST(&(xo[WS(os, 1)]), VADD(TF, T10), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 12)]), VSUB(T10, TF), ovs, &(xo[0]));
			 {
			      Var T15, T16, T1p, T1q;
			      T15 = VBYI(VFMS(LDK(KP2_000000000), T13, T14));
			      T16 = VSUB(TZ, TS);
			      ST(&(xo[WS(os, 5)]), VADD(T15, T16), ovs, &(xo[WS(os, 1)]));
			      ST(&(xo[WS(os, 8)]), VSUB(T16, T15), ovs, &(xo[0]));
			      T1p = VADD(T1n, T1m);
			      T1q = VBYI(VADD(T1j, T1k));
			      ST(&(xo[WS(os, 4)]), VSUB(T1p, T1q), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 9)]), VADD(T1q, T1p), ovs, &(xo[WS(os, 1)]));
			 }
			 T1l = VBYI(VSUB(T1j, T1k));
			 T1o = VSUB(T1m, T1n);
			 ST(&(xo[WS(os, 3)]), VADD(T1l, T1o), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 10)]), VSUB(T1o, T1l), ovs, &(xo[0]));
			 {
			      Var T1h, T1i, T1d, T1g;
			      T1h = VBYI(VSUB(T1e, T1f));
			      T1i = VSUB(T1c, T19);
			      ST(&(xo[WS(os, 6)]), VADD(T1h, T1i), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 7)]), VSUB(T1i, T1h), ovs, &(xo[WS(os, 1)]));
			      T1d = VADD(T19, T1c);
			      T1g = VBYI(VADD(T1e, T1f));
			      ST(&(xo[WS(os, 2)]), VSUB(T1d, T1g), ovs, &(xo[0]));
			      ST(&(xo[WS(os, 11)]), VADD(T1g, T1d), ovs, &(xo[WS(os, 1)]));
			 }
		    }
	       }
	  }
}
static void n1fv_14(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP900968867, +0.900968867902419126236102319507445051165919162);//1
     DVK(KP222520933, +0.222520933956314404288902564496794759466355569);//0
     DVK(KP623489801, +0.623489801858733530525004884004239810632274731);//2
     DVK(KP781831482, +0.781831482468029808708444526674057750232334519);//4
     DVK(KP974927912, +0.974927912181823607018131682993931217232785801);//5
     DVK(KP433883739, +0.433883739117558120475768332848358754609990728);//3
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(28, is), MAKE_VOLATILE_STRIDE(28, os)) {
	       Var T3, Ty, To, TK, Tr, TE, Ta, TJ, Tq, TB, Th, TL, Ts, TH, T1;
	       Var T2;
	       T1 = LD(&(xi[0]), ivs, &(xi[0]));
	       T2 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
	       T3 = VSUB(T1, T2);
	       Ty = VADD(T1, T2);
	       {
		    Var Tk, TC, Tn, TD;
		    {
			 Var Ti, Tj, Tl, Tm;
			 Ti = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Tj = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 Tk = VSUB(Ti, Tj);
			 TC = VADD(Ti, Tj);
			 Tl = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 Tm = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tn = VSUB(Tl, Tm);
			 TD = VADD(Tl, Tm);
		    }
		    To = VADD(Tk, Tn);
		    TK = VSUB(TC, TD);
		    Tr = VSUB(Tn, Tk);
		    TE = VADD(TC, TD);
	       }
	       {
		    Var T6, Tz, T9, TA;
		    {
			 Var T4, T5, T7, T8;
			 T4 = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 T6 = VSUB(T4, T5);
			 Tz = VADD(T4, T5);
			 T7 = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 T9 = VSUB(T7, T8);
			 TA = VADD(T7, T8);
		    }
		    Ta = VADD(T6, T9);
		    TJ = VSUB(TA, Tz);
		    Tq = VSUB(T9, T6);
		    TB = VADD(Tz, TA);
	       }
	       {
		    Var Td, TF, Tg, TG;
		    {
			 Var Tb, Tc, Te, Tf;
			 Tb = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 Tc = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Td = VSUB(Tb, Tc);
			 TF = VADD(Tb, Tc);
			 Te = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 Tf = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 Tg = VSUB(Te, Tf);
			 TG = VADD(Te, Tf);
		    }
		    Th = VADD(Td, Tg);
		    TL = VSUB(TF, TG);
		    Ts = VSUB(Tg, Td);
		    TH = VADD(TF, TG);
	       }
	       ST(&(xo[WS(os, 7)]), VADD(T3, VADD(Ta, VADD(Th, To))), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[0]), VADD(Ty, VADD(TB, VADD(TH, TE))), ovs, &(xo[0]));
	       {
		    Var Tt, Tp, TP, TQ;
		    Tt = VBYI(VFNMS(LDK(KP781831482), Tr, VFNMS(LDK(KP433883739), Ts, VMUL(LDK(KP974927912), Tq))));
		    Tp = VFMA(LDK(KP623489801), To, VFNMS(LDK(KP900968867), Th, VFNMS(LDK(KP222520933), Ta, T3)));
		    ST(&(xo[WS(os, 5)]), VSUB(Tp, Tt), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 9)]), VADD(Tp, Tt), ovs, &(xo[WS(os, 1)]));
		    TP = VBYI(VFMA(LDK(KP974927912), TJ, VFMA(LDK(KP433883739), TL, VMUL(LDK(KP781831482), TK))));
		    TQ = VFMA(LDK(KP623489801), TE, VFNMS(LDK(KP900968867), TH, VFNMS(LDK(KP222520933), TB, Ty)));
		    ST(&(xo[WS(os, 2)]), VADD(TP, TQ), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 12)]), VSUB(TQ, TP), ovs, &(xo[0]));
	       }
	       {
		    Var Tv, Tu, TM, TI;
		    Tv = VBYI(VFMA(LDK(KP781831482), Tq, VFMA(LDK(KP974927912), Ts, VMUL(LDK(KP433883739), Tr))));
		    Tu = VFMA(LDK(KP623489801), Ta, VFNMS(LDK(KP900968867), To, VFNMS(LDK(KP222520933), Th, T3)));
		    ST(&(xo[WS(os, 13)]), VSUB(Tu, Tv), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(Tu, Tv), ovs, &(xo[WS(os, 1)]));
		    TM = VBYI(VFNMS(LDK(KP433883739), TK, VFNMS(LDK(KP974927912), TL, VMUL(LDK(KP781831482), TJ))));
		    TI = VFMA(LDK(KP623489801), TB, VFNMS(LDK(KP900968867), TE, VFNMS(LDK(KP222520933), TH, Ty)));
		    ST(&(xo[WS(os, 6)]), VSUB(TI, TM), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VADD(TM, TI), ovs, &(xo[0]));
	       }
	       {
		    Var TO, TN, Tx, Tw;
		    TO = VBYI(VFMA(LDK(KP433883739), TJ, VFNMS(LDK(KP974927912), TK, VMUL(LDK(KP781831482), TL))));
		    TN = VFMA(LDK(KP623489801), TH, VFNMS(LDK(KP222520933), TE, VFNMS(LDK(KP900968867), TB, Ty)));
		    ST(&(xo[WS(os, 4)]), VSUB(TN, TO), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 10)]), VADD(TO, TN), ovs, &(xo[0]));
		    Tx = VBYI(VFMA(LDK(KP433883739), Tq, VFNMS(LDK(KP781831482), Ts, VMUL(LDK(KP974927912), Tr))));
		    Tw = VFMA(LDK(KP623489801), Th, VFNMS(LDK(KP222520933), To, VFNMS(LDK(KP900968867), Ta, T3)));
		    ST(&(xo[WS(os, 11)]), VSUB(Tw, Tx), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VADD(Tw, Tx), ovs, &(xo[WS(os, 1)]));
	       }
	  }
}
static void n1fv_15(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP509036960, +0.509036960455127183450980863393907648510733164);//0
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);//9
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);//7
     DVK(KP823639103, +0.823639103546331925877420039278190003029660514);//1
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);//5
     DVK(KP216506350, +0.216506350946109661690930792688234045867850657);//2
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);//3
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);//4
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);//6
     DVK(KP484122918, +0.484122918275927110647408174972799951354115213);//8
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(30, is), MAKE_VOLATILE_STRIDE(30, os)) {
	       Var T5, T10, TB, TO, TU, TV, TR, Ta, Tf, Tg, Tl, Tq, Tr, TE, TH;
	       Var TI, TZ, T11, T1f, T1g;
	       {
		    Var T1, T2, T3, T4;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    T3 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    T4 = VADD(T2, T3);
		    T5 = VADD(T1, T4);
		    T10 = VSUB(T3, T2);
		    TB = VFNMS(LDK(KP500000000), T4, T1);
	       }
	       {
		    Var T6, T9, TC, TP, Tm, Tp, TG, TN, Tb, Te, TD, TQ, Th, Tk, TF;
		    Var TM, TX, TY;
		    {
			 Var T7, T8, Tn, To;
			 T6 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T7 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 T9 = VADD(T7, T8);
			 TC = VFNMS(LDK(KP500000000), T9, T6);
			 TP = VSUB(T8, T7);
			 Tm = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 Tn = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 To = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 Tp = VADD(Tn, To);
			 TG = VFNMS(LDK(KP500000000), Tp, Tm);
			 TN = VSUB(To, Tn);
		    }
		    {
			 Var Tc, Td, Ti, Tj;
			 Tb = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 Tc = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 Td = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 Te = VADD(Tc, Td);
			 TD = VFNMS(LDK(KP500000000), Te, Tb);
			 TQ = VSUB(Td, Tc);
			 Th = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Ti = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Tj = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tk = VADD(Ti, Tj);
			 TF = VFNMS(LDK(KP500000000), Tk, Th);
			 TM = VSUB(Tj, Ti);
		    }
		    TO = VSUB(TM, TN);
		    TU = VSUB(TF, TG);
		    TV = VSUB(TC, TD);
		    TR = VSUB(TP, TQ);
		    Ta = VADD(T6, T9);
		    Tf = VADD(Tb, Te);
		    Tg = VADD(Ta, Tf);
		    Tl = VADD(Th, Tk);
		    Tq = VADD(Tm, Tp);
		    Tr = VADD(Tl, Tq);
		    TE = VADD(TC, TD);
		    TH = VADD(TF, TG);
		    TI = VADD(TE, TH);
		    TX = VADD(TP, TQ);
		    TY = VADD(TM, TN);
		    TZ = VMUL(LDK(KP484122918), VSUB(TX, TY));
		    T11 = VADD(TX, TY);
	       }
	       T1f = VADD(TB, TI);
	       T1g = VBYI(VMUL(LDK(KP866025403), VADD(T10, T11)));
	       ST(&(xo[WS(os, 5)]), VSUB(T1f, T1g), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 10)]), VADD(T1f, T1g), ovs, &(xo[0]));
	       {
		    Var Tu, Ts, Tt, Ty, TA, Tw, Tx, Tz, Tv;
		    Tu = VMUL(LDK(KP559016994), VSUB(Tg, Tr));
		    Ts = VADD(Tg, Tr);
		    Tt = VFNMS(LDK(KP250000000), Ts, T5);
		    Tw = VSUB(Tl, Tq);
		    Tx = VSUB(Ta, Tf);
		    Ty = VBYI(VFNMS(LDK(KP587785252), Tx, VMUL(LDK(KP951056516), Tw)));
		    TA = VBYI(VFMA(LDK(KP951056516), Tx, VMUL(LDK(KP587785252), Tw)));
		    ST(&(xo[0]), VADD(T5, Ts), ovs, &(xo[0]));
		    Tz = VADD(Tu, Tt);
		    ST(&(xo[WS(os, 6)]), VSUB(Tz, TA), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 9)]), VADD(TA, Tz), ovs, &(xo[WS(os, 1)]));
		    Tv = VSUB(Tt, Tu);
		    ST(&(xo[WS(os, 3)]), VSUB(Tv, Ty), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 12)]), VADD(Ty, Tv), ovs, &(xo[0]));
	       }
	       {
		    Var TS, TW, T1b, T18, T13, T1a, TL, T17, T12, TJ, TK;
		    TS = VFNMS(LDK(KP509036960), TR, VMUL(LDK(KP823639103), TO));
		    TW = VFNMS(LDK(KP587785252), TV, VMUL(LDK(KP951056516), TU));
		    T1b = VFMA(LDK(KP951056516), TV, VMUL(LDK(KP587785252), TU));
		    T18 = VFMA(LDK(KP823639103), TR, VMUL(LDK(KP509036960), TO));
		    T12 = VFNMS(LDK(KP216506350), T11, VMUL(LDK(KP866025403), T10));
		    T13 = VSUB(TZ, T12);
		    T1a = VADD(TZ, T12);
		    TJ = VFNMS(LDK(KP250000000), TI, TB);
		    TK = VMUL(LDK(KP559016994), VSUB(TE, TH));
		    TL = VSUB(TJ, TK);
		    T17 = VADD(TK, TJ);
		    {
			 Var TT, T14, T1d, T1e;
			 TT = VSUB(TL, TS);
			 T14 = VBYI(VSUB(TW, T13));
			 ST(&(xo[WS(os, 8)]), VSUB(TT, T14), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 7)]), VADD(TT, T14), ovs, &(xo[WS(os, 1)]));
			 T1d = VSUB(T17, T18);
			 T1e = VBYI(VADD(T1b, T1a));
			 ST(&(xo[WS(os, 11)]), VSUB(T1d, T1e), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 4)]), VADD(T1d, T1e), ovs, &(xo[0]));
		    }
		    {
			 Var T15, T16, T19, T1c;
			 T15 = VADD(TL, TS);
			 T16 = VBYI(VADD(TW, T13));
			 ST(&(xo[WS(os, 13)]), VSUB(T15, T16), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 2)]), VADD(T15, T16), ovs, &(xo[0]));
			 T19 = VADD(T17, T18);
			 T1c = VBYI(VSUB(T1a, T1b));
			 ST(&(xo[WS(os, 14)]), VSUB(T19, T1c), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 1)]), VADD(T19, T1c), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
static void n1fv_16(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(32, is), MAKE_VOLATILE_STRIDE(32, os)) {
	       Var Tp, T13, Tu, TN, Tm, T14, Tv, TY, T7, T17, Ty, TT, Te, T16, Tx;
	       Var TQ;
	       {
		    Var Tn, To, TM, Ts, Tt, TL;
		    Tn = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    To = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
		    TM = VADD(Tn, To);
		    Ts = LD(&(xi[0]), ivs, &(xi[0]));
		    Tt = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
		    TL = VADD(Ts, Tt);
		    Tp = VSUB(Tn, To);
		    T13 = VADD(TL, TM);
		    Tu = VSUB(Ts, Tt);
		    TN = VSUB(TL, TM);
	       }
	       {
		    Var Ti, TW, Tl, TX;
		    {
			 Var Tg, Th, Tj, Tk;
			 Tg = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 Th = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Ti = VSUB(Tg, Th);
			 TW = VADD(Tg, Th);
			 Tj = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 Tk = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 Tl = VSUB(Tj, Tk);
			 TX = VADD(Tj, Tk);
		    }
		    Tm = VMUL(LDK(KP707106781), VSUB(Ti, Tl));
		    T14 = VADD(TX, TW);
		    Tv = VMUL(LDK(KP707106781), VADD(Tl, Ti));
		    TY = VSUB(TW, TX);
	       }
	       {
		    Var T3, TR, T6, TS;
		    {
			 Var T1, T2, T4, T5;
			 T1 = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
			 T2 = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 T3 = VSUB(T1, T2);
			 TR = VADD(T1, T2);
			 T4 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T5 = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 T6 = VSUB(T4, T5);
			 TS = VADD(T4, T5);
		    }
		    T7 = VFNMS(LDK(KP923879532), T6, VMUL(LDK(KP382683432), T3));
		    T17 = VADD(TR, TS);
		    Ty = VFMA(LDK(KP923879532), T3, VMUL(LDK(KP382683432), T6));
		    TT = VSUB(TR, TS);
	       }
	       {
		    Var Ta, TO, Td, TP;
		    {
			 Var T8, T9, Tb, Tc;
			 T8 = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 T9 = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 Ta = VSUB(T8, T9);
			 TO = VADD(T8, T9);
			 Tb = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 Tc = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 Td = VSUB(Tb, Tc);
			 TP = VADD(Tb, Tc);
		    }
		    Te = VFMA(LDK(KP382683432), Ta, VMUL(LDK(KP923879532), Td));
		    T16 = VADD(TO, TP);
		    Tx = VFNMS(LDK(KP382683432), Td, VMUL(LDK(KP923879532), Ta));
		    TQ = VSUB(TO, TP);
	       }
	       {
		    Var T15, T18, T19, T1a;
		    T15 = VADD(T13, T14);
		    T18 = VADD(T16, T17);
		    ST(&(xo[WS(os, 8)]), VSUB(T15, T18), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(T15, T18), ovs, &(xo[0]));
		    T19 = VSUB(T13, T14);
		    T1a = VBYI(VSUB(T17, T16));
		    ST(&(xo[WS(os, 12)]), VSUB(T19, T1a), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(T19, T1a), ovs, &(xo[0]));
	       }
	       {
		    Var TV, T11, T10, T12, TU, TZ;
		    TU = VMUL(LDK(KP707106781), VADD(TQ, TT));
		    TV = VADD(TN, TU);
		    T11 = VSUB(TN, TU);
		    TZ = VMUL(LDK(KP707106781), VSUB(TT, TQ));
		    T10 = VBYI(VADD(TY, TZ));
		    T12 = VBYI(VSUB(TZ, TY));
		    ST(&(xo[WS(os, 14)]), VSUB(TV, T10), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 6)]), VADD(T11, T12), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(TV, T10), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 10)]), VSUB(T11, T12), ovs, &(xo[0]));
	       }
	       {
		    Var Tr, TB, TA, TC;
		    {
			 Var Tf, Tq, Tw, Tz;
			 Tf = VSUB(T7, Te);
			 Tq = VSUB(Tm, Tp);
			 Tr = VBYI(VSUB(Tf, Tq));
			 TB = VBYI(VADD(Tq, Tf));
			 Tw = VADD(Tu, Tv);
			 Tz = VADD(Tx, Ty);
			 TA = VSUB(Tw, Tz);
			 TC = VADD(Tw, Tz);
		    }
		    ST(&(xo[WS(os, 7)]), VADD(Tr, TA), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 15)]), VSUB(TC, TB), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 9)]), VSUB(TA, Tr), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 1)]), VADD(TB, TC), ovs, &(xo[WS(os, 1)]));
	       }
	       {
		    Var TF, TJ, TI, TK;
		    {
			 Var TD, TE, TG, TH;
			 TD = VSUB(Tu, Tv);
			 TE = VADD(Te, T7);
			 TF = VADD(TD, TE);
			 TJ = VSUB(TD, TE);
			 TG = VADD(Tp, Tm);
			 TH = VSUB(Ty, Tx);
			 TI = VBYI(VADD(TG, TH));
			 TK = VBYI(VSUB(TH, TG));
		    }
		    ST(&(xo[WS(os, 13)]), VSUB(TF, TI), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 5)]), VADD(TJ, TK), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 3)]), VADD(TF, TI), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 11)]), VSUB(TJ, TK), ovs, &(xo[WS(os, 1)]));
	       }
	  }
}
static void n1fv_20(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(40, is), MAKE_VOLATILE_STRIDE(40, os)) {
	       Var T3, T1B, Tm, T1i, TG, TN, TO, TH, T13, T16, T1k, T1u, T1v, T1z, T1r;
	       Var T1s, T1y, T1a, T1d, T1j, Ti, TD, TB, TL, Tj, TC;
	       {
		    Var T1, T2, T1g, Tk, Tl, T1h;
		    T1 = LD(&(xi[0]), ivs, &(xi[0]));
		    T2 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
		    T1g = VADD(T1, T2);
		    Tk = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
		    Tl = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
		    T1h = VADD(Tk, Tl);
		    T3 = VSUB(T1, T2);
		    T1B = VADD(T1g, T1h);
		    Tm = VSUB(Tk, Tl);
		    T1i = VSUB(T1g, T1h);
	       }
	       {
		    Var T6, T18, Tw, T12, Tz, T15, T9, T1b, Td, T11, Tp, T19, Ts, T1c, Tg;
		    Var T14;
		    {
			 Var T4, T5, Tu, Tv;
			 T4 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 T6 = VSUB(T4, T5);
			 T18 = VADD(T4, T5);
			 Tu = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 Tv = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 Tw = VSUB(Tu, Tv);
			 T12 = VADD(Tu, Tv);
		    }
		    {
			 Var Tx, Ty, T7, T8;
			 Tx = LD(&(xi[WS(is, 17)]), ivs, &(xi[WS(is, 1)]));
			 Ty = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 Tz = VSUB(Tx, Ty);
			 T15 = VADD(Tx, Ty);
			 T7 = LD(&(xi[WS(is, 16)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 T9 = VSUB(T7, T8);
			 T1b = VADD(T7, T8);
		    }
		    {
			 Var Tb, Tc, Tn, To;
			 Tb = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 Tc = LD(&(xi[WS(is, 18)]), ivs, &(xi[0]));
			 Td = VSUB(Tb, Tc);
			 T11 = VADD(Tb, Tc);
			 Tn = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 To = LD(&(xi[WS(is, 19)]), ivs, &(xi[WS(is, 1)]));
			 Tp = VSUB(Tn, To);
			 T19 = VADD(Tn, To);
		    }
		    {
			 Var Tq, Tr, Te, Tf;
			 Tq = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tr = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Ts = VSUB(Tq, Tr);
			 T1c = VADD(Tq, Tr);
			 Te = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 Tf = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 Tg = VSUB(Te, Tf);
			 T14 = VADD(Te, Tf);
		    }
		    TG = VSUB(Ts, Tp);
		    TN = VSUB(T6, T9);
		    TO = VSUB(Td, Tg);
		    TH = VSUB(Tz, Tw);
		    T13 = VSUB(T11, T12);
		    T16 = VSUB(T14, T15);
		    T1k = VADD(T13, T16);
		    T1u = VADD(T11, T12);
		    T1v = VADD(T14, T15);
		    T1z = VADD(T1u, T1v);
		    T1r = VADD(T18, T19);
		    T1s = VADD(T1b, T1c);
		    T1y = VADD(T1r, T1s);
		    T1a = VSUB(T18, T19);
		    T1d = VSUB(T1b, T1c);
		    T1j = VADD(T1a, T1d);
		    {
			 Var Ta, Th, Tt, TA;
			 Ta = VADD(T6, T9);
			 Th = VADD(Td, Tg);
			 Ti = VADD(Ta, Th);
			 TD = VMUL(LDK(KP559016994), VSUB(Ta, Th));
			 Tt = VADD(Tp, Ts);
			 TA = VADD(Tw, Tz);
			 TB = VADD(Tt, TA);
			 TL = VMUL(LDK(KP559016994), VSUB(TA, Tt));
		    }
	       }
	       Tj = VADD(T3, Ti);
	       TC = VBYI(VADD(Tm, TB));
	       ST(&(xo[WS(os, 5)]), VSUB(Tj, TC), ovs, &(xo[WS(os, 1)]));
	       ST(&(xo[WS(os, 15)]), VADD(Tj, TC), ovs, &(xo[WS(os, 1)]));
	       {
		    Var T1A, T1C, T1D, T1x, T1G, T1t, T1w, T1F, T1E;
		    T1A = VMUL(LDK(KP559016994), VSUB(T1y, T1z));
		    T1C = VADD(T1y, T1z);
		    T1D = VFNMS(LDK(KP250000000), T1C, T1B);
		    T1t = VSUB(T1r, T1s);
		    T1w = VSUB(T1u, T1v);
		    T1x = VBYI(VFMA(LDK(KP951056516), T1t, VMUL(LDK(KP587785252), T1w)));
		    T1G = VBYI(VFNMS(LDK(KP587785252), T1t, VMUL(LDK(KP951056516), T1w)));
		    ST(&(xo[0]), VADD(T1B, T1C), ovs, &(xo[0]));
		    T1F = VSUB(T1D, T1A);
		    ST(&(xo[WS(os, 8)]), VSUB(T1F, T1G), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 12)]), VADD(T1G, T1F), ovs, &(xo[0]));
		    T1E = VADD(T1A, T1D);
		    ST(&(xo[WS(os, 4)]), VADD(T1x, T1E), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 16)]), VSUB(T1E, T1x), ovs, &(xo[0]));
	       }
	       {
		    Var T1n, T1l, T1m, T1f, T1q, T17, T1e, T1p, T1o;
		    T1n = VMUL(LDK(KP559016994), VSUB(T1j, T1k));
		    T1l = VADD(T1j, T1k);
		    T1m = VFNMS(LDK(KP250000000), T1l, T1i);
		    T17 = VSUB(T13, T16);
		    T1e = VSUB(T1a, T1d);
		    T1f = VBYI(VFNMS(LDK(KP587785252), T1e, VMUL(LDK(KP951056516), T17)));
		    T1q = VBYI(VFMA(LDK(KP951056516), T1e, VMUL(LDK(KP587785252), T17)));
		    ST(&(xo[WS(os, 10)]), VADD(T1i, T1l), ovs, &(xo[0]));
		    T1p = VADD(T1n, T1m);
		    ST(&(xo[WS(os, 6)]), VSUB(T1p, T1q), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 14)]), VADD(T1q, T1p), ovs, &(xo[0]));
		    T1o = VSUB(T1m, T1n);
		    ST(&(xo[WS(os, 2)]), VADD(T1f, T1o), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 18)]), VSUB(T1o, T1f), ovs, &(xo[0]));
	       }
	       {
		    Var TI, TP, TX, TU, TM, TW, TF, TT, TK, TE;
		    TI = VFMA(LDK(KP951056516), TG, VMUL(LDK(KP587785252), TH));
		    TP = VFMA(LDK(KP951056516), TN, VMUL(LDK(KP587785252), TO));
		    TX = VFNMS(LDK(KP587785252), TN, VMUL(LDK(KP951056516), TO));
		    TU = VFNMS(LDK(KP587785252), TG, VMUL(LDK(KP951056516), TH));
		    TK = VFMS(LDK(KP250000000), TB, Tm);
		    TM = VADD(TK, TL);
		    TW = VSUB(TL, TK);
		    TE = VFNMS(LDK(KP250000000), Ti, T3);
		    TF = VADD(TD, TE);
		    TT = VSUB(TE, TD);
		    {
			 Var TJ, TQ, TZ, T10;
			 TJ = VADD(TF, TI);
			 TQ = VBYI(VSUB(TM, TP));
			 ST(&(xo[WS(os, 19)]), VSUB(TJ, TQ), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 1)]), VADD(TJ, TQ), ovs, &(xo[WS(os, 1)]));
			 TZ = VADD(TT, TU);
			 T10 = VBYI(VADD(TX, TW));
			 ST(&(xo[WS(os, 13)]), VSUB(TZ, T10), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 7)]), VADD(TZ, T10), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var TR, TS, TV, TY;
			 TR = VSUB(TF, TI);
			 TS = VBYI(VADD(TP, TM));
			 ST(&(xo[WS(os, 11)]), VSUB(TR, TS), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 9)]), VADD(TR, TS), ovs, &(xo[WS(os, 1)]));
			 TV = VSUB(TT, TU);
			 TY = VBYI(VSUB(TW, TX));
			 ST(&(xo[WS(os, 17)]), VSUB(TV, TY), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 3)]), VADD(TV, TY), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
static void n1fv_25(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
	// ..... ..... ..... .....
	// 01234 abcde klmno uvwx
	// 56789 fghij pqrst yzAB
	
	// 0t234 uA6ds khlBm qjrf
	// 5a789 ybci1 CnDop wvex

     DVK(KP809016994, +0.809016994374947424102293417182819058860154590);  //0
     DVK(KP309016994, +0.309016994374947424102293417182819058860154590);  //t
     DVK(KP992114701, +0.992114701314477831049793042785778521453036709);  //2
     DVK(KP968583161, +0.968583161128631119490168375464735813836012403);  //3
     DVK(KP876306680, +0.876306680043863587308115903922062583399064238);  //4
     DVK(KP587785252, +0.587785252292473129168705954639072768597652438);  //5
     DVK(KP951056516, +0.951056516295153572116439333379382143405698634);  //a
     DVK(KP125333233, +0.125333233564304245373118759816508793942918247);  //7
     DVK(KP248689887, +0.248689887164854788242283746006447968417567406);  //8
     DVK(KP481753674, +0.481753674101715274987191502872129653528542010);  //9
     DVK(KP728968627, +0.728968627421411523146730319055259111372571664);  //u
     DVK(KP637423989, +0.637423989748689710176712811676016195434917298);  //A
     DVK(KP535826794, +0.535826794978996618271308767867639978063575346);  //6
     DVK(KP425779291, +0.425779291565072648862502445744251703979973042);  //d
     DVK(KP062790519, +0.062790519529313376076178224565631133122484832);  //s
     DVK(KP684547105, +0.684547105928688673732283357621209269889519233);  //y
     DVK(KP770513242, +0.770513242775789230803009636396177847271667672);  //b
     DVK(KP844327925, +0.844327925502015078548558063966681505381659241);  //c
     DVK(KP904827052, +0.904827052466019527713668647932697593970413911);  //i
     DVK(KP998026728, +0.998026728428271561952336806863450553336905220);  //1
     DVK(KP1_984229402, +1.984229402628955662099586085571557042906073418);//k
     DVK(KP1_937166322, +1.937166322257262238980336750929471627672024806);//h
     DVK(KP1_752613360, +1.752613360087727174616231807844125166798128477);//l
     DVK(KP1_457937254, +1.457937254842823046293460638110518222745143328);//B
     DVK(KP1_274847979, +1.274847979497379420353425623352032390869834596);//m
     DVK(KP250666467, +0.250666467128608490746237519633017587885836494);  //C
     DVK(KP497379774, +0.497379774329709576484567492012895936835134813);  //n
     DVK(KP963507348, +0.963507348203430549974383005744259307057084020);  //D
     DVK(KP1_369094211, +1.369094211857377347464566715242418539779038465);//o
     DVK(KP1_541026485, +1.541026485551578461606019272792355694543335344);//p
     DVK(KP1_071653589, +1.071653589957993236542617535735279956127150691);//q
     DVK(KP851558583, +0.851558583130145297725004891488503407959946084);  //j
     DVK(KP125581039, +0.125581039058626752152356449131262266244969664);  //r
     DVK(KP293892626, +0.293892626146236564584352977319536384298826219);  //f
     DVK(KP1_688655851, +1.688655851004030157097116127933363010763318483);//w
     DVK(KP1_809654104, +1.809654104932039055427337295865395187940827822);//v
     DVK(KP1_996053456, +1.996053456856543123904673613726901106673810439);//e
     DVK(KP475528258, +0.475528258147576786058219666689691071702849317);  //x
     DVK(KP250000000, +0.250000000000000000000000000000000000000000000);  //g
     DVK(KP559016994, +0.559016994374947424102293417182819058860154590);  //z

	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(50, is), MAKE_VOLATILE_STRIDE(50, os)) {
	       Var T7, T1g, T26, Ta, T2R, T2N, T2O, T2P, T19, T1Y, T16, T1Z, T1a, T2v, T1l;
	       Var T2m, TU, T21, TR, T22, TV, T2u, T1k, T2l, T2K, T2L, T2M, TE, T1R, TB;
	       Var T1S, TF, T2r, T1i, T2j, Tp, T1U, Tm, T1V, Tq, T2s, T1h, T2i;
	       {
		    Var T8, T6, T1f, T3, T1e, T25, T9;
		    T8 = LD(&(xi[0]), ivs, &(xi[0]));
		    {
			 Var T4, T5, T1, T2;
			 T4 = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
			 T6 = VADD(T4, T5);
			 T1f = VSUB(T4, T5);
			 T1 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 T2 = LD(&(xi[WS(is, 20)]), ivs, &(xi[0]));
			 T3 = VADD(T1, T2);
			 T1e = VSUB(T1, T2);
		    }
		    T7 = VMUL(LDK(KP559016994), VSUB(T3, T6));
		    T1g = VFMA(LDK(KP951056516), T1e, VMUL(LDK(KP587785252), T1f));
		    T25 = VMUL(LDK(KP951056516), T1f);
		    T26 = VFNMS(LDK(KP587785252), T1e, T25);
		    T9 = VADD(T3, T6);
		    Ta = VFNMS(LDK(KP250000000), T9, T8);
		    T2R = VADD(T8, T9);
	       }
	       {
		    Var TO, T13, TN, TT, TP, TS, T12, T18, T14, T17, T15, TQ;
		    TO = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
		    T13 = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
		    {
			 Var TH, TI, TJ, TK, TL, TM;
			 TH = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 TI = LD(&(xi[WS(is, 22)]), ivs, &(xi[0]));
			 TJ = VADD(TH, TI);
			 TK = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 TL = LD(&(xi[WS(is, 17)]), ivs, &(xi[WS(is, 1)]));
			 TM = VADD(TK, TL);
			 TN = VMUL(LDK(KP559016994), VSUB(TJ, TM));
			 TT = VSUB(TK, TL);
			 TP = VADD(TJ, TM);
			 TS = VSUB(TH, TI);
		    }
		    {
			 Var TW, TX, TY, TZ, T10, T11;
			 TW = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 TX = LD(&(xi[WS(is, 23)]), ivs, &(xi[WS(is, 1)]));
			 TY = VADD(TW, TX);
			 TZ = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 T10 = LD(&(xi[WS(is, 18)]), ivs, &(xi[0]));
			 T11 = VADD(TZ, T10);
			 T12 = VMUL(LDK(KP559016994), VSUB(TY, T11));
			 T18 = VSUB(TZ, T10);
			 T14 = VADD(TY, T11);
			 T17 = VSUB(TW, TX);
		    }
		    T2N = VADD(TO, TP);
		    T2O = VADD(T13, T14);
		    T2P = VADD(T2N, T2O);
		    T19 = VFMA(LDK(KP475528258), T17, VMUL(LDK(KP293892626), T18));
		    T1Y = VFNMS(LDK(KP293892626), T17, VMUL(LDK(KP475528258), T18));
		    T15 = VFNMS(LDK(KP250000000), T14, T13);
		    T16 = VADD(T12, T15);
		    T1Z = VSUB(T15, T12);
		    T1a = VFNMS(LDK(KP1_369094211), T19, VMUL(LDK(KP728968627), T16));
		    T2v = VFMA(LDK(KP1_996053456), T1Y, VMUL(LDK(KP062790519), T1Z));
		    T1l = VFMA(LDK(KP1_457937254), T19, VMUL(LDK(KP684547105), T16));
		    T2m = VFNMS(LDK(KP998026728), T1Z, VMUL(LDK(KP125581039), T1Y));
		    TU = VFMA(LDK(KP475528258), TS, VMUL(LDK(KP293892626), TT));
		    T21 = VFNMS(LDK(KP293892626), TS, VMUL(LDK(KP475528258), TT));
		    TQ = VFNMS(LDK(KP250000000), TP, TO);
		    TR = VADD(TN, TQ);
		    T22 = VSUB(TQ, TN);
		    TV = VFNMS(LDK(KP963507348), TU, VMUL(LDK(KP876306680), TR));
		    T2u = VFMA(LDK(KP1_688655851), T21, VMUL(LDK(KP535826794), T22));
		    T1k = VFMA(LDK(KP1_752613360), TU, VMUL(LDK(KP481753674), TR));
		    T2l = VFNMS(LDK(KP844327925), T22, VMUL(LDK(KP1_071653589), T21));
	       }
	       {
		    Var Tj, Ty, Ti, To, Tk, Tn, Tx, TD, Tz, TC, TA, Tl;
		    Tj = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
		    Ty = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
		    {
			 Var Tc, Td, Te, Tf, Tg, Th;
			 Tc = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Td = LD(&(xi[WS(is, 21)]), ivs, &(xi[WS(is, 1)]));
			 Te = VADD(Tc, Td);
			 Tf = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 Tg = LD(&(xi[WS(is, 16)]), ivs, &(xi[0]));
			 Th = VADD(Tf, Tg);
			 Ti = VMUL(LDK(KP559016994), VSUB(Te, Th));
			 To = VSUB(Tf, Tg);
			 Tk = VADD(Te, Th);
			 Tn = VSUB(Tc, Td);
		    }
		    {
			 Var Tr, Ts, Tt, Tu, Tv, Tw;
			 Tr = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 Ts = LD(&(xi[WS(is, 24)]), ivs, &(xi[0]));
			 Tt = VADD(Tr, Ts);
			 Tu = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 Tv = LD(&(xi[WS(is, 19)]), ivs, &(xi[WS(is, 1)]));
			 Tw = VADD(Tu, Tv);
			 Tx = VMUL(LDK(KP559016994), VSUB(Tt, Tw));
			 TD = VSUB(Tu, Tv);
			 Tz = VADD(Tt, Tw);
			 TC = VSUB(Tr, Ts);
		    }
		    T2K = VADD(Tj, Tk);
		    T2L = VADD(Ty, Tz);
		    T2M = VADD(T2K, T2L);
		    TE = VFMA(LDK(KP475528258), TC, VMUL(LDK(KP293892626), TD));
		    T1R = VFNMS(LDK(KP293892626), TC, VMUL(LDK(KP475528258), TD));
		    TA = VFNMS(LDK(KP250000000), Tz, Ty);
		    TB = VADD(Tx, TA);
		    T1S = VSUB(TA, Tx);
		    TF = VFNMS(LDK(KP1_688655851), TE, VMUL(LDK(KP535826794), TB));
		    T2r = VFNMS(LDK(KP425779291), T1S, VMUL(LDK(KP1_809654104), T1R));
		    T1i = VFMA(LDK(KP1_071653589), TE, VMUL(LDK(KP844327925), TB));
		    T2j = VFMA(LDK(KP851558583), T1R, VMUL(LDK(KP904827052), T1S));
		    Tp = VFMA(LDK(KP475528258), Tn, VMUL(LDK(KP293892626), To));
		    T1U = VFNMS(LDK(KP293892626), Tn, VMUL(LDK(KP475528258), To));
		    Tl = VFNMS(LDK(KP250000000), Tk, Tj);
		    Tm = VADD(Ti, Tl);
		    T1V = VSUB(Tl, Ti);
		    Tq = VFNMS(LDK(KP497379774), Tp, VMUL(LDK(KP968583161), Tm));
		    T2s = VFMA(LDK(KP963507348), T1U, VMUL(LDK(KP876306680), T1V));
		    T1h = VFMA(LDK(KP1_937166322), Tp, VMUL(LDK(KP248689887), Tm));
		    T2i = VFNMS(LDK(KP481753674), T1V, VMUL(LDK(KP1_752613360), T1U));
	       }
	       {
		    Var T2Q, T2S, T2T, T2X, T2Y, T2V, T2W, T2Z, T2U;
		    T2Q = VMUL(LDK(KP559016994), VSUB(T2M, T2P));
		    T2S = VADD(T2M, T2P);
		    T2T = VFNMS(LDK(KP250000000), T2S, T2R);
		    T2V = VSUB(T2K, T2L);
		    T2W = VSUB(T2N, T2O);
		    T2X = VBYI(VFMA(LDK(KP951056516), T2V, VMUL(LDK(KP587785252), T2W)));
		    T2Y = VBYI(VFNMS(LDK(KP587785252), T2V, VMUL(LDK(KP951056516), T2W)));
		    ST(&(xo[0]), VADD(T2R, T2S), ovs, &(xo[0]));
		    T2Z = VSUB(T2T, T2Q);
		    ST(&(xo[WS(os, 10)]), VADD(T2Y, T2Z), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 15)]), VSUB(T2Z, T2Y), ovs, &(xo[WS(os, 1)]));
		    T2U = VADD(T2Q, T2T);
		    ST(&(xo[WS(os, 5)]), VSUB(T2U, T2X), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 20)]), VADD(T2X, T2U), ovs, &(xo[0]));
	       }
	       {
		    Var T2t, T2y, T2z, T2w, T1T, T1W, T1X, T2c, T2d, T2e, T29, T2a, T2b, T20, T23;
		    Var T24, T2p, T2o, T2q, T28, T2D, T2C, T2E, T2x, T2F;
		    T2t = VSUB(T2r, T2s);
		    T2y = VADD(T2i, T2j);
		    T2z = VSUB(T2l, T2m);
		    T2w = VSUB(T2u, T2v);
		    T1T = VFNMS(LDK(KP125333233), T1S, VMUL(LDK(KP1_984229402), T1R));
		    T1W = VFMA(LDK(KP1_457937254), T1U, VMUL(LDK(KP684547105), T1V));
		    T1X = VSUB(T1T, T1W);
		    T2c = VFNMS(LDK(KP1_996053456), T21, VMUL(LDK(KP062790519), T22));
		    T2d = VFMA(LDK(KP1_541026485), T1Y, VMUL(LDK(KP637423989), T1Z));
		    T2e = VSUB(T2c, T2d);
		    T29 = VFNMS(LDK(KP1_369094211), T1U, VMUL(LDK(KP728968627), T1V));
		    T2a = VFMA(LDK(KP250666467), T1R, VMUL(LDK(KP992114701), T1S));
		    T2b = VSUB(T29, T2a);
		    T20 = VFNMS(LDK(KP770513242), T1Z, VMUL(LDK(KP1_274847979), T1Y));
		    T23 = VFMA(LDK(KP125581039), T21, VMUL(LDK(KP998026728), T22));
		    T24 = VSUB(T20, T23);
		    {
			 Var T2k, T2n, T2A, T2B;
			 T2k = VSUB(T2i, T2j);
			 T2n = VADD(T2l, T2m);
			 T2p = VADD(T2k, T2n);
			 T2o = VMUL(LDK(KP559016994), VSUB(T2k, T2n));
			 T2q = VFNMS(LDK(KP250000000), T2p, T26);
			 T28 = VSUB(Ta, T7);
			 T2A = VADD(T2s, T2r);
			 T2B = VADD(T2u, T2v);
			 T2D = VADD(T2A, T2B);
			 T2C = VMUL(LDK(KP559016994), VSUB(T2A, T2B));
			 T2E = VFNMS(LDK(KP250000000), T2D, T28);
		    }
		    {
			 Var T2I, T2J, T27, T2f;
			 T2I = VBYI(VADD(T26, T2p));
			 T2J = VADD(T28, T2D);
			 ST(&(xo[WS(os, 2)]), VADD(T2I, T2J), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 23)]), VSUB(T2J, T2I), ovs, &(xo[WS(os, 1)]));
			 T27 = VBYI(VSUB(VADD(T1X, T24), T26));
			 T2f = VADD(T28, VADD(T2b, T2e));
			 ST(&(xo[WS(os, 3)]), VADD(T27, T2f), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 22)]), VSUB(T2f, T27), ovs, &(xo[0]));
		    }
		    T2x = VBYI(VADD(T2o, VADD(T2q, VFNMS(LDK(KP587785252), T2w, VMUL(LDK(KP951056516), T2t)))));
		    T2F = VFMA(LDK(KP951056516), T2y, VFMA(LDK(KP587785252), T2z, VADD(T2C, T2E)));
		    ST(&(xo[WS(os, 7)]), VADD(T2x, T2F), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 18)]), VSUB(T2F, T2x), ovs, &(xo[0]));
		    {
			 Var T2G, T2H, T2g, T2h;
			 T2G = VBYI(VADD(T2q, VSUB(VFMA(LDK(KP587785252), T2t, VMUL(LDK(KP951056516), T2w)), T2o)));
			 T2H = VFMA(LDK(KP587785252), T2y, VSUB(VFNMS(LDK(KP951056516), T2z, T2E), T2C));
			 ST(&(xo[WS(os, 12)]), VADD(T2G, T2H), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 13)]), VSUB(T2H, T2G), ovs, &(xo[WS(os, 1)]));
			 T2g = VFMA(LDK(KP309016994), T2b, VFNMS(LDK(KP809016994), T2e, VFNMS(LDK(KP587785252), VADD(T23, T20), VFNMS(LDK(KP951056516), VADD(T1W, T1T), T28))));
			 T2h = VBYI(VSUB(VFNMS(LDK(KP587785252), VADD(T2c, T2d), VFNMS(LDK(KP809016994), T24, VFNMS(LDK(KP951056516), VADD(T29, T2a), VMUL(LDK(KP309016994), T1X)))), T26));
			 ST(&(xo[WS(os, 17)]), VSUB(T2g, T2h), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 8)]), VADD(T2g, T2h), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T1p, T1u, T1w, T1q, T1B, T1C, T1D, T1L, T1M, T1N, T1I, T1J, T1K, T1E, T1F;
		    Var T1G, T1n, T1r, T1s, Tb, T1c, T1v, T1x, T1t, T1y;
		    T1p = VSUB(TF, Tq);
		    T1u = VSUB(T1i, T1h);
		    T1w = VSUB(T1l, T1k);
		    T1q = VSUB(TV, T1a);
		    T1B = VFMA(LDK(KP1_688655851), Tp, VMUL(LDK(KP535826794), Tm));
		    T1C = VFMA(LDK(KP1_541026485), TE, VMUL(LDK(KP637423989), TB));
		    T1D = VSUB(T1B, T1C);
		    T1L = VFMA(LDK(KP851558583), TU, VMUL(LDK(KP904827052), TR));
		    T1M = VFMA(LDK(KP1_984229402), T19, VMUL(LDK(KP125333233), T16));
		    T1N = VADD(T1L, T1M);
		    T1I = VFNMS(LDK(KP844327925), Tm, VMUL(LDK(KP1_071653589), Tp));
		    T1J = VFNMS(LDK(KP1_274847979), TE, VMUL(LDK(KP770513242), TB));
		    T1K = VADD(T1I, T1J);
		    T1E = VFNMS(LDK(KP425779291), TR, VMUL(LDK(KP1_809654104), TU));
		    T1F = VFNMS(LDK(KP992114701), T16, VMUL(LDK(KP250666467), T19));
		    T1G = VADD(T1E, T1F);
		    {
			 Var T1j, T1m, TG, T1b;
			 T1j = VADD(T1h, T1i);
			 T1m = VADD(T1k, T1l);
			 T1n = VADD(T1j, T1m);
			 T1r = VFMS(LDK(KP250000000), T1n, T1g);
			 T1s = VMUL(LDK(KP559016994), VSUB(T1m, T1j));
			 Tb = VADD(T7, Ta);
			 TG = VADD(Tq, TF);
			 T1b = VADD(TV, T1a);
			 T1c = VADD(TG, T1b);
			 T1v = VFNMS(LDK(KP250000000), T1c, Tb);
			 T1x = VMUL(LDK(KP559016994), VSUB(TG, T1b));
		    }
		    {
			 Var T1d, T1o, T1H, T1O;
			 T1d = VADD(Tb, T1c);
			 T1o = VBYI(VADD(T1g, T1n));
			 ST(&(xo[WS(os, 1)]), VSUB(T1d, T1o), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 24)]), VADD(T1d, T1o), ovs, &(xo[0]));
			 T1H = VADD(Tb, VADD(T1D, T1G));
			 T1O = VBYI(VADD(T1g, VSUB(T1K, T1N)));
			 ST(&(xo[WS(os, 21)]), VSUB(T1H, T1O), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 4)]), VADD(T1H, T1O), ovs, &(xo[0]));
		    }
		    T1t = VBYI(VADD(VFMA(LDK(KP587785252), T1p, VMUL(LDK(KP951056516), T1q)), VSUB(T1r, T1s)));
		    T1y = VFMA(LDK(KP587785252), T1u, VFNMS(LDK(KP951056516), T1w, VSUB(T1v, T1x)));
		    ST(&(xo[WS(os, 11)]), VADD(T1t, T1y), ovs, &(xo[WS(os, 1)]));
		    ST(&(xo[WS(os, 14)]), VSUB(T1y, T1t), ovs, &(xo[0]));
		    {
			 Var T1z, T1A, T1P, T1Q;
			 T1z = VBYI(VADD(VFNMS(LDK(KP587785252), T1q, VMUL(LDK(KP951056516), T1p)), VADD(T1r, T1s)));
			 T1A = VFMA(LDK(KP951056516), T1u, VADD(T1x, VFMA(LDK(KP587785252), T1w, T1v)));
			 ST(&(xo[WS(os, 6)]), VADD(T1z, T1A), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 19)]), VSUB(T1A, T1z), ovs, &(xo[WS(os, 1)]));
			 T1P = VBYI(VADD(T1g, VFMA(LDK(KP309016994), T1K, VFMA(LDK(KP587785252), VSUB(T1F, T1E), VFNMS(LDK(KP951056516), VADD(T1B, T1C), VMUL(LDK(KP809016994), T1N))))));
			 T1Q = VFMA(LDK(KP309016994), T1D, VFMA(LDK(KP951056516), VSUB(T1I, T1J), VFMA(LDK(KP587785252), VSUB(T1M, T1L), VFNMS(LDK(KP809016994), T1G, Tb))));
			 ST(&(xo[WS(os, 9)]), VADD(T1P, T1Q), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 16)]), VSUB(T1Q, T1P), ovs, &(xo[0]));
		    }
	       }
	  }
}
static void n1fv_32(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);//5
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);//1
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);//3
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);//6
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);//4
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);//0
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);//2

	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(64, is), MAKE_VOLATILE_STRIDE(64, os)) {
	       Var T1T, T1W, T2K, T2x, T16, T1A, Tb, T1p, TT, T1v, TY, T1w, T27, T2a, T2b;
	       Var T2H, T2O, TC, T1s, TH, T1t, T20, T23, T24, T2E, T2N, T2g, T2j, Tq, T1B;
	       Var T19, T1q, T2A, T2L;
	       {
		    Var T3, T1R, T15, T1S, T6, T1U, T9, T1V, T12, Ta;
		    {
			 Var T1, T2, T13, T14;
			 T1 = LD(&(xi[0]), ivs, &(xi[0]));
			 T2 = LD(&(xi[WS(is, 16)]), ivs, &(xi[0]));
			 T3 = VSUB(T1, T2);
			 T1R = VADD(T1, T2);
			 T13 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T14 = LD(&(xi[WS(is, 24)]), ivs, &(xi[0]));
			 T15 = VSUB(T13, T14);
			 T1S = VADD(T13, T14);
		    }
		    {
			 Var T4, T5, T7, T8;
			 T4 = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 20)]), ivs, &(xi[0]));
			 T6 = VSUB(T4, T5);
			 T1U = VADD(T4, T5);
			 T7 = LD(&(xi[WS(is, 28)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 T9 = VSUB(T7, T8);
			 T1V = VADD(T7, T8);
		    }
		    T1T = VADD(T1R, T1S);
		    T1W = VADD(T1U, T1V);
		    T2K = VSUB(T1V, T1U);
		    T2x = VSUB(T1R, T1S);
		    T12 = VMUL(LDK(KP707106781), VSUB(T9, T6));
		    T16 = VSUB(T12, T15);
		    T1A = VADD(T15, T12);
		    Ta = VMUL(LDK(KP707106781), VADD(T6, T9));
		    Tb = VADD(T3, Ta);
		    T1p = VSUB(T3, Ta);
	       }
	       {
		    Var TL, T25, TX, T26, TO, T28, TR, T29;
		    {
			 Var TJ, TK, TV, TW;
			 TJ = LD(&(xi[WS(is, 31)]), ivs, &(xi[WS(is, 1)]));
			 TK = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
			 TL = VSUB(TJ, TK);
			 T25 = VADD(TJ, TK);
			 TV = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 TW = LD(&(xi[WS(is, 23)]), ivs, &(xi[WS(is, 1)]));
			 TX = VSUB(TV, TW);
			 T26 = VADD(TV, TW);
		    }
		    {
			 Var TM, TN, TP, TQ;
			 TM = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 TN = LD(&(xi[WS(is, 19)]), ivs, &(xi[WS(is, 1)]));
			 TO = VSUB(TM, TN);
			 T28 = VADD(TM, TN);
			 TP = LD(&(xi[WS(is, 27)]), ivs, &(xi[WS(is, 1)]));
			 TQ = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			 TR = VSUB(TP, TQ);
			 T29 = VADD(TP, TQ);
		    }
		    {
			 Var TS, TU, T2F, T2G;
			 TS = VMUL(LDK(KP707106781), VADD(TO, TR));
			 TT = VADD(TL, TS);
			 T1v = VSUB(TL, TS);
			 TU = VMUL(LDK(KP707106781), VSUB(TR, TO));
			 TY = VSUB(TU, TX);
			 T1w = VADD(TX, TU);
			 T27 = VADD(T25, T26);
			 T2a = VADD(T28, T29);
			 T2b = VSUB(T27, T2a);
			 T2F = VSUB(T25, T26);
			 T2G = VSUB(T29, T28);
			 T2H = VFNMS(LDK(KP382683432), T2G, VMUL(LDK(KP923879532), T2F));
			 T2O = VFMA(LDK(KP382683432), T2F, VMUL(LDK(KP923879532), T2G));
		    }
	       }
	       {
		    Var Tu, T1Y, TG, T1Z, Tx, T21, TA, T22;
		    {
			 Var Ts, Tt, TE, TF;
			 Ts = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 Tt = LD(&(xi[WS(is, 17)]), ivs, &(xi[WS(is, 1)]));
			 Tu = VSUB(Ts, Tt);
			 T1Y = VADD(Ts, Tt);
			 TE = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 TF = LD(&(xi[WS(is, 25)]), ivs, &(xi[WS(is, 1)]));
			 TG = VSUB(TE, TF);
			 T1Z = VADD(TE, TF);
		    }
		    {
			 Var Tv, Tw, Ty, Tz;
			 Tv = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 Tw = LD(&(xi[WS(is, 21)]), ivs, &(xi[WS(is, 1)]));
			 Tx = VSUB(Tv, Tw);
			 T21 = VADD(Tv, Tw);
			 Ty = LD(&(xi[WS(is, 29)]), ivs, &(xi[WS(is, 1)]));
			 Tz = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			 TA = VSUB(Ty, Tz);
			 T22 = VADD(Ty, Tz);
		    }
		    {
			 Var TB, TD, T2C, T2D;
			 TB = VMUL(LDK(KP707106781), VADD(Tx, TA));
			 TC = VADD(Tu, TB);
			 T1s = VSUB(Tu, TB);
			 TD = VMUL(LDK(KP707106781), VSUB(TA, Tx));
			 TH = VSUB(TD, TG);
			 T1t = VADD(TG, TD);
			 T20 = VADD(T1Y, T1Z);
			 T23 = VADD(T21, T22);
			 T24 = VSUB(T20, T23);
			 T2C = VSUB(T1Y, T1Z);
			 T2D = VSUB(T22, T21);
			 T2E = VFMA(LDK(KP923879532), T2C, VMUL(LDK(KP382683432), T2D));
			 T2N = VFNMS(LDK(KP382683432), T2C, VMUL(LDK(KP923879532), T2D));
		    }
	       }
	       {
		    Var Te, T2h, To, T2f, Th, T2i, Tl, T2e, Ti, Tp;
		    {
			 Var Tc, Td, Tm, Tn;
			 Tc = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 Td = LD(&(xi[WS(is, 18)]), ivs, &(xi[0]));
			 Te = VSUB(Tc, Td);
			 T2h = VADD(Tc, Td);
			 Tm = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			 Tn = LD(&(xi[WS(is, 22)]), ivs, &(xi[0]));
			 To = VSUB(Tm, Tn);
			 T2f = VADD(Tm, Tn);
		    }
		    {
			 Var Tf, Tg, Tj, Tk;
			 Tf = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			 Tg = LD(&(xi[WS(is, 26)]), ivs, &(xi[0]));
			 Th = VSUB(Tf, Tg);
			 T2i = VADD(Tf, Tg);
			 Tj = LD(&(xi[WS(is, 30)]), ivs, &(xi[0]));
			 Tk = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 Tl = VSUB(Tj, Tk);
			 T2e = VADD(Tj, Tk);
		    }
		    T2g = VADD(T2e, T2f);
		    T2j = VADD(T2h, T2i);
		    Ti = VFNMS(LDK(KP382683432), Th, VMUL(LDK(KP923879532), Te));
		    Tp = VFMA(LDK(KP923879532), Tl, VMUL(LDK(KP382683432), To));
		    Tq = VADD(Ti, Tp);
		    T1B = VSUB(Tp, Ti);
		    {
			 Var T17, T18, T2y, T2z;
			 T17 = VFNMS(LDK(KP923879532), To, VMUL(LDK(KP382683432), Tl));
			 T18 = VFMA(LDK(KP382683432), Te, VMUL(LDK(KP923879532), Th));
			 T19 = VSUB(T17, T18);
			 T1q = VADD(T18, T17);
			 T2y = VSUB(T2h, T2i);
			 T2z = VSUB(T2e, T2f);
			 T2A = VMUL(LDK(KP707106781), VADD(T2y, T2z));
			 T2L = VMUL(LDK(KP707106781), VSUB(T2z, T2y));
		    }
	       }
	       {
		    Var T2d, T2n, T2m, T2o;
		    {
			 Var T1X, T2c, T2k, T2l;
			 T1X = VSUB(T1T, T1W);
			 T2c = VMUL(LDK(KP707106781), VADD(T24, T2b));
			 T2d = VADD(T1X, T2c);
			 T2n = VSUB(T1X, T2c);
			 T2k = VSUB(T2g, T2j);
			 T2l = VMUL(LDK(KP707106781), VSUB(T2b, T24));
			 T2m = VBYI(VADD(T2k, T2l));
			 T2o = VBYI(VSUB(T2l, T2k));
		    }
		    ST(&(xo[WS(os, 28)]), VSUB(T2d, T2m), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 12)]), VADD(T2n, T2o), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 4)]), VADD(T2d, T2m), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 20)]), VSUB(T2n, T2o), ovs, &(xo[0]));
	       }
	       {
		    Var T2r, T2v, T2u, T2w;
		    {
			 Var T2p, T2q, T2s, T2t;
			 T2p = VADD(T1T, T1W);
			 T2q = VADD(T2j, T2g);
			 T2r = VADD(T2p, T2q);
			 T2v = VSUB(T2p, T2q);
			 T2s = VADD(T20, T23);
			 T2t = VADD(T27, T2a);
			 T2u = VADD(T2s, T2t);
			 T2w = VBYI(VSUB(T2t, T2s));
		    }
		    ST(&(xo[WS(os, 16)]), VSUB(T2r, T2u), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 8)]), VADD(T2v, T2w), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(T2r, T2u), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 24)]), VSUB(T2v, T2w), ovs, &(xo[0]));
	       }
	       {
		    Var T2V, T2Z, T2Y, T30;
		    {
			 Var T2T, T2U, T2W, T2X;
			 T2T = VSUB(T2H, T2E);
			 T2U = VSUB(T2L, T2K);
			 T2V = VBYI(VSUB(T2T, T2U));
			 T2Z = VBYI(VADD(T2U, T2T));
			 T2W = VSUB(T2x, T2A);
			 T2X = VSUB(T2O, T2N);
			 T2Y = VSUB(T2W, T2X);
			 T30 = VADD(T2W, T2X);
		    }
		    ST(&(xo[WS(os, 10)]), VADD(T2V, T2Y), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 26)]), VSUB(T30, T2Z), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 22)]), VSUB(T2Y, T2V), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 6)]), VADD(T2Z, T30), ovs, &(xo[0]));
	       }
	       {
		    Var T2J, T2R, T2Q, T2S;
		    {
			 Var T2B, T2I, T2M, T2P;
			 T2B = VADD(T2x, T2A);
			 T2I = VADD(T2E, T2H);
			 T2J = VADD(T2B, T2I);
			 T2R = VSUB(T2B, T2I);
			 T2M = VADD(T2K, T2L);
			 T2P = VADD(T2N, T2O);
			 T2Q = VBYI(VADD(T2M, T2P));
			 T2S = VBYI(VSUB(T2P, T2M));
		    }
		    ST(&(xo[WS(os, 30)]), VSUB(T2J, T2Q), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 14)]), VADD(T2R, T2S), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 2)]), VADD(T2J, T2Q), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 18)]), VSUB(T2R, T2S), ovs, &(xo[0]));
	       }
	       {
		    Var T1r, T1C, T1M, T1K, T1F, T1N, T1y, T1J;
		    T1r = VADD(T1p, T1q);
		    T1C = VADD(T1A, T1B);
		    T1M = VSUB(T1p, T1q);
		    T1K = VSUB(T1B, T1A);
		    {
			 Var T1D, T1E, T1u, T1x;
			 T1D = VFNMS(LDK(KP555570233), T1s, VMUL(LDK(KP831469612), T1t));
			 T1E = VFMA(LDK(KP555570233), T1v, VMUL(LDK(KP831469612), T1w));
			 T1F = VADD(T1D, T1E);
			 T1N = VSUB(T1E, T1D);
			 T1u = VFMA(LDK(KP831469612), T1s, VMUL(LDK(KP555570233), T1t));
			 T1x = VFNMS(LDK(KP555570233), T1w, VMUL(LDK(KP831469612), T1v));
			 T1y = VADD(T1u, T1x);
			 T1J = VSUB(T1x, T1u);
		    }
		    {
			 Var T1z, T1G, T1P, T1Q;
			 T1z = VADD(T1r, T1y);
			 T1G = VBYI(VADD(T1C, T1F));
			 ST(&(xo[WS(os, 29)]), VSUB(T1z, T1G), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 3)]), VADD(T1z, T1G), ovs, &(xo[WS(os, 1)]));
			 T1P = VBYI(VADD(T1K, T1J));
			 T1Q = VADD(T1M, T1N);
			 ST(&(xo[WS(os, 5)]), VADD(T1P, T1Q), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 27)]), VSUB(T1Q, T1P), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T1H, T1I, T1L, T1O;
			 T1H = VSUB(T1r, T1y);
			 T1I = VBYI(VSUB(T1F, T1C));
			 ST(&(xo[WS(os, 19)]), VSUB(T1H, T1I), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 13)]), VADD(T1H, T1I), ovs, &(xo[WS(os, 1)]));
			 T1L = VBYI(VSUB(T1J, T1K));
			 T1O = VSUB(T1M, T1N);
			 ST(&(xo[WS(os, 11)]), VADD(T1L, T1O), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 21)]), VSUB(T1O, T1L), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var Tr, T1a, T1k, T1i, T1d, T1l, T10, T1h;
		    Tr = VADD(Tb, Tq);
		    T1a = VADD(T16, T19);
		    T1k = VSUB(Tb, Tq);
		    T1i = VSUB(T19, T16);
		    {
			 Var T1b, T1c, TI, TZ;
			 T1b = VFNMS(LDK(KP195090322), TC, VMUL(LDK(KP980785280), TH));
			 T1c = VFMA(LDK(KP195090322), TT, VMUL(LDK(KP980785280), TY));
			 T1d = VADD(T1b, T1c);
			 T1l = VSUB(T1c, T1b);
			 TI = VFMA(LDK(KP980785280), TC, VMUL(LDK(KP195090322), TH));
			 TZ = VFNMS(LDK(KP195090322), TY, VMUL(LDK(KP980785280), TT));
			 T10 = VADD(TI, TZ);
			 T1h = VSUB(TZ, TI);
		    }
		    {
			 Var T11, T1e, T1n, T1o;
			 T11 = VADD(Tr, T10);
			 T1e = VBYI(VADD(T1a, T1d));
			 ST(&(xo[WS(os, 31)]), VSUB(T11, T1e), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 1)]), VADD(T11, T1e), ovs, &(xo[WS(os, 1)]));
			 T1n = VBYI(VADD(T1i, T1h));
			 T1o = VADD(T1k, T1l);
			 ST(&(xo[WS(os, 7)]), VADD(T1n, T1o), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 25)]), VSUB(T1o, T1n), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T1f, T1g, T1j, T1m;
			 T1f = VSUB(Tr, T10);
			 T1g = VBYI(VSUB(T1d, T1a));
			 ST(&(xo[WS(os, 17)]), VSUB(T1f, T1g), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 15)]), VADD(T1f, T1g), ovs, &(xo[WS(os, 1)]));
			 T1j = VBYI(VSUB(T1h, T1i));
			 T1m = VSUB(T1k, T1l);
			 ST(&(xo[WS(os, 9)]), VADD(T1j, T1m), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 23)]), VSUB(T1m, T1j), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
static void n1fv_64(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP471396736, +0.471396736825997648556387625905254377657460319);
     DVK(KP881921264, +0.881921264348355029712756863660388349508442621);
     DVK(KP098017140, +0.098017140329560601994195563888641845861136673);
     DVK(KP995184726, +0.995184726672196886244836953109479921575474869);
     DVK(KP290284677, +0.290284677254462367636192375817395274691476278);
     DVK(KP956940335, +0.956940335732208864935797886980269969482849206);
     DVK(KP634393284, +0.634393284163645498215171613225493370675687095);
     DVK(KP773010453, +0.773010453362736960810906609758469800971041293);
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(128, is), MAKE_VOLATILE_STRIDE(128, os)) {
	       Var T4p, T5q, Tb, T39, T2n, T3A, T6f, T6T, Tq, T3B, T6i, T76, T2i, T3a, T4w;
	       Var T5r, TI, T2p, T6C, T6V, T3h, T3E, T4L, T5u, TZ, T2q, T6F, T6U, T3e, T3D;
	       Var T4E, T5t, T23, T2N, T6t, T71, T6w, T72, T2c, T2O, T3t, T41, T5f, T5R, T5k;
	       Var T5S, T3w, T42, T1s, T2K, T6m, T6Y, T6p, T6Z, T1B, T2L, T3m, T3Y, T4Y, T5O;
	       Var T53, T5P, T3p, T3Z;
	       {
		    Var T3, T4n, T2m, T4o, T6, T5p, T9, T5o;
		    {
			 Var T1, T2, T2k, T2l;
			 T1 = LD(&(xi[0]), ivs, &(xi[0]));
			 T2 = LD(&(xi[WS(is, 32)]), ivs, &(xi[0]));
			 T3 = VSUB(T1, T2);
			 T4n = VADD(T1, T2);
			 T2k = LD(&(xi[WS(is, 16)]), ivs, &(xi[0]));
			 T2l = LD(&(xi[WS(is, 48)]), ivs, &(xi[0]));
			 T2m = VSUB(T2k, T2l);
			 T4o = VADD(T2k, T2l);
		    }
		    {
			 Var T4, T5, T7, T8;
			 T4 = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 40)]), ivs, &(xi[0]));
			 T6 = VSUB(T4, T5);
			 T5p = VADD(T4, T5);
			 T7 = LD(&(xi[WS(is, 56)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 24)]), ivs, &(xi[0]));
			 T9 = VSUB(T7, T8);
			 T5o = VADD(T7, T8);
		    }
		    T4p = VSUB(T4n, T4o);
		    T5q = VSUB(T5o, T5p);
		    {
			 Var Ta, T2j, T6d, T6e;
			 Ta = VMUL(LDK(KP707106781), VADD(T6, T9));
			 Tb = VADD(T3, Ta);
			 T39 = VSUB(T3, Ta);
			 T2j = VMUL(LDK(KP707106781), VSUB(T9, T6));
			 T2n = VSUB(T2j, T2m);
			 T3A = VADD(T2m, T2j);
			 T6d = VADD(T4n, T4o);
			 T6e = VADD(T5p, T5o);
			 T6f = VADD(T6d, T6e);
			 T6T = VSUB(T6d, T6e);
		    }
	       }
	       {
		    Var Te, T4q, To, T4u, Th, T4r, Tl, T4t;
		    {
			 Var Tc, Td, Tm, Tn;
			 Tc = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 Td = LD(&(xi[WS(is, 36)]), ivs, &(xi[0]));
			 Te = VSUB(Tc, Td);
			 T4q = VADD(Tc, Td);
			 Tm = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			 Tn = LD(&(xi[WS(is, 44)]), ivs, &(xi[0]));
			 To = VSUB(Tm, Tn);
			 T4u = VADD(Tm, Tn);
		    }
		    {
			 Var Tf, Tg, Tj, Tk;
			 Tf = LD(&(xi[WS(is, 20)]), ivs, &(xi[0]));
			 Tg = LD(&(xi[WS(is, 52)]), ivs, &(xi[0]));
			 Th = VSUB(Tf, Tg);
			 T4r = VADD(Tf, Tg);
			 Tj = LD(&(xi[WS(is, 60)]), ivs, &(xi[0]));
			 Tk = LD(&(xi[WS(is, 28)]), ivs, &(xi[0]));
			 Tl = VSUB(Tj, Tk);
			 T4t = VADD(Tj, Tk);
		    }
		    {
			 Var Ti, Tp, T6g, T6h;
			 Ti = VFNMS(LDK(KP382683432), Th, VMUL(LDK(KP923879532), Te));
			 Tp = VFMA(LDK(KP923879532), Tl, VMUL(LDK(KP382683432), To));
			 Tq = VADD(Ti, Tp);
			 T3B = VSUB(Tp, Ti);
			 T6g = VADD(T4q, T4r);
			 T6h = VADD(T4t, T4u);
			 T6i = VADD(T6g, T6h);
			 T76 = VSUB(T6h, T6g);
		    }
		    {
			 Var T2g, T2h, T4s, T4v;
			 T2g = VFNMS(LDK(KP923879532), To, VMUL(LDK(KP382683432), Tl));
			 T2h = VFMA(LDK(KP382683432), Te, VMUL(LDK(KP923879532), Th));
			 T2i = VSUB(T2g, T2h);
			 T3a = VADD(T2h, T2g);
			 T4s = VSUB(T4q, T4r);
			 T4v = VSUB(T4t, T4u);
			 T4w = VMUL(LDK(KP707106781), VADD(T4s, T4v));
			 T5r = VMUL(LDK(KP707106781), VSUB(T4v, T4s));
		    }
	       }
	       {
		    Var Tu, T4F, TG, T4G, TB, T4J, TD, T4I;
		    {
			 Var Ts, Tt, TE, TF;
			 Ts = LD(&(xi[WS(is, 62)]), ivs, &(xi[0]));
			 Tt = LD(&(xi[WS(is, 30)]), ivs, &(xi[0]));
			 Tu = VSUB(Ts, Tt);
			 T4F = VADD(Ts, Tt);
			 TE = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 TF = LD(&(xi[WS(is, 46)]), ivs, &(xi[0]));
			 TG = VSUB(TE, TF);
			 T4G = VADD(TE, TF);
			 {
			      Var Tv, Tw, Tx, Ty, Tz, TA;
			      Tv = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			      Tw = LD(&(xi[WS(is, 38)]), ivs, &(xi[0]));
			      Tx = VSUB(Tv, Tw);
			      Ty = LD(&(xi[WS(is, 54)]), ivs, &(xi[0]));
			      Tz = LD(&(xi[WS(is, 22)]), ivs, &(xi[0]));
			      TA = VSUB(Ty, Tz);
			      TB = VMUL(LDK(KP707106781), VADD(Tx, TA));
			      T4J = VADD(Tv, Tw);
			      TD = VMUL(LDK(KP707106781), VSUB(TA, Tx));
			      T4I = VADD(Ty, Tz);
			 }
		    }
		    {
			 Var TC, TH, T6A, T6B;
			 TC = VADD(Tu, TB);
			 TH = VSUB(TD, TG);
			 TI = VFMA(LDK(KP195090322), TC, VMUL(LDK(KP980785280), TH));
			 T2p = VFNMS(LDK(KP195090322), TH, VMUL(LDK(KP980785280), TC));
			 T6A = VADD(T4F, T4G);
			 T6B = VADD(T4J, T4I);
			 T6C = VADD(T6A, T6B);
			 T6V = VSUB(T6A, T6B);
		    }
		    {
			 Var T3f, T3g, T4H, T4K;
			 T3f = VSUB(Tu, TB);
			 T3g = VADD(TG, TD);
			 T3h = VFNMS(LDK(KP555570233), T3g, VMUL(LDK(KP831469612), T3f));
			 T3E = VFMA(LDK(KP555570233), T3f, VMUL(LDK(KP831469612), T3g));
			 T4H = VSUB(T4F, T4G);
			 T4K = VSUB(T4I, T4J);
			 T4L = VFNMS(LDK(KP382683432), T4K, VMUL(LDK(KP923879532), T4H));
			 T5u = VFMA(LDK(KP382683432), T4H, VMUL(LDK(KP923879532), T4K));
		    }
	       }
	       {
		    Var TS, T4z, TW, T4y, TP, T4C, TX, T4B;
		    {
			 Var TQ, TR, TU, TV;
			 TQ = LD(&(xi[WS(is, 18)]), ivs, &(xi[0]));
			 TR = LD(&(xi[WS(is, 50)]), ivs, &(xi[0]));
			 TS = VSUB(TQ, TR);
			 T4z = VADD(TQ, TR);
			 TU = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 TV = LD(&(xi[WS(is, 34)]), ivs, &(xi[0]));
			 TW = VSUB(TU, TV);
			 T4y = VADD(TU, TV);
			 {
			      Var TJ, TK, TL, TM, TN, TO;
			      TJ = LD(&(xi[WS(is, 58)]), ivs, &(xi[0]));
			      TK = LD(&(xi[WS(is, 26)]), ivs, &(xi[0]));
			      TL = VSUB(TJ, TK);
			      TM = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			      TN = LD(&(xi[WS(is, 42)]), ivs, &(xi[0]));
			      TO = VSUB(TM, TN);
			      TP = VMUL(LDK(KP707106781), VSUB(TL, TO));
			      T4C = VADD(TM, TN);
			      TX = VMUL(LDK(KP707106781), VADD(TO, TL));
			      T4B = VADD(TJ, TK);
			 }
		    }
		    {
			 Var TT, TY, T6D, T6E;
			 TT = VSUB(TP, TS);
			 TY = VADD(TW, TX);
			 TZ = VFNMS(LDK(KP195090322), TY, VMUL(LDK(KP980785280), TT));
			 T2q = VFMA(LDK(KP980785280), TY, VMUL(LDK(KP195090322), TT));
			 T6D = VADD(T4y, T4z);
			 T6E = VADD(T4C, T4B);
			 T6F = VADD(T6D, T6E);
			 T6U = VSUB(T6D, T6E);
		    }
		    {
			 Var T3c, T3d, T4A, T4D;
			 T3c = VSUB(TW, TX);
			 T3d = VADD(TS, TP);
			 T3e = VFMA(LDK(KP831469612), T3c, VMUL(LDK(KP555570233), T3d));
			 T3D = VFNMS(LDK(KP555570233), T3c, VMUL(LDK(KP831469612), T3d));
			 T4A = VSUB(T4y, T4z);
			 T4D = VSUB(T4B, T4C);
			 T4E = VFMA(LDK(KP923879532), T4A, VMUL(LDK(KP382683432), T4D));
			 T5t = VFNMS(LDK(KP382683432), T4A, VMUL(LDK(KP923879532), T4D));
		    }
	       }
	       {
		    Var T1F, T55, T2a, T56, T1M, T5h, T27, T5g, T58, T59, T1U, T5a, T25, T5b, T5c;
		    Var T21, T5d, T24;
		    {
			 Var T1D, T1E, T28, T29;
			 T1D = LD(&(xi[WS(is, 63)]), ivs, &(xi[WS(is, 1)]));
			 T1E = LD(&(xi[WS(is, 31)]), ivs, &(xi[WS(is, 1)]));
			 T1F = VSUB(T1D, T1E);
			 T55 = VADD(T1D, T1E);
			 T28 = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
			 T29 = LD(&(xi[WS(is, 47)]), ivs, &(xi[WS(is, 1)]));
			 T2a = VSUB(T28, T29);
			 T56 = VADD(T28, T29);
		    }
		    {
			 Var T1G, T1H, T1I, T1J, T1K, T1L;
			 T1G = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			 T1H = LD(&(xi[WS(is, 39)]), ivs, &(xi[WS(is, 1)]));
			 T1I = VSUB(T1G, T1H);
			 T1J = LD(&(xi[WS(is, 55)]), ivs, &(xi[WS(is, 1)]));
			 T1K = LD(&(xi[WS(is, 23)]), ivs, &(xi[WS(is, 1)]));
			 T1L = VSUB(T1J, T1K);
			 T1M = VMUL(LDK(KP707106781), VADD(T1I, T1L));
			 T5h = VADD(T1G, T1H);
			 T27 = VMUL(LDK(KP707106781), VSUB(T1L, T1I));
			 T5g = VADD(T1J, T1K);
		    }
		    {
			 Var T1Q, T1T, T1X, T20;
			 {
			      Var T1O, T1P, T1R, T1S;
			      T1O = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			      T1P = LD(&(xi[WS(is, 35)]), ivs, &(xi[WS(is, 1)]));
			      T1Q = VSUB(T1O, T1P);
			      T58 = VADD(T1O, T1P);
			      T1R = LD(&(xi[WS(is, 19)]), ivs, &(xi[WS(is, 1)]));
			      T1S = LD(&(xi[WS(is, 51)]), ivs, &(xi[WS(is, 1)]));
			      T1T = VSUB(T1R, T1S);
			      T59 = VADD(T1R, T1S);
			 }
			 T1U = VFNMS(LDK(KP382683432), T1T, VMUL(LDK(KP923879532), T1Q));
			 T5a = VSUB(T58, T59);
			 T25 = VFMA(LDK(KP382683432), T1Q, VMUL(LDK(KP923879532), T1T));
			 {
			      Var T1V, T1W, T1Y, T1Z;
			      T1V = LD(&(xi[WS(is, 59)]), ivs, &(xi[WS(is, 1)]));
			      T1W = LD(&(xi[WS(is, 27)]), ivs, &(xi[WS(is, 1)]));
			      T1X = VSUB(T1V, T1W);
			      T5b = VADD(T1V, T1W);
			      T1Y = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			      T1Z = LD(&(xi[WS(is, 43)]), ivs, &(xi[WS(is, 1)]));
			      T20 = VSUB(T1Y, T1Z);
			      T5c = VADD(T1Y, T1Z);
			 }
			 T21 = VFMA(LDK(KP923879532), T1X, VMUL(LDK(KP382683432), T20));
			 T5d = VSUB(T5b, T5c);
			 T24 = VFNMS(LDK(KP923879532), T20, VMUL(LDK(KP382683432), T1X));
		    }
		    {
			 Var T1N, T22, T6r, T6s;
			 T1N = VADD(T1F, T1M);
			 T22 = VADD(T1U, T21);
			 T23 = VSUB(T1N, T22);
			 T2N = VADD(T1N, T22);
			 T6r = VADD(T55, T56);
			 T6s = VADD(T5h, T5g);
			 T6t = VADD(T6r, T6s);
			 T71 = VSUB(T6r, T6s);
		    }
		    {
			 Var T6u, T6v, T26, T2b;
			 T6u = VADD(T58, T59);
			 T6v = VADD(T5b, T5c);
			 T6w = VADD(T6u, T6v);
			 T72 = VSUB(T6v, T6u);
			 T26 = VSUB(T24, T25);
			 T2b = VSUB(T27, T2a);
			 T2c = VSUB(T26, T2b);
			 T2O = VADD(T2b, T26);
		    }
		    {
			 Var T3r, T3s, T57, T5e;
			 T3r = VSUB(T1F, T1M);
			 T3s = VADD(T25, T24);
			 T3t = VADD(T3r, T3s);
			 T41 = VSUB(T3r, T3s);
			 T57 = VSUB(T55, T56);
			 T5e = VMUL(LDK(KP707106781), VADD(T5a, T5d));
			 T5f = VADD(T57, T5e);
			 T5R = VSUB(T57, T5e);
		    }
		    {
			 Var T5i, T5j, T3u, T3v;
			 T5i = VSUB(T5g, T5h);
			 T5j = VMUL(LDK(KP707106781), VSUB(T5d, T5a));
			 T5k = VADD(T5i, T5j);
			 T5S = VSUB(T5j, T5i);
			 T3u = VADD(T2a, T27);
			 T3v = VSUB(T21, T1U);
			 T3w = VADD(T3u, T3v);
			 T42 = VSUB(T3v, T3u);
		    }
	       }
	       {
		    Var T1q, T4P, T1v, T4O, T1n, T50, T1w, T4Z, T4U, T4V, T18, T4W, T1z, T4R, T4S;
		    Var T1f, T4T, T1y;
		    {
			 Var T1o, T1p, T1t, T1u;
			 T1o = LD(&(xi[WS(is, 17)]), ivs, &(xi[WS(is, 1)]));
			 T1p = LD(&(xi[WS(is, 49)]), ivs, &(xi[WS(is, 1)]));
			 T1q = VSUB(T1o, T1p);
			 T4P = VADD(T1o, T1p);
			 T1t = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 T1u = LD(&(xi[WS(is, 33)]), ivs, &(xi[WS(is, 1)]));
			 T1v = VSUB(T1t, T1u);
			 T4O = VADD(T1t, T1u);
		    }
		    {
			 Var T1h, T1i, T1j, T1k, T1l, T1m;
			 T1h = LD(&(xi[WS(is, 57)]), ivs, &(xi[WS(is, 1)]));
			 T1i = LD(&(xi[WS(is, 25)]), ivs, &(xi[WS(is, 1)]));
			 T1j = VSUB(T1h, T1i);
			 T1k = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			 T1l = LD(&(xi[WS(is, 41)]), ivs, &(xi[WS(is, 1)]));
			 T1m = VSUB(T1k, T1l);
			 T1n = VMUL(LDK(KP707106781), VSUB(T1j, T1m));
			 T50 = VADD(T1k, T1l);
			 T1w = VMUL(LDK(KP707106781), VADD(T1m, T1j));
			 T4Z = VADD(T1h, T1i);
		    }
		    {
			 Var T14, T17, T1b, T1e;
			 {
			      Var T12, T13, T15, T16;
			      T12 = LD(&(xi[WS(is, 61)]), ivs, &(xi[WS(is, 1)]));
			      T13 = LD(&(xi[WS(is, 29)]), ivs, &(xi[WS(is, 1)]));
			      T14 = VSUB(T12, T13);
			      T4U = VADD(T12, T13);
			      T15 = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			      T16 = LD(&(xi[WS(is, 45)]), ivs, &(xi[WS(is, 1)]));
			      T17 = VSUB(T15, T16);
			      T4V = VADD(T15, T16);
			 }
			 T18 = VFNMS(LDK(KP923879532), T17, VMUL(LDK(KP382683432), T14));
			 T4W = VSUB(T4U, T4V);
			 T1z = VFMA(LDK(KP923879532), T14, VMUL(LDK(KP382683432), T17));
			 {
			      Var T19, T1a, T1c, T1d;
			      T19 = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			      T1a = LD(&(xi[WS(is, 37)]), ivs, &(xi[WS(is, 1)]));
			      T1b = VSUB(T19, T1a);
			      T4R = VADD(T19, T1a);
			      T1c = LD(&(xi[WS(is, 21)]), ivs, &(xi[WS(is, 1)]));
			      T1d = LD(&(xi[WS(is, 53)]), ivs, &(xi[WS(is, 1)]));
			      T1e = VSUB(T1c, T1d);
			      T4S = VADD(T1c, T1d);
			 }
			 T1f = VFMA(LDK(KP382683432), T1b, VMUL(LDK(KP923879532), T1e));
			 T4T = VSUB(T4R, T4S);
			 T1y = VFNMS(LDK(KP382683432), T1e, VMUL(LDK(KP923879532), T1b));
		    }
		    {
			 Var T1g, T1r, T6k, T6l;
			 T1g = VSUB(T18, T1f);
			 T1r = VSUB(T1n, T1q);
			 T1s = VSUB(T1g, T1r);
			 T2K = VADD(T1r, T1g);
			 T6k = VADD(T4O, T4P);
			 T6l = VADD(T50, T4Z);
			 T6m = VADD(T6k, T6l);
			 T6Y = VSUB(T6k, T6l);
		    }
		    {
			 Var T6n, T6o, T1x, T1A;
			 T6n = VADD(T4R, T4S);
			 T6o = VADD(T4U, T4V);
			 T6p = VADD(T6n, T6o);
			 T6Z = VSUB(T6o, T6n);
			 T1x = VADD(T1v, T1w);
			 T1A = VADD(T1y, T1z);
			 T1B = VSUB(T1x, T1A);
			 T2L = VADD(T1x, T1A);
		    }
		    {
			 Var T3k, T3l, T4Q, T4X;
			 T3k = VSUB(T1v, T1w);
			 T3l = VADD(T1f, T18);
			 T3m = VADD(T3k, T3l);
			 T3Y = VSUB(T3k, T3l);
			 T4Q = VSUB(T4O, T4P);
			 T4X = VMUL(LDK(KP707106781), VADD(T4T, T4W));
			 T4Y = VADD(T4Q, T4X);
			 T5O = VSUB(T4Q, T4X);
		    }
		    {
			 Var T51, T52, T3n, T3o;
			 T51 = VSUB(T4Z, T50);
			 T52 = VMUL(LDK(KP707106781), VSUB(T4W, T4T));
			 T53 = VADD(T51, T52);
			 T5P = VSUB(T52, T51);
			 T3n = VADD(T1q, T1n);
			 T3o = VSUB(T1z, T1y);
			 T3p = VADD(T3n, T3o);
			 T3Z = VSUB(T3o, T3n);
		    }
	       }
	       {
		    Var T6N, T6R, T6Q, T6S;
		    {
			 Var T6L, T6M, T6O, T6P;
			 T6L = VADD(T6f, T6i);
			 T6M = VADD(T6F, T6C);
			 T6N = VADD(T6L, T6M);
			 T6R = VSUB(T6L, T6M);
			 T6O = VADD(T6m, T6p);
			 T6P = VADD(T6t, T6w);
			 T6Q = VADD(T6O, T6P);
			 T6S = VBYI(VSUB(T6P, T6O));
		    }
		    ST(&(xo[WS(os, 32)]), VSUB(T6N, T6Q), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 16)]), VADD(T6R, T6S), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(T6N, T6Q), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 48)]), VSUB(T6R, T6S), ovs, &(xo[0]));
	       }
	       {
		    Var T6j, T6G, T6y, T6H, T6q, T6x;
		    T6j = VSUB(T6f, T6i);
		    T6G = VSUB(T6C, T6F);
		    T6q = VSUB(T6m, T6p);
		    T6x = VSUB(T6t, T6w);
		    T6y = VMUL(LDK(KP707106781), VADD(T6q, T6x));
		    T6H = VMUL(LDK(KP707106781), VSUB(T6x, T6q));
		    {
			 Var T6z, T6I, T6J, T6K;
			 T6z = VADD(T6j, T6y);
			 T6I = VBYI(VADD(T6G, T6H));
			 ST(&(xo[WS(os, 56)]), VSUB(T6z, T6I), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 8)]), VADD(T6z, T6I), ovs, &(xo[0]));
			 T6J = VSUB(T6j, T6y);
			 T6K = VBYI(VSUB(T6H, T6G));
			 ST(&(xo[WS(os, 40)]), VSUB(T6J, T6K), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 24)]), VADD(T6J, T6K), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T6X, T7i, T78, T7g, T74, T7f, T7b, T7j, T6W, T77;
		    T6W = VMUL(LDK(KP707106781), VADD(T6U, T6V));
		    T6X = VADD(T6T, T6W);
		    T7i = VSUB(T6T, T6W);
		    T77 = VMUL(LDK(KP707106781), VSUB(T6V, T6U));
		    T78 = VADD(T76, T77);
		    T7g = VSUB(T77, T76);
		    {
			 Var T70, T73, T79, T7a;
			 T70 = VFMA(LDK(KP923879532), T6Y, VMUL(LDK(KP382683432), T6Z));
			 T73 = VFNMS(LDK(KP382683432), T72, VMUL(LDK(KP923879532), T71));
			 T74 = VADD(T70, T73);
			 T7f = VSUB(T73, T70);
			 T79 = VFNMS(LDK(KP382683432), T6Y, VMUL(LDK(KP923879532), T6Z));
			 T7a = VFMA(LDK(KP382683432), T71, VMUL(LDK(KP923879532), T72));
			 T7b = VADD(T79, T7a);
			 T7j = VSUB(T7a, T79);
		    }
		    {
			 Var T75, T7c, T7l, T7m;
			 T75 = VADD(T6X, T74);
			 T7c = VBYI(VADD(T78, T7b));
			 ST(&(xo[WS(os, 60)]), VSUB(T75, T7c), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 4)]), VADD(T75, T7c), ovs, &(xo[0]));
			 T7l = VBYI(VADD(T7g, T7f));
			 T7m = VADD(T7i, T7j);
			 ST(&(xo[WS(os, 12)]), VADD(T7l, T7m), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 52)]), VSUB(T7m, T7l), ovs, &(xo[0]));
		    }
		    {
			 Var T7d, T7e, T7h, T7k;
			 T7d = VSUB(T6X, T74);
			 T7e = VBYI(VSUB(T7b, T78));
			 ST(&(xo[WS(os, 36)]), VSUB(T7d, T7e), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 28)]), VADD(T7d, T7e), ovs, &(xo[0]));
			 T7h = VBYI(VSUB(T7f, T7g));
			 T7k = VSUB(T7i, T7j);
			 ST(&(xo[WS(os, 20)]), VADD(T7h, T7k), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 44)]), VSUB(T7k, T7h), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T5N, T68, T61, T69, T5U, T65, T5Y, T66;
		    {
			 Var T5L, T5M, T5Z, T60;
			 T5L = VSUB(T4p, T4w);
			 T5M = VSUB(T5u, T5t);
			 T5N = VADD(T5L, T5M);
			 T68 = VSUB(T5L, T5M);
			 T5Z = VFNMS(LDK(KP555570233), T5O, VMUL(LDK(KP831469612), T5P));
			 T60 = VFMA(LDK(KP555570233), T5R, VMUL(LDK(KP831469612), T5S));
			 T61 = VADD(T5Z, T60);
			 T69 = VSUB(T60, T5Z);
		    }
		    {
			 Var T5Q, T5T, T5W, T5X;
			 T5Q = VFMA(LDK(KP831469612), T5O, VMUL(LDK(KP555570233), T5P));
			 T5T = VFNMS(LDK(KP555570233), T5S, VMUL(LDK(KP831469612), T5R));
			 T5U = VADD(T5Q, T5T);
			 T65 = VSUB(T5T, T5Q);
			 T5W = VSUB(T5r, T5q);
			 T5X = VSUB(T4L, T4E);
			 T5Y = VADD(T5W, T5X);
			 T66 = VSUB(T5X, T5W);
		    }
		    {
			 Var T5V, T62, T6b, T6c;
			 T5V = VADD(T5N, T5U);
			 T62 = VBYI(VADD(T5Y, T61));
			 ST(&(xo[WS(os, 58)]), VSUB(T5V, T62), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 6)]), VADD(T5V, T62), ovs, &(xo[0]));
			 T6b = VBYI(VADD(T66, T65));
			 T6c = VADD(T68, T69);
			 ST(&(xo[WS(os, 10)]), VADD(T6b, T6c), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 54)]), VSUB(T6c, T6b), ovs, &(xo[0]));
		    }
		    {
			 Var T63, T64, T67, T6a;
			 T63 = VSUB(T5N, T5U);
			 T64 = VBYI(VSUB(T61, T5Y));
			 ST(&(xo[WS(os, 38)]), VSUB(T63, T64), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 26)]), VADD(T63, T64), ovs, &(xo[0]));
			 T67 = VBYI(VSUB(T65, T66));
			 T6a = VSUB(T68, T69);
			 ST(&(xo[WS(os, 22)]), VADD(T67, T6a), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 42)]), VSUB(T6a, T67), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T11, T2C, T2v, T2D, T2e, T2z, T2s, T2A;
		    {
			 Var Tr, T10, T2t, T2u;
			 Tr = VSUB(Tb, Tq);
			 T10 = VSUB(TI, TZ);
			 T11 = VADD(Tr, T10);
			 T2C = VSUB(Tr, T10);
			 T2t = VFNMS(LDK(KP634393284), T1B, VMUL(LDK(KP773010453), T1s));
			 T2u = VFMA(LDK(KP773010453), T2c, VMUL(LDK(KP634393284), T23));
			 T2v = VADD(T2t, T2u);
			 T2D = VSUB(T2u, T2t);
		    }
		    {
			 Var T1C, T2d, T2o, T2r;
			 T1C = VFMA(LDK(KP634393284), T1s, VMUL(LDK(KP773010453), T1B));
			 T2d = VFNMS(LDK(KP634393284), T2c, VMUL(LDK(KP773010453), T23));
			 T2e = VADD(T1C, T2d);
			 T2z = VSUB(T2d, T1C);
			 T2o = VSUB(T2i, T2n);
			 T2r = VSUB(T2p, T2q);
			 T2s = VADD(T2o, T2r);
			 T2A = VSUB(T2r, T2o);
		    }
		    {
			 Var T2f, T2w, T2F, T2G;
			 T2f = VADD(T11, T2e);
			 T2w = VBYI(VADD(T2s, T2v));
			 ST(&(xo[WS(os, 57)]), VSUB(T2f, T2w), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 7)]), VADD(T2f, T2w), ovs, &(xo[WS(os, 1)]));
			 T2F = VBYI(VADD(T2A, T2z));
			 T2G = VADD(T2C, T2D);
			 ST(&(xo[WS(os, 9)]), VADD(T2F, T2G), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 55)]), VSUB(T2G, T2F), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T2x, T2y, T2B, T2E;
			 T2x = VSUB(T11, T2e);
			 T2y = VBYI(VSUB(T2v, T2s));
			 ST(&(xo[WS(os, 39)]), VSUB(T2x, T2y), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 25)]), VADD(T2x, T2y), ovs, &(xo[WS(os, 1)]));
			 T2B = VBYI(VSUB(T2z, T2A));
			 T2E = VSUB(T2C, T2D);
			 ST(&(xo[WS(os, 23)]), VADD(T2B, T2E), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 41)]), VSUB(T2E, T2B), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T3j, T3Q, T3J, T3R, T3y, T3N, T3G, T3O;
		    {
			 Var T3b, T3i, T3H, T3I;
			 T3b = VADD(T39, T3a);
			 T3i = VADD(T3e, T3h);
			 T3j = VADD(T3b, T3i);
			 T3Q = VSUB(T3b, T3i);
			 T3H = VFNMS(LDK(KP290284677), T3m, VMUL(LDK(KP956940335), T3p));
			 T3I = VFMA(LDK(KP290284677), T3t, VMUL(LDK(KP956940335), T3w));
			 T3J = VADD(T3H, T3I);
			 T3R = VSUB(T3I, T3H);
		    }
		    {
			 Var T3q, T3x, T3C, T3F;
			 T3q = VFMA(LDK(KP956940335), T3m, VMUL(LDK(KP290284677), T3p));
			 T3x = VFNMS(LDK(KP290284677), T3w, VMUL(LDK(KP956940335), T3t));
			 T3y = VADD(T3q, T3x);
			 T3N = VSUB(T3x, T3q);
			 T3C = VADD(T3A, T3B);
			 T3F = VADD(T3D, T3E);
			 T3G = VADD(T3C, T3F);
			 T3O = VSUB(T3F, T3C);
		    }
		    {
			 Var T3z, T3K, T3T, T3U;
			 T3z = VADD(T3j, T3y);
			 T3K = VBYI(VADD(T3G, T3J));
			 ST(&(xo[WS(os, 61)]), VSUB(T3z, T3K), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 3)]), VADD(T3z, T3K), ovs, &(xo[WS(os, 1)]));
			 T3T = VBYI(VADD(T3O, T3N));
			 T3U = VADD(T3Q, T3R);
			 ST(&(xo[WS(os, 13)]), VADD(T3T, T3U), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 51)]), VSUB(T3U, T3T), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T3L, T3M, T3P, T3S;
			 T3L = VSUB(T3j, T3y);
			 T3M = VBYI(VSUB(T3J, T3G));
			 ST(&(xo[WS(os, 35)]), VSUB(T3L, T3M), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 29)]), VADD(T3L, T3M), ovs, &(xo[WS(os, 1)]));
			 T3P = VBYI(VSUB(T3N, T3O));
			 T3S = VSUB(T3Q, T3R);
			 ST(&(xo[WS(os, 19)]), VADD(T3P, T3S), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 45)]), VSUB(T3S, T3P), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T4N, T5G, T5z, T5H, T5m, T5D, T5w, T5E;
		    {
			 Var T4x, T4M, T5x, T5y;
			 T4x = VADD(T4p, T4w);
			 T4M = VADD(T4E, T4L);
			 T4N = VADD(T4x, T4M);
			 T5G = VSUB(T4x, T4M);
			 T5x = VFNMS(LDK(KP195090322), T4Y, VMUL(LDK(KP980785280), T53));
			 T5y = VFMA(LDK(KP195090322), T5f, VMUL(LDK(KP980785280), T5k));
			 T5z = VADD(T5x, T5y);
			 T5H = VSUB(T5y, T5x);
		    }
		    {
			 Var T54, T5l, T5s, T5v;
			 T54 = VFMA(LDK(KP980785280), T4Y, VMUL(LDK(KP195090322), T53));
			 T5l = VFNMS(LDK(KP195090322), T5k, VMUL(LDK(KP980785280), T5f));
			 T5m = VADD(T54, T5l);
			 T5D = VSUB(T5l, T54);
			 T5s = VADD(T5q, T5r);
			 T5v = VADD(T5t, T5u);
			 T5w = VADD(T5s, T5v);
			 T5E = VSUB(T5v, T5s);
		    }
		    {
			 Var T5n, T5A, T5J, T5K;
			 T5n = VADD(T4N, T5m);
			 T5A = VBYI(VADD(T5w, T5z));
			 ST(&(xo[WS(os, 62)]), VSUB(T5n, T5A), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 2)]), VADD(T5n, T5A), ovs, &(xo[0]));
			 T5J = VBYI(VADD(T5E, T5D));
			 T5K = VADD(T5G, T5H);
			 ST(&(xo[WS(os, 14)]), VADD(T5J, T5K), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 50)]), VSUB(T5K, T5J), ovs, &(xo[0]));
		    }
		    {
			 Var T5B, T5C, T5F, T5I;
			 T5B = VSUB(T4N, T5m);
			 T5C = VBYI(VSUB(T5z, T5w));
			 ST(&(xo[WS(os, 34)]), VSUB(T5B, T5C), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 30)]), VADD(T5B, T5C), ovs, &(xo[0]));
			 T5F = VBYI(VSUB(T5D, T5E));
			 T5I = VSUB(T5G, T5H);
			 ST(&(xo[WS(os, 18)]), VADD(T5F, T5I), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 46)]), VSUB(T5I, T5F), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T2J, T34, T2X, T35, T2Q, T31, T2U, T32;
		    {
			 Var T2H, T2I, T2V, T2W;
			 T2H = VADD(Tb, Tq);
			 T2I = VADD(T2q, T2p);
			 T2J = VADD(T2H, T2I);
			 T34 = VSUB(T2H, T2I);
			 T2V = VFNMS(LDK(KP098017140), T2L, VMUL(LDK(KP995184726), T2K));
			 T2W = VFMA(LDK(KP995184726), T2O, VMUL(LDK(KP098017140), T2N));
			 T2X = VADD(T2V, T2W);
			 T35 = VSUB(T2W, T2V);
		    }
		    {
			 Var T2M, T2P, T2S, T2T;
			 T2M = VFMA(LDK(KP098017140), T2K, VMUL(LDK(KP995184726), T2L));
			 T2P = VFNMS(LDK(KP098017140), T2O, VMUL(LDK(KP995184726), T2N));
			 T2Q = VADD(T2M, T2P);
			 T31 = VSUB(T2P, T2M);
			 T2S = VADD(T2n, T2i);
			 T2T = VADD(TZ, TI);
			 T2U = VADD(T2S, T2T);
			 T32 = VSUB(T2T, T2S);
		    }
		    {
			 Var T2R, T2Y, T37, T38;
			 T2R = VADD(T2J, T2Q);
			 T2Y = VBYI(VADD(T2U, T2X));
			 ST(&(xo[WS(os, 63)]), VSUB(T2R, T2Y), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 1)]), VADD(T2R, T2Y), ovs, &(xo[WS(os, 1)]));
			 T37 = VBYI(VADD(T32, T31));
			 T38 = VADD(T34, T35);
			 ST(&(xo[WS(os, 15)]), VADD(T37, T38), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 49)]), VSUB(T38, T37), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T2Z, T30, T33, T36;
			 T2Z = VSUB(T2J, T2Q);
			 T30 = VBYI(VSUB(T2X, T2U));
			 ST(&(xo[WS(os, 33)]), VSUB(T2Z, T30), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 31)]), VADD(T2Z, T30), ovs, &(xo[WS(os, 1)]));
			 T33 = VBYI(VSUB(T31, T32));
			 T36 = VSUB(T34, T35);
			 ST(&(xo[WS(os, 17)]), VADD(T33, T36), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 47)]), VSUB(T36, T33), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T3X, T4i, T4b, T4j, T44, T4f, T48, T4g;
		    {
			 Var T3V, T3W, T49, T4a;
			 T3V = VSUB(T39, T3a);
			 T3W = VSUB(T3E, T3D);
			 T3X = VADD(T3V, T3W);
			 T4i = VSUB(T3V, T3W);
			 T49 = VFNMS(LDK(KP471396736), T3Y, VMUL(LDK(KP881921264), T3Z));
			 T4a = VFMA(LDK(KP471396736), T41, VMUL(LDK(KP881921264), T42));
			 T4b = VADD(T49, T4a);
			 T4j = VSUB(T4a, T49);
		    }
		    {
			 Var T40, T43, T46, T47;
			 T40 = VFMA(LDK(KP881921264), T3Y, VMUL(LDK(KP471396736), T3Z));
			 T43 = VFNMS(LDK(KP471396736), T42, VMUL(LDK(KP881921264), T41));
			 T44 = VADD(T40, T43);
			 T4f = VSUB(T43, T40);
			 T46 = VSUB(T3B, T3A);
			 T47 = VSUB(T3h, T3e);
			 T48 = VADD(T46, T47);
			 T4g = VSUB(T47, T46);
		    }
		    {
			 Var T45, T4c, T4l, T4m;
			 T45 = VADD(T3X, T44);
			 T4c = VBYI(VADD(T48, T4b));
			 ST(&(xo[WS(os, 59)]), VSUB(T45, T4c), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 5)]), VADD(T45, T4c), ovs, &(xo[WS(os, 1)]));
			 T4l = VBYI(VADD(T4g, T4f));
			 T4m = VADD(T4i, T4j);
			 ST(&(xo[WS(os, 11)]), VADD(T4l, T4m), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 53)]), VSUB(T4m, T4l), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T4d, T4e, T4h, T4k;
			 T4d = VSUB(T3X, T44);
			 T4e = VBYI(VSUB(T4b, T48));
			 ST(&(xo[WS(os, 37)]), VSUB(T4d, T4e), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 27)]), VADD(T4d, T4e), ovs, &(xo[WS(os, 1)]));
			 T4h = VBYI(VSUB(T4f, T4g));
			 T4k = VSUB(T4i, T4j);
			 ST(&(xo[WS(os, 21)]), VADD(T4h, T4k), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 43)]), VSUB(T4k, T4h), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
static void n1fv_128(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DVK(KP941544065, +0.941544065183020778412509402599502357185589796);
     DVK(KP336889853, +0.336889853392220050689253212619147570477766780);
     DVK(KP903989293, +0.903989293123443331586200297230537048710132025);
     DVK(KP427555093, +0.427555093430282094320966856888798534304578629);
     DVK(KP970031253, +0.970031253194543992603984207286100251456865962);
     DVK(KP242980179, +0.242980179903263889948274162077471118320990783);
     DVK(KP857728610, +0.857728610000272069902269984284770137042490799);
     DVK(KP514102744, +0.514102744193221726593693838968815772608049120);
     DVK(KP671558954, +0.671558954847018400625376850427421803228750632);
     DVK(KP740951125, +0.740951125354959091175616897495162729728955309);
     DVK(KP049067674, +0.049067674327418014254954976942682658314745363);
     DVK(KP998795456, +0.998795456205172392714771604759100694443203615);
     DVK(KP595699304, +0.595699304492433343467036528829969889511926338);
     DVK(KP803207531, +0.803207531480644909806676512963141923879569427);
     DVK(KP146730474, +0.146730474455361751658850129646717819706215317);
     DVK(KP989176509, +0.989176509964780973451673738016243063983689533);
     DVK(KP956940335, +0.956940335732208864935797886980269969482849206);
     DVK(KP290284677, +0.290284677254462367636192375817395274691476278);
     DVK(KP881921264, +0.881921264348355029712756863660388349508442621);
     DVK(KP471396736, +0.471396736825997648556387625905254377657460319);
     DVK(KP634393284, +0.634393284163645498215171613225493370675687095);
     DVK(KP773010453, +0.773010453362736960810906609758469800971041293);
     DVK(KP098017140, +0.098017140329560601994195563888641845861136673);
     DVK(KP995184726, +0.995184726672196886244836953109479921575474869);
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
	  INT i;
	  const R *xi;
	  R *xo;
	  xi = ri;
	  xo = ro;
	  for (i = v; i > 0; i = i - VL, xi = xi + (VL * ivs), xo = xo + (VL * ovs), MAKE_VOLATILE_STRIDE(256, is), MAKE_VOLATILE_STRIDE(256, os)) {
	       Var Tr, T5J, Ted, Tgf, Tfq, TgH, T4U, T6b, T6Z, T8T, Tad, TcZ, Tcc, Td0, T84;
	       Var T9l, Tb6, Tbt, T2G, T5X, TeV, Tgr, T3p, T5V, T7B, T95, TeK, Tgt, T7q, T97;
	       Var Td8, TdK, TbD, Tc0, T3V, T61, Tfg, TgB, T4E, T65, T7U, T9f, Tf5, Tgx, T7J;
	       Var T9b, Tdf, TdN, Td2, Td3, TI, T4V, Tft, Tgg, TZ, T4W, T75, T86, Tek, TgG;
	       Var T72, T85, Tas, Tcd, Tdp, Tdq, TdG, Teq, Tgm, Tet, Tgl, T1s, T5P, T1B, T5Q;
	       Var T7d, T8Z, TaI, Tcf, T7a, T90, Tdm, Tdn, TdH, Tez, Tgi, TeC, Tgj, T23, T5N;
	       Var T2c, T5M, T7k, T8X, TaX, Tcg, T7h, T8W, Tbl, Tbu, Tdb, TdL, TeY, Tgu, TeR;
	       Var Tgq, T7x, T98, T7E, T94, T3f, T5Y, T3s, T5U, TbS, Tc1, Tdi, TdO, Tfj, Tgy;
	       Var Tfc, TgA, T7Q, T9e, T7X, T9c, T4u, T64, T4H, T62;
	       {
		    Var T3, Ta7, T4P, Ta8, Ta, Tab, T4M, Taa, Tc9, Tca, Ti, Tea, T4S, Tc6, Tc7;
		    Var Tp, Teb, T4R;
		    {
			 Var T1, T2, T4N, T4O;
			 T1 = LD(&(xi[0]), ivs, &(xi[0]));
			 T2 = LD(&(xi[WS(is, 64)]), ivs, &(xi[0]));
			 T3 = VSUB(T1, T2);
			 Ta7 = VADD(T1, T2);
			 T4N = LD(&(xi[WS(is, 32)]), ivs, &(xi[0]));
			 T4O = LD(&(xi[WS(is, 96)]), ivs, &(xi[0]));
			 T4P = VSUB(T4N, T4O);
			 Ta8 = VADD(T4N, T4O);
		    }
		    {
			 Var T4, T5, T6, T7, T8, T9;
			 T4 = LD(&(xi[WS(is, 16)]), ivs, &(xi[0]));
			 T5 = LD(&(xi[WS(is, 80)]), ivs, &(xi[0]));
			 T6 = VSUB(T4, T5);
			 T7 = LD(&(xi[WS(is, 112)]), ivs, &(xi[0]));
			 T8 = LD(&(xi[WS(is, 48)]), ivs, &(xi[0]));
			 T9 = VSUB(T7, T8);
			 Ta = VMUL(LDK(KP707106781), VADD(T6, T9));
			 Tab = VADD(T7, T8);
			 T4M = VMUL(LDK(KP707106781), VSUB(T9, T6));
			 Taa = VADD(T4, T5);
		    }
		    {
			 Var Te, Th, Tl, To;
			 {
			      Var Tc, Td, Tf, Tg;
			      Tc = LD(&(xi[WS(is, 8)]), ivs, &(xi[0]));
			      Td = LD(&(xi[WS(is, 72)]), ivs, &(xi[0]));
			      Te = VSUB(Tc, Td);
			      Tc9 = VADD(Tc, Td);
			      Tf = LD(&(xi[WS(is, 40)]), ivs, &(xi[0]));
			      Tg = LD(&(xi[WS(is, 104)]), ivs, &(xi[0]));
			      Th = VSUB(Tf, Tg);
			      Tca = VADD(Tf, Tg);
			 }
			 Ti = VFNMS(LDK(KP382683432), Th, VMUL(LDK(KP923879532), Te));
			 Tea = VSUB(Tc9, Tca);
			 T4S = VFMA(LDK(KP382683432), Te, VMUL(LDK(KP923879532), Th));
			 {
			      Var Tj, Tk, Tm, Tn;
			      Tj = LD(&(xi[WS(is, 120)]), ivs, &(xi[0]));
			      Tk = LD(&(xi[WS(is, 56)]), ivs, &(xi[0]));
			      Tl = VSUB(Tj, Tk);
			      Tc6 = VADD(Tj, Tk);
			      Tm = LD(&(xi[WS(is, 24)]), ivs, &(xi[0]));
			      Tn = LD(&(xi[WS(is, 88)]), ivs, &(xi[0]));
			      To = VSUB(Tm, Tn);
			      Tc7 = VADD(Tm, Tn);
			 }
			 Tp = VFMA(LDK(KP923879532), Tl, VMUL(LDK(KP382683432), To));
			 Teb = VSUB(Tc6, Tc7);
			 T4R = VFNMS(LDK(KP923879532), To, VMUL(LDK(KP382683432), Tl));
		    }
		    {
			 Var Tb, Tq, Te9, Tec;
			 Tb = VADD(T3, Ta);
			 Tq = VADD(Ti, Tp);
			 Tr = VADD(Tb, Tq);
			 T5J = VSUB(Tb, Tq);
			 Te9 = VSUB(Ta7, Ta8);
			 Tec = VMUL(LDK(KP707106781), VADD(Tea, Teb));
			 Ted = VADD(Te9, Tec);
			 Tgf = VSUB(Te9, Tec);
		    }
		    {
			 Var Tfo, Tfp, T4Q, T4T;
			 Tfo = VSUB(Tab, Taa);
			 Tfp = VMUL(LDK(KP707106781), VSUB(Teb, Tea));
			 Tfq = VADD(Tfo, Tfp);
			 TgH = VSUB(Tfp, Tfo);
			 T4Q = VSUB(T4M, T4P);
			 T4T = VSUB(T4R, T4S);
			 T4U = VADD(T4Q, T4T);
			 T6b = VSUB(T4T, T4Q);
		    }
		    {
			 Var T6X, T6Y, Ta9, Tac;
			 T6X = VSUB(T3, Ta);
			 T6Y = VADD(T4S, T4R);
			 T6Z = VADD(T6X, T6Y);
			 T8T = VSUB(T6X, T6Y);
			 Ta9 = VADD(Ta7, Ta8);
			 Tac = VADD(Taa, Tab);
			 Tad = VSUB(Ta9, Tac);
			 TcZ = VADD(Ta9, Tac);
		    }
		    {
			 Var Tc8, Tcb, T82, T83;
			 Tc8 = VADD(Tc6, Tc7);
			 Tcb = VADD(Tc9, Tca);
			 Tcc = VSUB(Tc8, Tcb);
			 Td0 = VADD(Tcb, Tc8);
			 T82 = VADD(T4P, T4M);
			 T83 = VSUB(Tp, Ti);
			 T84 = VADD(T82, T83);
			 T9l = VSUB(T83, T82);
		    }
	       }
	       {
		    Var Tb0, Tb1, T2i, Tb2, T3k, Tb3, Tb4, T2p, Tb5, T3h, T2x, TeH, T3n, Tbs, T2E;
		    Var TeI, T3m, Tbp, T2l, T2o, TeG, TeJ;
		    {
			 Var T2g, T2h, T3i, T3j;
			 T2g = LD(&(xi[WS(is, 1)]), ivs, &(xi[WS(is, 1)]));
			 T2h = LD(&(xi[WS(is, 65)]), ivs, &(xi[WS(is, 1)]));
			 Tb0 = VADD(T2g, T2h);
			 T3i = LD(&(xi[WS(is, 33)]), ivs, &(xi[WS(is, 1)]));
			 T3j = LD(&(xi[WS(is, 97)]), ivs, &(xi[WS(is, 1)]));
			 Tb1 = VADD(T3i, T3j);
			 T2i = VSUB(T2g, T2h);
			 Tb2 = VADD(Tb0, Tb1);
			 T3k = VSUB(T3i, T3j);
		    }
		    {
			 Var T2j, T2k, T2m, T2n;
			 T2j = LD(&(xi[WS(is, 17)]), ivs, &(xi[WS(is, 1)]));
			 T2k = LD(&(xi[WS(is, 81)]), ivs, &(xi[WS(is, 1)]));
			 T2l = VSUB(T2j, T2k);
			 Tb3 = VADD(T2j, T2k);
			 T2m = LD(&(xi[WS(is, 113)]), ivs, &(xi[WS(is, 1)]));
			 T2n = LD(&(xi[WS(is, 49)]), ivs, &(xi[WS(is, 1)]));
			 T2o = VSUB(T2m, T2n);
			 Tb4 = VADD(T2m, T2n);
		    }
		    T2p = VMUL(LDK(KP707106781), VADD(T2l, T2o));
		    Tb5 = VADD(Tb3, Tb4);
		    T3h = VMUL(LDK(KP707106781), VSUB(T2o, T2l));
		    {
			 Var T2t, Tbq, T2w, Tbr;
			 {
			      Var T2r, T2s, T2u, T2v;
			      T2r = LD(&(xi[WS(is, 9)]), ivs, &(xi[WS(is, 1)]));
			      T2s = LD(&(xi[WS(is, 73)]), ivs, &(xi[WS(is, 1)]));
			      T2t = VSUB(T2r, T2s);
			      Tbq = VADD(T2r, T2s);
			      T2u = LD(&(xi[WS(is, 41)]), ivs, &(xi[WS(is, 1)]));
			      T2v = LD(&(xi[WS(is, 105)]), ivs, &(xi[WS(is, 1)]));
			      T2w = VSUB(T2u, T2v);
			      Tbr = VADD(T2u, T2v);
			 }
			 T2x = VFNMS(LDK(KP382683432), T2w, VMUL(LDK(KP923879532), T2t));
			 TeH = VSUB(Tbq, Tbr);
			 T3n = VFMA(LDK(KP382683432), T2t, VMUL(LDK(KP923879532), T2w));
			 Tbs = VADD(Tbq, Tbr);
		    }
		    {
			 Var T2A, Tbn, T2D, Tbo;
			 {
			      Var T2y, T2z, T2B, T2C;
			      T2y = LD(&(xi[WS(is, 121)]), ivs, &(xi[WS(is, 1)]));
			      T2z = LD(&(xi[WS(is, 57)]), ivs, &(xi[WS(is, 1)]));
			      T2A = VSUB(T2y, T2z);
			      Tbn = VADD(T2y, T2z);
			      T2B = LD(&(xi[WS(is, 25)]), ivs, &(xi[WS(is, 1)]));
			      T2C = LD(&(xi[WS(is, 89)]), ivs, &(xi[WS(is, 1)]));
			      T2D = VSUB(T2B, T2C);
			      Tbo = VADD(T2B, T2C);
			 }
			 T2E = VFMA(LDK(KP923879532), T2A, VMUL(LDK(KP382683432), T2D));
			 TeI = VSUB(Tbn, Tbo);
			 T3m = VFNMS(LDK(KP923879532), T2D, VMUL(LDK(KP382683432), T2A));
			 Tbp = VADD(Tbn, Tbo);
		    }
		    Tb6 = VSUB(Tb2, Tb5);
		    Tbt = VSUB(Tbp, Tbs);
		    {
			 Var T2q, T2F, TeT, TeU;
			 T2q = VADD(T2i, T2p);
			 T2F = VADD(T2x, T2E);
			 T2G = VADD(T2q, T2F);
			 T5X = VSUB(T2q, T2F);
			 TeT = VSUB(Tb4, Tb3);
			 TeU = VMUL(LDK(KP707106781), VSUB(TeI, TeH));
			 TeV = VADD(TeT, TeU);
			 Tgr = VSUB(TeU, TeT);
		    }
		    {
			 Var T3l, T3o, T7z, T7A;
			 T3l = VSUB(T3h, T3k);
			 T3o = VSUB(T3m, T3n);
			 T3p = VADD(T3l, T3o);
			 T5V = VSUB(T3o, T3l);
			 T7z = VADD(T3k, T3h);
			 T7A = VSUB(T2E, T2x);
			 T7B = VADD(T7z, T7A);
			 T95 = VSUB(T7A, T7z);
		    }
		    TeG = VSUB(Tb0, Tb1);
		    TeJ = VMUL(LDK(KP707106781), VADD(TeH, TeI));
		    TeK = VADD(TeG, TeJ);
		    Tgt = VSUB(TeG, TeJ);
		    {
			 Var T7o, T7p, Td6, Td7;
			 T7o = VSUB(T2i, T2p);
			 T7p = VADD(T3n, T3m);
			 T7q = VADD(T7o, T7p);
			 T97 = VSUB(T7o, T7p);
			 Td6 = VADD(Tb2, Tb5);
			 Td7 = VADD(Tbs, Tbp);
			 Td8 = VADD(Td6, Td7);
			 TdK = VSUB(Td6, Td7);
		    }
	       }
	       {
		    Var Tbx, Tby, T3x, Tbz, T4z, TbA, TbB, T3E, TbC, T4w, T3M, Tf2, T4C, TbZ, T3T;
		    Var Tf3, T4B, TbW, T3A, T3D, Tf1, Tf4;
		    {
			 Var T3v, T3w, T4x, T4y;
			 T3v = LD(&(xi[WS(is, 127)]), ivs, &(xi[WS(is, 1)]));
			 T3w = LD(&(xi[WS(is, 63)]), ivs, &(xi[WS(is, 1)]));
			 Tbx = VADD(T3v, T3w);
			 T4x = LD(&(xi[WS(is, 31)]), ivs, &(xi[WS(is, 1)]));
			 T4y = LD(&(xi[WS(is, 95)]), ivs, &(xi[WS(is, 1)]));
			 Tby = VADD(T4x, T4y);
			 T3x = VSUB(T3v, T3w);
			 Tbz = VADD(Tbx, Tby);
			 T4z = VSUB(T4x, T4y);
		    }
		    {
			 Var T3y, T3z, T3B, T3C;
			 T3y = LD(&(xi[WS(is, 15)]), ivs, &(xi[WS(is, 1)]));
			 T3z = LD(&(xi[WS(is, 79)]), ivs, &(xi[WS(is, 1)]));
			 T3A = VSUB(T3y, T3z);
			 TbA = VADD(T3y, T3z);
			 T3B = LD(&(xi[WS(is, 111)]), ivs, &(xi[WS(is, 1)]));
			 T3C = LD(&(xi[WS(is, 47)]), ivs, &(xi[WS(is, 1)]));
			 T3D = VSUB(T3B, T3C);
			 TbB = VADD(T3B, T3C);
		    }
		    T3E = VMUL(LDK(KP707106781), VADD(T3A, T3D));
		    TbC = VADD(TbA, TbB);
		    T4w = VMUL(LDK(KP707106781), VSUB(T3D, T3A));
		    {
			 Var T3I, TbX, T3L, TbY;
			 {
			      Var T3G, T3H, T3J, T3K;
			      T3G = LD(&(xi[WS(is, 7)]), ivs, &(xi[WS(is, 1)]));
			      T3H = LD(&(xi[WS(is, 71)]), ivs, &(xi[WS(is, 1)]));
			      T3I = VSUB(T3G, T3H);
			      TbX = VADD(T3G, T3H);
			      T3J = LD(&(xi[WS(is, 39)]), ivs, &(xi[WS(is, 1)]));
			      T3K = LD(&(xi[WS(is, 103)]), ivs, &(xi[WS(is, 1)]));
			      T3L = VSUB(T3J, T3K);
			      TbY = VADD(T3J, T3K);
			 }
			 T3M = VFNMS(LDK(KP382683432), T3L, VMUL(LDK(KP923879532), T3I));
			 Tf2 = VSUB(TbX, TbY);
			 T4C = VFMA(LDK(KP382683432), T3I, VMUL(LDK(KP923879532), T3L));
			 TbZ = VADD(TbX, TbY);
		    }
		    {
			 Var T3P, TbU, T3S, TbV;
			 {
			      Var T3N, T3O, T3Q, T3R;
			      T3N = LD(&(xi[WS(is, 119)]), ivs, &(xi[WS(is, 1)]));
			      T3O = LD(&(xi[WS(is, 55)]), ivs, &(xi[WS(is, 1)]));
			      T3P = VSUB(T3N, T3O);
			      TbU = VADD(T3N, T3O);
			      T3Q = LD(&(xi[WS(is, 23)]), ivs, &(xi[WS(is, 1)]));
			      T3R = LD(&(xi[WS(is, 87)]), ivs, &(xi[WS(is, 1)]));
			      T3S = VSUB(T3Q, T3R);
			      TbV = VADD(T3Q, T3R);
			 }
			 T3T = VFMA(LDK(KP923879532), T3P, VMUL(LDK(KP382683432), T3S));
			 Tf3 = VSUB(TbU, TbV);
			 T4B = VFNMS(LDK(KP923879532), T3S, VMUL(LDK(KP382683432), T3P));
			 TbW = VADD(TbU, TbV);
		    }
		    TbD = VSUB(Tbz, TbC);
		    Tc0 = VSUB(TbW, TbZ);
		    {
			 Var T3F, T3U, Tfe, Tff;
			 T3F = VADD(T3x, T3E);
			 T3U = VADD(T3M, T3T);
			 T3V = VADD(T3F, T3U);
			 T61 = VSUB(T3F, T3U);
			 Tfe = VSUB(TbB, TbA);
			 Tff = VMUL(LDK(KP707106781), VSUB(Tf3, Tf2));
			 Tfg = VADD(Tfe, Tff);
			 TgB = VSUB(Tff, Tfe);
		    }
		    {
			 Var T4A, T4D, T7S, T7T;
			 T4A = VSUB(T4w, T4z);
			 T4D = VSUB(T4B, T4C);
			 T4E = VADD(T4A, T4D);
			 T65 = VSUB(T4D, T4A);
			 T7S = VADD(T4z, T4w);
			 T7T = VSUB(T3T, T3M);
			 T7U = VADD(T7S, T7T);
			 T9f = VSUB(T7T, T7S);
		    }
		    Tf1 = VSUB(Tbx, Tby);
		    Tf4 = VMUL(LDK(KP707106781), VADD(Tf2, Tf3));
		    Tf5 = VADD(Tf1, Tf4);
		    Tgx = VSUB(Tf1, Tf4);
		    {
			 Var T7H, T7I, Tdd, Tde;
			 T7H = VSUB(T3x, T3E);
			 T7I = VADD(T4C, T4B);
			 T7J = VADD(T7H, T7I);
			 T9b = VSUB(T7H, T7I);
			 Tdd = VADD(Tbz, TbC);
			 Tde = VADD(TbZ, TbW);
			 Tdf = VADD(Tdd, Tde);
			 TdN = VSUB(Tdd, Tde);
		    }
	       }
	       {
		    Var Tu, Tee, TG, Tag, TL, Teh, TX, Tan, TB, Tef, TD, Taj, TS, Tei, TU;
		    Var Taq, Teg, Tej;
		    {
			 Var Ts, Tt, Tae, TE, TF, Taf;
			 Ts = LD(&(xi[WS(is, 4)]), ivs, &(xi[0]));
			 Tt = LD(&(xi[WS(is, 68)]), ivs, &(xi[0]));
			 Tae = VADD(Ts, Tt);
			 TE = LD(&(xi[WS(is, 36)]), ivs, &(xi[0]));
			 TF = LD(&(xi[WS(is, 100)]), ivs, &(xi[0]));
			 Taf = VADD(TE, TF);
			 Tu = VSUB(Ts, Tt);
			 Tee = VSUB(Tae, Taf);
			 TG = VSUB(TE, TF);
			 Tag = VADD(Tae, Taf);
		    }
		    {
			 Var TJ, TK, Tal, TV, TW, Tam;
			 TJ = LD(&(xi[WS(is, 124)]), ivs, &(xi[0]));
			 TK = LD(&(xi[WS(is, 60)]), ivs, &(xi[0]));
			 Tal = VADD(TJ, TK);
			 TV = LD(&(xi[WS(is, 28)]), ivs, &(xi[0]));
			 TW = LD(&(xi[WS(is, 92)]), ivs, &(xi[0]));
			 Tam = VADD(TV, TW);
			 TL = VSUB(TJ, TK);
			 Teh = VSUB(Tal, Tam);
			 TX = VSUB(TV, TW);
			 Tan = VADD(Tal, Tam);
		    }
		    {
			 Var Tx, Tah, TA, Tai;
			 {
			      Var Tv, Tw, Ty, Tz;
			      Tv = LD(&(xi[WS(is, 20)]), ivs, &(xi[0]));
			      Tw = LD(&(xi[WS(is, 84)]), ivs, &(xi[0]));
			      Tx = VSUB(Tv, Tw);
			      Tah = VADD(Tv, Tw);
			      Ty = LD(&(xi[WS(is, 116)]), ivs, &(xi[0]));
			      Tz = LD(&(xi[WS(is, 52)]), ivs, &(xi[0]));
			      TA = VSUB(Ty, Tz);
			      Tai = VADD(Ty, Tz);
			 }
			 TB = VMUL(LDK(KP707106781), VADD(Tx, TA));
			 Tef = VSUB(Tai, Tah);
			 TD = VMUL(LDK(KP707106781), VSUB(TA, Tx));
			 Taj = VADD(Tah, Tai);
		    }
		    {
			 Var TO, Tao, TR, Tap;
			 {
			      Var TM, TN, TP, TQ;
			      TM = LD(&(xi[WS(is, 12)]), ivs, &(xi[0]));
			      TN = LD(&(xi[WS(is, 76)]), ivs, &(xi[0]));
			      TO = VSUB(TM, TN);
			      Tao = VADD(TM, TN);
			      TP = LD(&(xi[WS(is, 108)]), ivs, &(xi[0]));
			      TQ = LD(&(xi[WS(is, 44)]), ivs, &(xi[0]));
			      TR = VSUB(TP, TQ);
			      Tap = VADD(TP, TQ);
			 }
			 TS = VMUL(LDK(KP707106781), VADD(TO, TR));
			 Tei = VSUB(Tap, Tao);
			 TU = VMUL(LDK(KP707106781), VSUB(TR, TO));
			 Taq = VADD(Tao, Tap);
		    }
		    Td2 = VADD(Tag, Taj);
		    Td3 = VADD(Tan, Taq);
		    {
			 Var TC, TH, Tfr, Tfs;
			 TC = VADD(Tu, TB);
			 TH = VSUB(TD, TG);
			 TI = VFMA(LDK(KP980785280), TC, VMUL(LDK(KP195090322), TH));
			 T4V = VFNMS(LDK(KP195090322), TC, VMUL(LDK(KP980785280), TH));
			 Tfr = VFNMS(LDK(KP382683432), Tee, VMUL(LDK(KP923879532), Tef));
			 Tfs = VFMA(LDK(KP382683432), Teh, VMUL(LDK(KP923879532), Tei));
			 Tft = VADD(Tfr, Tfs);
			 Tgg = VSUB(Tfs, Tfr);
		    }
		    {
			 Var TT, TY, T73, T74;
			 TT = VADD(TL, TS);
			 TY = VSUB(TU, TX);
			 TZ = VFNMS(LDK(KP195090322), TY, VMUL(LDK(KP980785280), TT));
			 T4W = VFMA(LDK(KP195090322), TT, VMUL(LDK(KP980785280), TY));
			 T73 = VSUB(TL, TS);
			 T74 = VADD(TX, TU);
			 T75 = VFNMS(LDK(KP555570233), T74, VMUL(LDK(KP831469612), T73));
			 T86 = VFMA(LDK(KP555570233), T73, VMUL(LDK(KP831469612), T74));
		    }
		    Teg = VFMA(LDK(KP923879532), Tee, VMUL(LDK(KP382683432), Tef));
		    Tej = VFNMS(LDK(KP382683432), Tei, VMUL(LDK(KP923879532), Teh));
		    Tek = VADD(Teg, Tej);
		    TgG = VSUB(Tej, Teg);
		    {
			 Var T70, T71, Tak, Tar;
			 T70 = VSUB(Tu, TB);
			 T71 = VADD(TG, TD);
			 T72 = VFMA(LDK(KP831469612), T70, VMUL(LDK(KP555570233), T71));
			 T85 = VFNMS(LDK(KP555570233), T70, VMUL(LDK(KP831469612), T71));
			 Tak = VSUB(Tag, Taj);
			 Tar = VSUB(Tan, Taq);
			 Tas = VMUL(LDK(KP707106781), VADD(Tak, Tar));
			 Tcd = VMUL(LDK(KP707106781), VSUB(Tar, Tak));
		    }
	       }
	       {
		    Var Tav, Tau, T1b, Taw, T1v, Tay, Tax, T18, Taz, T1w, T1j, Teo, T1z, TaD, T1q;
		    Var Ten, T1y, TaG, T14, T17, Tem, Tep;
		    {
			 Var T19, T1a, T1t, T1u;
			 T19 = LD(&(xi[WS(is, 34)]), ivs, &(xi[0]));
			 T1a = LD(&(xi[WS(is, 98)]), ivs, &(xi[0]));
			 Tav = VADD(T19, T1a);
			 T1t = LD(&(xi[WS(is, 2)]), ivs, &(xi[0]));
			 T1u = LD(&(xi[WS(is, 66)]), ivs, &(xi[0]));
			 Tau = VADD(T1t, T1u);
			 T1b = VSUB(T19, T1a);
			 Taw = VADD(Tau, Tav);
			 T1v = VSUB(T1t, T1u);
		    }
		    {
			 Var T12, T13, T15, T16;
			 T12 = LD(&(xi[WS(is, 114)]), ivs, &(xi[0]));
			 T13 = LD(&(xi[WS(is, 50)]), ivs, &(xi[0]));
			 T14 = VSUB(T12, T13);
			 Tay = VADD(T12, T13);
			 T15 = LD(&(xi[WS(is, 18)]), ivs, &(xi[0]));
			 T16 = LD(&(xi[WS(is, 82)]), ivs, &(xi[0]));
			 T17 = VSUB(T15, T16);
			 Tax = VADD(T15, T16);
		    }
		    T18 = VMUL(LDK(KP707106781), VSUB(T14, T17));
		    Taz = VADD(Tax, Tay);
		    T1w = VMUL(LDK(KP707106781), VADD(T17, T14));
		    {
			 Var T1f, TaB, T1i, TaC;
			 {
			      Var T1d, T1e, T1g, T1h;
			      T1d = LD(&(xi[WS(is, 122)]), ivs, &(xi[0]));
			      T1e = LD(&(xi[WS(is, 58)]), ivs, &(xi[0]));
			      T1f = VSUB(T1d, T1e);
			      TaB = VADD(T1d, T1e);
			      T1g = LD(&(xi[WS(is, 26)]), ivs, &(xi[0]));
			      T1h = LD(&(xi[WS(is, 90)]), ivs, &(xi[0]));
			      T1i = VSUB(T1g, T1h);
			      TaC = VADD(T1g, T1h);
			 }
			 T1j = VFNMS(LDK(KP923879532), T1i, VMUL(LDK(KP382683432), T1f));
			 Teo = VSUB(TaB, TaC);
			 T1z = VFMA(LDK(KP923879532), T1f, VMUL(LDK(KP382683432), T1i));
			 TaD = VADD(TaB, TaC);
		    }
		    {
			 Var T1m, TaE, T1p, TaF;
			 {
			      Var T1k, T1l, T1n, T1o;
			      T1k = LD(&(xi[WS(is, 10)]), ivs, &(xi[0]));
			      T1l = LD(&(xi[WS(is, 74)]), ivs, &(xi[0]));
			      T1m = VSUB(T1k, T1l);
			      TaE = VADD(T1k, T1l);
			      T1n = LD(&(xi[WS(is, 42)]), ivs, &(xi[0]));
			      T1o = LD(&(xi[WS(is, 106)]), ivs, &(xi[0]));
			      T1p = VSUB(T1n, T1o);
			      TaF = VADD(T1n, T1o);
			 }
			 T1q = VFMA(LDK(KP382683432), T1m, VMUL(LDK(KP923879532), T1p));
			 Ten = VSUB(TaE, TaF);
			 T1y = VFNMS(LDK(KP382683432), T1p, VMUL(LDK(KP923879532), T1m));
			 TaG = VADD(TaE, TaF);
		    }
		    Tdp = VADD(Taw, Taz);
		    Tdq = VADD(TaG, TaD);
		    TdG = VSUB(Tdp, Tdq);
		    Tem = VSUB(Tau, Tav);
		    Tep = VMUL(LDK(KP707106781), VADD(Ten, Teo));
		    Teq = VADD(Tem, Tep);
		    Tgm = VSUB(Tem, Tep);
		    {
			 Var Ter, Tes, T1c, T1r;
			 Ter = VSUB(Tay, Tax);
			 Tes = VMUL(LDK(KP707106781), VSUB(Teo, Ten));
			 Tet = VADD(Ter, Tes);
			 Tgl = VSUB(Tes, Ter);
			 T1c = VSUB(T18, T1b);
			 T1r = VSUB(T1j, T1q);
			 T1s = VADD(T1c, T1r);
			 T5P = VSUB(T1r, T1c);
		    }
		    {
			 Var T1x, T1A, T7b, T7c;
			 T1x = VADD(T1v, T1w);
			 T1A = VADD(T1y, T1z);
			 T1B = VADD(T1x, T1A);
			 T5Q = VSUB(T1x, T1A);
			 T7b = VADD(T1b, T18);
			 T7c = VSUB(T1z, T1y);
			 T7d = VADD(T7b, T7c);
			 T8Z = VSUB(T7c, T7b);
		    }
		    {
			 Var TaA, TaH, T78, T79;
			 TaA = VSUB(Taw, Taz);
			 TaH = VSUB(TaD, TaG);
			 TaI = VFMA(LDK(KP923879532), TaA, VMUL(LDK(KP382683432), TaH));
			 Tcf = VFNMS(LDK(KP382683432), TaA, VMUL(LDK(KP923879532), TaH));
			 T78 = VSUB(T1v, T1w);
			 T79 = VADD(T1q, T1j);
			 T7a = VADD(T78, T79);
			 T90 = VSUB(T78, T79);
		    }
	       }
	       {
		    Var TaJ, TaK, T1F, TaL, T27, TaM, TaN, T1M, TaO, T24, T1U, Tew, T2a, TaV, T21;
		    Var Tex, T29, TaS, T1I, T1L, Tev, Tey;
		    {
			 Var T1D, T1E, T25, T26;
			 T1D = LD(&(xi[WS(is, 126)]), ivs, &(xi[0]));
			 T1E = LD(&(xi[WS(is, 62)]), ivs, &(xi[0]));
			 TaJ = VADD(T1D, T1E);
			 T25 = LD(&(xi[WS(is, 30)]), ivs, &(xi[0]));
			 T26 = LD(&(xi[WS(is, 94)]), ivs, &(xi[0]));
			 TaK = VADD(T25, T26);
			 T1F = VSUB(T1D, T1E);
			 TaL = VADD(TaJ, TaK);
			 T27 = VSUB(T25, T26);
		    }
		    {
			 Var T1G, T1H, T1J, T1K;
			 T1G = LD(&(xi[WS(is, 14)]), ivs, &(xi[0]));
			 T1H = LD(&(xi[WS(is, 78)]), ivs, &(xi[0]));
			 T1I = VSUB(T1G, T1H);
			 TaM = VADD(T1G, T1H);
			 T1J = LD(&(xi[WS(is, 110)]), ivs, &(xi[0]));
			 T1K = LD(&(xi[WS(is, 46)]), ivs, &(xi[0]));
			 T1L = VSUB(T1J, T1K);
			 TaN = VADD(T1J, T1K);
		    }
		    T1M = VMUL(LDK(KP707106781), VADD(T1I, T1L));
		    TaO = VADD(TaM, TaN);
		    T24 = VMUL(LDK(KP707106781), VSUB(T1L, T1I));
		    {
			 Var T1Q, TaT, T1T, TaU;
			 {
			      Var T1O, T1P, T1R, T1S;
			      T1O = LD(&(xi[WS(is, 6)]), ivs, &(xi[0]));
			      T1P = LD(&(xi[WS(is, 70)]), ivs, &(xi[0]));
			      T1Q = VSUB(T1O, T1P);
			      TaT = VADD(T1O, T1P);
			      T1R = LD(&(xi[WS(is, 38)]), ivs, &(xi[0]));
			      T1S = LD(&(xi[WS(is, 102)]), ivs, &(xi[0]));
			      T1T = VSUB(T1R, T1S);
			      TaU = VADD(T1R, T1S);
			 }
			 T1U = VFNMS(LDK(KP382683432), T1T, VMUL(LDK(KP923879532), T1Q));
			 Tew = VSUB(TaT, TaU);
			 T2a = VFMA(LDK(KP382683432), T1Q, VMUL(LDK(KP923879532), T1T));
			 TaV = VADD(TaT, TaU);
		    }
		    {
			 Var T1X, TaQ, T20, TaR;
			 {
			      Var T1V, T1W, T1Y, T1Z;
			      T1V = LD(&(xi[WS(is, 118)]), ivs, &(xi[0]));
			      T1W = LD(&(xi[WS(is, 54)]), ivs, &(xi[0]));
			      T1X = VSUB(T1V, T1W);
			      TaQ = VADD(T1V, T1W);
			      T1Y = LD(&(xi[WS(is, 22)]), ivs, &(xi[0]));
			      T1Z = LD(&(xi[WS(is, 86)]), ivs, &(xi[0]));
			      T20 = VSUB(T1Y, T1Z);
			      TaR = VADD(T1Y, T1Z);
			 }
			 T21 = VFMA(LDK(KP923879532), T1X, VMUL(LDK(KP382683432), T20));
			 Tex = VSUB(TaQ, TaR);
			 T29 = VFNMS(LDK(KP923879532), T20, VMUL(LDK(KP382683432), T1X));
			 TaS = VADD(TaQ, TaR);
		    }
		    Tdm = VADD(TaL, TaO);
		    Tdn = VADD(TaV, TaS);
		    TdH = VSUB(Tdm, Tdn);
		    Tev = VSUB(TaJ, TaK);
		    Tey = VMUL(LDK(KP707106781), VADD(Tew, Tex));
		    Tez = VADD(Tev, Tey);
		    Tgi = VSUB(Tev, Tey);
		    {
			 Var TeA, TeB, T1N, T22;
			 TeA = VSUB(TaN, TaM);
			 TeB = VMUL(LDK(KP707106781), VSUB(Tex, Tew));
			 TeC = VADD(TeA, TeB);
			 Tgj = VSUB(TeB, TeA);
			 T1N = VADD(T1F, T1M);
			 T22 = VADD(T1U, T21);
			 T23 = VADD(T1N, T22);
			 T5N = VSUB(T1N, T22);
		    }
		    {
			 Var T28, T2b, T7i, T7j;
			 T28 = VSUB(T24, T27);
			 T2b = VSUB(T29, T2a);
			 T2c = VADD(T28, T2b);
			 T5M = VSUB(T2b, T28);
			 T7i = VADD(T27, T24);
			 T7j = VSUB(T21, T1U);
			 T7k = VADD(T7i, T7j);
			 T8X = VSUB(T7j, T7i);
		    }
		    {
			 Var TaP, TaW, T7f, T7g;
			 TaP = VSUB(TaL, TaO);
			 TaW = VSUB(TaS, TaV);
			 TaX = VFNMS(LDK(KP382683432), TaW, VMUL(LDK(KP923879532), TaP));
			 Tcg = VFMA(LDK(KP382683432), TaP, VMUL(LDK(KP923879532), TaW));
			 T7f = VSUB(T1F, T1M);
			 T7g = VADD(T2a, T29);
			 T7h = VADD(T7f, T7g);
			 T8W = VSUB(T7f, T7g);
		    }
	       }
	       {
		    Var T2J, TeL, T2V, Tb9, T30, TeO, T3c, Tbg, T2Q, TeM, T2S, Tbc, T37, TeP, T39;
		    Var Tbj;
		    {
			 Var T2H, T2I, Tb7, T2T, T2U, Tb8;
			 T2H = LD(&(xi[WS(is, 5)]), ivs, &(xi[WS(is, 1)]));
			 T2I = LD(&(xi[WS(is, 69)]), ivs, &(xi[WS(is, 1)]));
			 Tb7 = VADD(T2H, T2I);
			 T2T = LD(&(xi[WS(is, 37)]), ivs, &(xi[WS(is, 1)]));
			 T2U = LD(&(xi[WS(is, 101)]), ivs, &(xi[WS(is, 1)]));
			 Tb8 = VADD(T2T, T2U);
			 T2J = VSUB(T2H, T2I);
			 TeL = VSUB(Tb7, Tb8);
			 T2V = VSUB(T2T, T2U);
			 Tb9 = VADD(Tb7, Tb8);
		    }
		    {
			 Var T2Y, T2Z, Tbe, T3a, T3b, Tbf;
			 T2Y = LD(&(xi[WS(is, 125)]), ivs, &(xi[WS(is, 1)]));
			 T2Z = LD(&(xi[WS(is, 61)]), ivs, &(xi[WS(is, 1)]));
			 Tbe = VADD(T2Y, T2Z);
			 T3a = LD(&(xi[WS(is, 29)]), ivs, &(xi[WS(is, 1)]));
			 T3b = LD(&(xi[WS(is, 93)]), ivs, &(xi[WS(is, 1)]));
			 Tbf = VADD(T3a, T3b);
			 T30 = VSUB(T2Y, T2Z);
			 TeO = VSUB(Tbe, Tbf);
			 T3c = VSUB(T3a, T3b);
			 Tbg = VADD(Tbe, Tbf);
		    }
		    {
			 Var T2M, Tba, T2P, Tbb;
			 {
			      Var T2K, T2L, T2N, T2O;
			      T2K = LD(&(xi[WS(is, 21)]), ivs, &(xi[WS(is, 1)]));
			      T2L = LD(&(xi[WS(is, 85)]), ivs, &(xi[WS(is, 1)]));
			      T2M = VSUB(T2K, T2L);
			      Tba = VADD(T2K, T2L);
			      T2N = LD(&(xi[WS(is, 117)]), ivs, &(xi[WS(is, 1)]));
			      T2O = LD(&(xi[WS(is, 53)]), ivs, &(xi[WS(is, 1)]));
			      T2P = VSUB(T2N, T2O);
			      Tbb = VADD(T2N, T2O);
			 }
			 T2Q = VMUL(LDK(KP707106781), VADD(T2M, T2P));
			 TeM = VSUB(Tbb, Tba);
			 T2S = VMUL(LDK(KP707106781), VSUB(T2P, T2M));
			 Tbc = VADD(Tba, Tbb);
		    }
		    {
			 Var T33, Tbh, T36, Tbi;
			 {
			      Var T31, T32, T34, T35;
			      T31 = LD(&(xi[WS(is, 13)]), ivs, &(xi[WS(is, 1)]));
			      T32 = LD(&(xi[WS(is, 77)]), ivs, &(xi[WS(is, 1)]));
			      T33 = VSUB(T31, T32);
			      Tbh = VADD(T31, T32);
			      T34 = LD(&(xi[WS(is, 109)]), ivs, &(xi[WS(is, 1)]));
			      T35 = LD(&(xi[WS(is, 45)]), ivs, &(xi[WS(is, 1)]));
			      T36 = VSUB(T34, T35);
			      Tbi = VADD(T34, T35);
			 }
			 T37 = VMUL(LDK(KP707106781), VADD(T33, T36));
			 TeP = VSUB(Tbi, Tbh);
			 T39 = VMUL(LDK(KP707106781), VSUB(T36, T33));
			 Tbj = VADD(Tbh, Tbi);
		    }
		    {
			 Var Tbd, Tbk, TeN, TeQ;
			 Tbd = VSUB(Tb9, Tbc);
			 Tbk = VSUB(Tbg, Tbj);
			 Tbl = VMUL(LDK(KP707106781), VADD(Tbd, Tbk));
			 Tbu = VMUL(LDK(KP707106781), VSUB(Tbk, Tbd));
			 {
			      Var Td9, Tda, TeW, TeX;
			      Td9 = VADD(Tb9, Tbc);
			      Tda = VADD(Tbg, Tbj);
			      Tdb = VADD(Td9, Tda);
			      TdL = VSUB(Tda, Td9);
			      TeW = VFNMS(LDK(KP382683432), TeL, VMUL(LDK(KP923879532), TeM));
			      TeX = VFMA(LDK(KP382683432), TeO, VMUL(LDK(KP923879532), TeP));
			      TeY = VADD(TeW, TeX);
			      Tgu = VSUB(TeX, TeW);
			 }
			 TeN = VFMA(LDK(KP923879532), TeL, VMUL(LDK(KP382683432), TeM));
			 TeQ = VFNMS(LDK(KP382683432), TeP, VMUL(LDK(KP923879532), TeO));
			 TeR = VADD(TeN, TeQ);
			 Tgq = VSUB(TeQ, TeN);
			 {
			      Var T7t, T7C, T7w, T7D;
			      {
				   Var T7r, T7s, T7u, T7v;
				   T7r = VSUB(T2J, T2Q);
				   T7s = VADD(T2V, T2S);
				   T7t = VFMA(LDK(KP831469612), T7r, VMUL(LDK(KP555570233), T7s));
				   T7C = VFNMS(LDK(KP555570233), T7r, VMUL(LDK(KP831469612), T7s));
				   T7u = VSUB(T30, T37);
				   T7v = VADD(T3c, T39);
				   T7w = VFNMS(LDK(KP555570233), T7v, VMUL(LDK(KP831469612), T7u));
				   T7D = VFMA(LDK(KP555570233), T7u, VMUL(LDK(KP831469612), T7v));
			      }
			      T7x = VADD(T7t, T7w);
			      T98 = VSUB(T7D, T7C);
			      T7E = VADD(T7C, T7D);
			      T94 = VSUB(T7w, T7t);
			 }
			 {
			      Var T2X, T3q, T3e, T3r;
			      {
				   Var T2R, T2W, T38, T3d;
				   T2R = VADD(T2J, T2Q);
				   T2W = VSUB(T2S, T2V);
				   T2X = VFMA(LDK(KP980785280), T2R, VMUL(LDK(KP195090322), T2W));
				   T3q = VFNMS(LDK(KP195090322), T2R, VMUL(LDK(KP980785280), T2W));
				   T38 = VADD(T30, T37);
				   T3d = VSUB(T39, T3c);
				   T3e = VFNMS(LDK(KP195090322), T3d, VMUL(LDK(KP980785280), T38));
				   T3r = VFMA(LDK(KP195090322), T38, VMUL(LDK(KP980785280), T3d));
			      }
			      T3f = VADD(T2X, T3e);
			      T5Y = VSUB(T3r, T3q);
			      T3s = VADD(T3q, T3r);
			      T5U = VSUB(T3e, T2X);
			 }
		    }
	       }
	       {
		    Var T3Y, Tf6, T4a, TbG, T4f, Tf9, T4r, TbN, T45, Tf7, T47, TbJ, T4m, Tfa, T4o;
		    Var TbQ;
		    {
			 Var T3W, T3X, TbE, T48, T49, TbF;
			 T3W = LD(&(xi[WS(is, 3)]), ivs, &(xi[WS(is, 1)]));
			 T3X = LD(&(xi[WS(is, 67)]), ivs, &(xi[WS(is, 1)]));
			 TbE = VADD(T3W, T3X);
			 T48 = LD(&(xi[WS(is, 35)]), ivs, &(xi[WS(is, 1)]));
			 T49 = LD(&(xi[WS(is, 99)]), ivs, &(xi[WS(is, 1)]));
			 TbF = VADD(T48, T49);
			 T3Y = VSUB(T3W, T3X);
			 Tf6 = VSUB(TbE, TbF);
			 T4a = VSUB(T48, T49);
			 TbG = VADD(TbE, TbF);
		    }
		    {
			 Var T4d, T4e, TbL, T4p, T4q, TbM;
			 T4d = LD(&(xi[WS(is, 123)]), ivs, &(xi[WS(is, 1)]));
			 T4e = LD(&(xi[WS(is, 59)]), ivs, &(xi[WS(is, 1)]));
			 TbL = VADD(T4d, T4e);
			 T4p = LD(&(xi[WS(is, 27)]), ivs, &(xi[WS(is, 1)]));
			 T4q = LD(&(xi[WS(is, 91)]), ivs, &(xi[WS(is, 1)]));
			 TbM = VADD(T4p, T4q);
			 T4f = VSUB(T4d, T4e);
			 Tf9 = VSUB(TbL, TbM);
			 T4r = VSUB(T4p, T4q);
			 TbN = VADD(TbL, TbM);
		    }
		    {
			 Var T41, TbH, T44, TbI;
			 {
			      Var T3Z, T40, T42, T43;
			      T3Z = LD(&(xi[WS(is, 19)]), ivs, &(xi[WS(is, 1)]));
			      T40 = LD(&(xi[WS(is, 83)]), ivs, &(xi[WS(is, 1)]));
			      T41 = VSUB(T3Z, T40);
			      TbH = VADD(T3Z, T40);
			      T42 = LD(&(xi[WS(is, 115)]), ivs, &(xi[WS(is, 1)]));
			      T43 = LD(&(xi[WS(is, 51)]), ivs, &(xi[WS(is, 1)]));
			      T44 = VSUB(T42, T43);
			      TbI = VADD(T42, T43);
			 }
			 T45 = VMUL(LDK(KP707106781), VADD(T41, T44));
			 Tf7 = VSUB(TbI, TbH);
			 T47 = VMUL(LDK(KP707106781), VSUB(T44, T41));
			 TbJ = VADD(TbH, TbI);
		    }
		    {
			 Var T4i, TbO, T4l, TbP;
			 {
			      Var T4g, T4h, T4j, T4k;
			      T4g = LD(&(xi[WS(is, 11)]), ivs, &(xi[WS(is, 1)]));
			      T4h = LD(&(xi[WS(is, 75)]), ivs, &(xi[WS(is, 1)]));
			      T4i = VSUB(T4g, T4h);
			      TbO = VADD(T4g, T4h);
			      T4j = LD(&(xi[WS(is, 107)]), ivs, &(xi[WS(is, 1)]));
			      T4k = LD(&(xi[WS(is, 43)]), ivs, &(xi[WS(is, 1)]));
			      T4l = VSUB(T4j, T4k);
			      TbP = VADD(T4j, T4k);
			 }
			 T4m = VMUL(LDK(KP707106781), VADD(T4i, T4l));
			 Tfa = VSUB(TbP, TbO);
			 T4o = VMUL(LDK(KP707106781), VSUB(T4l, T4i));
			 TbQ = VADD(TbO, TbP);
		    }
		    {
			 Var TbK, TbR, Tf8, Tfb;
			 TbK = VSUB(TbG, TbJ);
			 TbR = VSUB(TbN, TbQ);
			 TbS = VMUL(LDK(KP707106781), VADD(TbK, TbR));
			 Tc1 = VMUL(LDK(KP707106781), VSUB(TbR, TbK));
			 {
			      Var Tdg, Tdh, Tfh, Tfi;
			      Tdg = VADD(TbG, TbJ);
			      Tdh = VADD(TbN, TbQ);
			      Tdi = VADD(Tdg, Tdh);
			      TdO = VSUB(Tdh, Tdg);
			      Tfh = VFNMS(LDK(KP382683432), Tf6, VMUL(LDK(KP923879532), Tf7));
			      Tfi = VFMA(LDK(KP382683432), Tf9, VMUL(LDK(KP923879532), Tfa));
			      Tfj = VADD(Tfh, Tfi);
			      Tgy = VSUB(Tfi, Tfh);
			 }
			 Tf8 = VFMA(LDK(KP923879532), Tf6, VMUL(LDK(KP382683432), Tf7));
			 Tfb = VFNMS(LDK(KP382683432), Tfa, VMUL(LDK(KP923879532), Tf9));
			 Tfc = VADD(Tf8, Tfb);
			 TgA = VSUB(Tfb, Tf8);
			 {
			      Var T7M, T7V, T7P, T7W;
			      {
				   Var T7K, T7L, T7N, T7O;
				   T7K = VSUB(T3Y, T45);
				   T7L = VADD(T4a, T47);
				   T7M = VFMA(LDK(KP831469612), T7K, VMUL(LDK(KP555570233), T7L));
				   T7V = VFNMS(LDK(KP555570233), T7K, VMUL(LDK(KP831469612), T7L));
				   T7N = VSUB(T4f, T4m);
				   T7O = VADD(T4r, T4o);
				   T7P = VFNMS(LDK(KP555570233), T7O, VMUL(LDK(KP831469612), T7N));
				   T7W = VFMA(LDK(KP555570233), T7N, VMUL(LDK(KP831469612), T7O));
			      }
			      T7Q = VADD(T7M, T7P);
			      T9e = VSUB(T7P, T7M);
			      T7X = VADD(T7V, T7W);
			      T9c = VSUB(T7W, T7V);
			 }
			 {
			      Var T4c, T4F, T4t, T4G;
			      {
				   Var T46, T4b, T4n, T4s;
				   T46 = VADD(T3Y, T45);
				   T4b = VSUB(T47, T4a);
				   T4c = VFMA(LDK(KP980785280), T46, VMUL(LDK(KP195090322), T4b));
				   T4F = VFNMS(LDK(KP195090322), T46, VMUL(LDK(KP980785280), T4b));
				   T4n = VADD(T4f, T4m);
				   T4s = VSUB(T4o, T4r);
				   T4t = VFNMS(LDK(KP195090322), T4s, VMUL(LDK(KP980785280), T4n));
				   T4G = VFMA(LDK(KP195090322), T4n, VMUL(LDK(KP980785280), T4s));
			      }
			      T4u = VADD(T4c, T4t);
			      T64 = VSUB(T4t, T4c);
			      T4H = VADD(T4F, T4G);
			      T62 = VSUB(T4G, T4F);
			 }
		    }
	       }
	       {
		    Var Td5, Tdx, TdC, TdE, Tdk, Tdt, Tds, Tdy, Tdz, TdD;
		    {
			 Var Td1, Td4, TdA, TdB;
			 Td1 = VADD(TcZ, Td0);
			 Td4 = VADD(Td2, Td3);
			 Td5 = VSUB(Td1, Td4);
			 Tdx = VADD(Td1, Td4);
			 TdA = VADD(Td8, Tdb);
			 TdB = VADD(Tdf, Tdi);
			 TdC = VADD(TdA, TdB);
			 TdE = VBYI(VSUB(TdB, TdA));
		    }
		    {
			 Var Tdc, Tdj, Tdo, Tdr;
			 Tdc = VSUB(Td8, Tdb);
			 Tdj = VSUB(Tdf, Tdi);
			 Tdk = VMUL(LDK(KP707106781), VADD(Tdc, Tdj));
			 Tdt = VMUL(LDK(KP707106781), VSUB(Tdj, Tdc));
			 Tdo = VADD(Tdm, Tdn);
			 Tdr = VADD(Tdp, Tdq);
			 Tds = VSUB(Tdo, Tdr);
			 Tdy = VADD(Tdr, Tdo);
		    }
		    Tdz = VADD(Tdx, Tdy);
		    ST(&(xo[WS(os, 64)]), VSUB(Tdz, TdC), ovs, &(xo[0]));
		    ST(&(xo[0]), VADD(Tdz, TdC), ovs, &(xo[0]));
		    TdD = VSUB(Tdx, Tdy);
		    ST(&(xo[WS(os, 96)]), VSUB(TdD, TdE), ovs, &(xo[0]));
		    ST(&(xo[WS(os, 32)]), VADD(TdD, TdE), ovs, &(xo[0]));
		    {
			 Var Tdl, Tdu, Tdv, Tdw;
			 Tdl = VADD(Td5, Tdk);
			 Tdu = VBYI(VADD(Tds, Tdt));
			 ST(&(xo[WS(os, 112)]), VSUB(Tdl, Tdu), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 16)]), VADD(Tdl, Tdu), ovs, &(xo[0]));
			 Tdv = VSUB(Td5, Tdk);
			 Tdw = VBYI(VSUB(Tdt, Tds));
			 ST(&(xo[WS(os, 80)]), VSUB(Tdv, Tdw), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 48)]), VADD(Tdv, Tdw), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var TdJ, Te4, TdX, Te5, TdQ, Te1, TdU, Te2;
		    {
			 Var TdF, TdI, TdV, TdW;
			 TdF = VSUB(TcZ, Td0);
			 TdI = VMUL(LDK(KP707106781), VADD(TdG, TdH));
			 TdJ = VADD(TdF, TdI);
			 Te4 = VSUB(TdF, TdI);
			 TdV = VFNMS(LDK(KP382683432), TdK, VMUL(LDK(KP923879532), TdL));
			 TdW = VFMA(LDK(KP382683432), TdN, VMUL(LDK(KP923879532), TdO));
			 TdX = VADD(TdV, TdW);
			 Te5 = VSUB(TdW, TdV);
		    }
		    {
			 Var TdM, TdP, TdS, TdT;
			 TdM = VFMA(LDK(KP923879532), TdK, VMUL(LDK(KP382683432), TdL));
			 TdP = VFNMS(LDK(KP382683432), TdO, VMUL(LDK(KP923879532), TdN));
			 TdQ = VADD(TdM, TdP);
			 Te1 = VSUB(TdP, TdM);
			 TdS = VSUB(Td3, Td2);
			 TdT = VMUL(LDK(KP707106781), VSUB(TdH, TdG));
			 TdU = VADD(TdS, TdT);
			 Te2 = VSUB(TdT, TdS);
		    }
		    {
			 Var TdR, TdY, Te7, Te8;
			 TdR = VADD(TdJ, TdQ);
			 TdY = VBYI(VADD(TdU, TdX));
			 ST(&(xo[WS(os, 120)]), VSUB(TdR, TdY), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 8)]), VADD(TdR, TdY), ovs, &(xo[0]));
			 Te7 = VBYI(VADD(Te2, Te1));
			 Te8 = VADD(Te4, Te5);
			 ST(&(xo[WS(os, 24)]), VADD(Te7, Te8), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 104)]), VSUB(Te8, Te7), ovs, &(xo[0]));
		    }
		    {
			 Var TdZ, Te0, Te3, Te6;
			 TdZ = VSUB(TdJ, TdQ);
			 Te0 = VBYI(VSUB(TdX, TdU));
			 ST(&(xo[WS(os, 72)]), VSUB(TdZ, Te0), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 56)]), VADD(TdZ, Te0), ovs, &(xo[0]));
			 Te3 = VBYI(VSUB(Te1, Te2));
			 Te6 = VSUB(Te4, Te5);
			 ST(&(xo[WS(os, 40)]), VADD(Te3, Te6), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 88)]), VSUB(Te6, Te3), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var TaZ, Tcs, Tci, Tcq, Tc4, Tct, Tcl, Tcp;
		    {
			 Var Tat, TaY, Tce, Tch;
			 Tat = VADD(Tad, Tas);
			 TaY = VADD(TaI, TaX);
			 TaZ = VADD(Tat, TaY);
			 Tcs = VSUB(Tat, TaY);
			 Tce = VADD(Tcc, Tcd);
			 Tch = VADD(Tcf, Tcg);
			 Tci = VADD(Tce, Tch);
			 Tcq = VSUB(Tch, Tce);
			 {
			      Var Tbw, Tcj, Tc3, Tck;
			      {
				   Var Tbm, Tbv, TbT, Tc2;
				   Tbm = VADD(Tb6, Tbl);
				   Tbv = VADD(Tbt, Tbu);
				   Tbw = VFMA(LDK(KP980785280), Tbm, VMUL(LDK(KP195090322), Tbv));
				   Tcj = VFNMS(LDK(KP195090322), Tbm, VMUL(LDK(KP980785280), Tbv));
				   TbT = VADD(TbD, TbS);
				   Tc2 = VADD(Tc0, Tc1);
				   Tc3 = VFNMS(LDK(KP195090322), Tc2, VMUL(LDK(KP980785280), TbT));
				   Tck = VFMA(LDK(KP195090322), TbT, VMUL(LDK(KP980785280), Tc2));
			      }
			      Tc4 = VADD(Tbw, Tc3);
			      Tct = VSUB(Tck, Tcj);
			      Tcl = VADD(Tcj, Tck);
			      Tcp = VSUB(Tc3, Tbw);
			 }
		    }
		    {
			 Var Tc5, Tcm, Tcv, Tcw;
			 Tc5 = VADD(TaZ, Tc4);
			 Tcm = VBYI(VADD(Tci, Tcl));
			 ST(&(xo[WS(os, 124)]), VSUB(Tc5, Tcm), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 4)]), VADD(Tc5, Tcm), ovs, &(xo[0]));
			 Tcv = VBYI(VADD(Tcq, Tcp));
			 Tcw = VADD(Tcs, Tct);
			 ST(&(xo[WS(os, 28)]), VADD(Tcv, Tcw), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 100)]), VSUB(Tcw, Tcv), ovs, &(xo[0]));
		    }
		    {
			 Var Tcn, Tco, Tcr, Tcu;
			 Tcn = VSUB(TaZ, Tc4);
			 Tco = VBYI(VSUB(Tcl, Tci));
			 ST(&(xo[WS(os, 68)]), VSUB(Tcn, Tco), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 60)]), VADD(Tcn, Tco), ovs, &(xo[0]));
			 Tcr = VBYI(VSUB(Tcp, Tcq));
			 Tcu = VSUB(Tcs, Tct);
			 ST(&(xo[WS(os, 36)]), VADD(Tcr, Tcu), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 92)]), VSUB(Tcu, Tcr), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var Tcz, TcU, TcK, TcS, TcG, TcV, TcN, TcR;
		    {
			 Var Tcx, Tcy, TcI, TcJ;
			 Tcx = VSUB(Tad, Tas);
			 Tcy = VSUB(Tcg, Tcf);
			 Tcz = VADD(Tcx, Tcy);
			 TcU = VSUB(Tcx, Tcy);
			 TcI = VSUB(Tcd, Tcc);
			 TcJ = VSUB(TaX, TaI);
			 TcK = VADD(TcI, TcJ);
			 TcS = VSUB(TcJ, TcI);
			 {
			      Var TcC, TcL, TcF, TcM;
			      {
				   Var TcA, TcB, TcD, TcE;
				   TcA = VSUB(Tb6, Tbl);
				   TcB = VSUB(Tbu, Tbt);
				   TcC = VFMA(LDK(KP831469612), TcA, VMUL(LDK(KP555570233), TcB));
				   TcL = VFNMS(LDK(KP555570233), TcA, VMUL(LDK(KP831469612), TcB));
				   TcD = VSUB(TbD, TbS);
				   TcE = VSUB(Tc1, Tc0);
				   TcF = VFNMS(LDK(KP555570233), TcE, VMUL(LDK(KP831469612), TcD));
				   TcM = VFMA(LDK(KP555570233), TcD, VMUL(LDK(KP831469612), TcE));
			      }
			      TcG = VADD(TcC, TcF);
			      TcV = VSUB(TcM, TcL);
			      TcN = VADD(TcL, TcM);
			      TcR = VSUB(TcF, TcC);
			 }
		    }
		    {
			 Var TcH, TcO, TcX, TcY;
			 TcH = VADD(Tcz, TcG);
			 TcO = VBYI(VADD(TcK, TcN));
			 ST(&(xo[WS(os, 116)]), VSUB(TcH, TcO), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 12)]), VADD(TcH, TcO), ovs, &(xo[0]));
			 TcX = VBYI(VADD(TcS, TcR));
			 TcY = VADD(TcU, TcV);
			 ST(&(xo[WS(os, 20)]), VADD(TcX, TcY), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 108)]), VSUB(TcY, TcX), ovs, &(xo[0]));
		    }
		    {
			 Var TcP, TcQ, TcT, TcW;
			 TcP = VSUB(Tcz, TcG);
			 TcQ = VBYI(VSUB(TcN, TcK));
			 ST(&(xo[WS(os, 76)]), VSUB(TcP, TcQ), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 52)]), VADD(TcP, TcQ), ovs, &(xo[0]));
			 TcT = VBYI(VSUB(TcR, TcS));
			 TcW = VSUB(TcU, TcV);
			 ST(&(xo[WS(os, 44)]), VADD(TcT, TcW), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 84)]), VSUB(TcW, TcT), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var TeF, Tg8, TfI, Tg0, Tfy, Tga, TfG, TfP, Tfm, TfJ, TfB, TfF, TfW, Tgb, Tg3;
		    Var Tg7;
		    {
			 Var Tel, TfY, TeE, TfZ, Teu, TeD;
			 Tel = VADD(Ted, Tek);
			 TfY = VSUB(Tft, Tfq);
			 Teu = VFMA(LDK(KP980785280), Teq, VMUL(LDK(KP195090322), Tet));
			 TeD = VFNMS(LDK(KP195090322), TeC, VMUL(LDK(KP980785280), Tez));
			 TeE = VADD(Teu, TeD);
			 TfZ = VSUB(TeD, Teu);
			 TeF = VADD(Tel, TeE);
			 Tg8 = VSUB(TfZ, TfY);
			 TfI = VSUB(Tel, TeE);
			 Tg0 = VADD(TfY, TfZ);
		    }
		    {
			 Var Tfu, TfN, Tfx, TfO, Tfv, Tfw;
			 Tfu = VADD(Tfq, Tft);
			 TfN = VSUB(Ted, Tek);
			 Tfv = VFNMS(LDK(KP195090322), Teq, VMUL(LDK(KP980785280), Tet));
			 Tfw = VFMA(LDK(KP195090322), Tez, VMUL(LDK(KP980785280), TeC));
			 Tfx = VADD(Tfv, Tfw);
			 TfO = VSUB(Tfw, Tfv);
			 Tfy = VADD(Tfu, Tfx);
			 Tga = VSUB(TfN, TfO);
			 TfG = VSUB(Tfx, Tfu);
			 TfP = VADD(TfN, TfO);
		    }
		    {
			 Var Tf0, Tfz, Tfl, TfA;
			 {
			      Var TeS, TeZ, Tfd, Tfk;
			      TeS = VADD(TeK, TeR);
			      TeZ = VADD(TeV, TeY);
			      Tf0 = VFMA(LDK(KP995184726), TeS, VMUL(LDK(KP098017140), TeZ));
			      Tfz = VFNMS(LDK(KP098017140), TeS, VMUL(LDK(KP995184726), TeZ));
			      Tfd = VADD(Tf5, Tfc);
			      Tfk = VADD(Tfg, Tfj);
			      Tfl = VFNMS(LDK(KP098017140), Tfk, VMUL(LDK(KP995184726), Tfd));
			      TfA = VFMA(LDK(KP098017140), Tfd, VMUL(LDK(KP995184726), Tfk));
			 }
			 Tfm = VADD(Tf0, Tfl);
			 TfJ = VSUB(TfA, Tfz);
			 TfB = VADD(Tfz, TfA);
			 TfF = VSUB(Tfl, Tf0);
		    }
		    {
			 Var TfS, Tg1, TfV, Tg2;
			 {
			      Var TfQ, TfR, TfT, TfU;
			      TfQ = VSUB(TeK, TeR);
			      TfR = VSUB(TeY, TeV);
			      TfS = VFMA(LDK(KP773010453), TfQ, VMUL(LDK(KP634393284), TfR));
			      Tg1 = VFNMS(LDK(KP634393284), TfQ, VMUL(LDK(KP773010453), TfR));
			      TfT = VSUB(Tf5, Tfc);
			      TfU = VSUB(Tfj, Tfg);
			      TfV = VFNMS(LDK(KP634393284), TfU, VMUL(LDK(KP773010453), TfT));
			      Tg2 = VFMA(LDK(KP634393284), TfT, VMUL(LDK(KP773010453), TfU));
			 }
			 TfW = VADD(TfS, TfV);
			 Tgb = VSUB(Tg2, Tg1);
			 Tg3 = VADD(Tg1, Tg2);
			 Tg7 = VSUB(TfV, TfS);
		    }
		    {
			 Var Tfn, TfC, Tg9, Tgc;
			 Tfn = VADD(TeF, Tfm);
			 TfC = VBYI(VADD(Tfy, TfB));
			 ST(&(xo[WS(os, 126)]), VSUB(Tfn, TfC), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 2)]), VADD(Tfn, TfC), ovs, &(xo[0]));
			 Tg9 = VBYI(VSUB(Tg7, Tg8));
			 Tgc = VSUB(Tga, Tgb);
			 ST(&(xo[WS(os, 46)]), VADD(Tg9, Tgc), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 82)]), VSUB(Tgc, Tg9), ovs, &(xo[0]));
		    }
		    {
			 Var Tgd, Tge, TfD, TfE;
			 Tgd = VBYI(VADD(Tg8, Tg7));
			 Tge = VADD(Tga, Tgb);
			 ST(&(xo[WS(os, 18)]), VADD(Tgd, Tge), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 110)]), VSUB(Tge, Tgd), ovs, &(xo[0]));
			 TfD = VSUB(TeF, Tfm);
			 TfE = VBYI(VSUB(TfB, Tfy));
			 ST(&(xo[WS(os, 66)]), VSUB(TfD, TfE), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 62)]), VADD(TfD, TfE), ovs, &(xo[0]));
		    }
		    {
			 Var TfH, TfK, TfX, Tg4;
			 TfH = VBYI(VSUB(TfF, TfG));
			 TfK = VSUB(TfI, TfJ);
			 ST(&(xo[WS(os, 34)]), VADD(TfH, TfK), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 94)]), VSUB(TfK, TfH), ovs, &(xo[0]));
			 TfX = VADD(TfP, TfW);
			 Tg4 = VBYI(VADD(Tg0, Tg3));
			 ST(&(xo[WS(os, 114)]), VSUB(TfX, Tg4), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 14)]), VADD(TfX, Tg4), ovs, &(xo[0]));
		    }
		    {
			 Var Tg5, Tg6, TfL, TfM;
			 Tg5 = VSUB(TfP, TfW);
			 Tg6 = VBYI(VSUB(Tg3, Tg0));
			 ST(&(xo[WS(os, 78)]), VSUB(Tg5, Tg6), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 50)]), VADD(Tg5, Tg6), ovs, &(xo[0]));
			 TfL = VBYI(VADD(TfG, TfF));
			 TfM = VADD(TfI, TfJ);
			 ST(&(xo[WS(os, 30)]), VADD(TfL, TfM), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 98)]), VSUB(TfM, TfL), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var Tgp, Thm, TgW, The, TgM, Tho, TgU, Th3, TgE, TgX, TgP, TgT, Tha, Thp, Thh;
		    Var Thl;
		    {
			 Var Tgh, Thc, Tgo, Thd, Tgk, Tgn;
			 Tgh = VSUB(Tgf, Tgg);
			 Thc = VADD(TgH, TgG);
			 Tgk = VFMA(LDK(KP555570233), Tgi, VMUL(LDK(KP831469612), Tgj));
			 Tgn = VFNMS(LDK(KP555570233), Tgm, VMUL(LDK(KP831469612), Tgl));
			 Tgo = VSUB(Tgk, Tgn);
			 Thd = VADD(Tgn, Tgk);
			 Tgp = VADD(Tgh, Tgo);
			 Thm = VSUB(Thd, Thc);
			 TgW = VSUB(Tgh, Tgo);
			 The = VADD(Thc, Thd);
		    }
		    {
			 Var TgI, Th1, TgL, Th2, TgJ, TgK;
			 TgI = VSUB(TgG, TgH);
			 Th1 = VADD(Tgf, Tgg);
			 TgJ = VFNMS(LDK(KP555570233), Tgj, VMUL(LDK(KP831469612), Tgi));
			 TgK = VFMA(LDK(KP831469612), Tgm, VMUL(LDK(KP555570233), Tgl));
			 TgL = VSUB(TgJ, TgK);
			 Th2 = VADD(TgK, TgJ);
			 TgM = VADD(TgI, TgL);
			 Tho = VSUB(Th1, Th2);
			 TgU = VSUB(TgL, TgI);
			 Th3 = VADD(Th1, Th2);
		    }
		    {
			 Var Tgw, TgN, TgD, TgO;
			 {
			      Var Tgs, Tgv, Tgz, TgC;
			      Tgs = VSUB(Tgq, Tgr);
			      Tgv = VSUB(Tgt, Tgu);
			      Tgw = VFMA(LDK(KP471396736), Tgs, VMUL(LDK(KP881921264), Tgv));
			      TgN = VFNMS(LDK(KP471396736), Tgv, VMUL(LDK(KP881921264), Tgs));
			      Tgz = VSUB(Tgx, Tgy);
			      TgC = VSUB(TgA, TgB);
			      TgD = VFNMS(LDK(KP471396736), TgC, VMUL(LDK(KP881921264), Tgz));
			      TgO = VFMA(LDK(KP881921264), TgC, VMUL(LDK(KP471396736), Tgz));
			 }
			 TgE = VADD(Tgw, TgD);
			 TgX = VSUB(TgO, TgN);
			 TgP = VADD(TgN, TgO);
			 TgT = VSUB(TgD, Tgw);
		    }
		    {
			 Var Th6, Thf, Th9, Thg;
			 {
			      Var Th4, Th5, Th7, Th8;
			      Th4 = VADD(Tgr, Tgq);
			      Th5 = VADD(Tgt, Tgu);
			      Th6 = VFMA(LDK(KP290284677), Th4, VMUL(LDK(KP956940335), Th5));
			      Thf = VFNMS(LDK(KP290284677), Th5, VMUL(LDK(KP956940335), Th4));
			      Th7 = VADD(Tgx, Tgy);
			      Th8 = VADD(TgB, TgA);
			      Th9 = VFNMS(LDK(KP290284677), Th8, VMUL(LDK(KP956940335), Th7));
			      Thg = VFMA(LDK(KP956940335), Th8, VMUL(LDK(KP290284677), Th7));
			 }
			 Tha = VADD(Th6, Th9);
			 Thp = VSUB(Thg, Thf);
			 Thh = VADD(Thf, Thg);
			 Thl = VSUB(Th9, Th6);
		    }
		    {
			 Var TgF, TgQ, Thn, Thq;
			 TgF = VADD(Tgp, TgE);
			 TgQ = VBYI(VADD(TgM, TgP));
			 ST(&(xo[WS(os, 118)]), VSUB(TgF, TgQ), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 10)]), VADD(TgF, TgQ), ovs, &(xo[0]));
			 Thn = VBYI(VSUB(Thl, Thm));
			 Thq = VSUB(Tho, Thp);
			 ST(&(xo[WS(os, 38)]), VADD(Thn, Thq), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 90)]), VSUB(Thq, Thn), ovs, &(xo[0]));
		    }
		    {
			 Var Thr, Ths, TgR, TgS;
			 Thr = VBYI(VADD(Thm, Thl));
			 Ths = VADD(Tho, Thp);
			 ST(&(xo[WS(os, 26)]), VADD(Thr, Ths), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 102)]), VSUB(Ths, Thr), ovs, &(xo[0]));
			 TgR = VSUB(Tgp, TgE);
			 TgS = VBYI(VSUB(TgP, TgM));
			 ST(&(xo[WS(os, 74)]), VSUB(TgR, TgS), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 54)]), VADD(TgR, TgS), ovs, &(xo[0]));
		    }
		    {
			 Var TgV, TgY, Thb, Thi;
			 TgV = VBYI(VSUB(TgT, TgU));
			 TgY = VSUB(TgW, TgX);
			 ST(&(xo[WS(os, 42)]), VADD(TgV, TgY), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 86)]), VSUB(TgY, TgV), ovs, &(xo[0]));
			 Thb = VADD(Th3, Tha);
			 Thi = VBYI(VADD(The, Thh));
			 ST(&(xo[WS(os, 122)]), VSUB(Thb, Thi), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 6)]), VADD(Thb, Thi), ovs, &(xo[0]));
		    }
		    {
			 Var Thj, Thk, TgZ, Th0;
			 Thj = VSUB(Th3, Tha);
			 Thk = VBYI(VSUB(Thh, The));
			 ST(&(xo[WS(os, 70)]), VSUB(Thj, Thk), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 58)]), VADD(Thj, Thk), ovs, &(xo[0]));
			 TgZ = VBYI(VADD(TgU, TgT));
			 Th0 = VADD(TgW, TgX);
			 ST(&(xo[WS(os, 22)]), VADD(TgZ, Th0), ovs, &(xo[0]));
			 ST(&(xo[WS(os, 106)]), VSUB(Th0, TgZ), ovs, &(xo[0]));
		    }
	       }
	       {
		    Var T80, T8n, T8f, T8j, T8A, T8P, T8H, T8L, T7n, T8M, T8O, T8c, T8k, T8t, T8E;
		    Var T8m;
		    {
			 Var T7G, T8d, T7Z, T8e;
			 {
			      Var T7y, T7F, T7R, T7Y;
			      T7y = VADD(T7q, T7x);
			      T7F = VADD(T7B, T7E);
			      T7G = VFMA(LDK(KP989176509), T7y, VMUL(LDK(KP146730474), T7F));
			      T8d = VFNMS(LDK(KP146730474), T7y, VMUL(LDK(KP989176509), T7F));
			      T7R = VADD(T7J, T7Q);
			      T7Y = VADD(T7U, T7X);
			      T7Z = VFNMS(LDK(KP146730474), T7Y, VMUL(LDK(KP989176509), T7R));
			      T8e = VFMA(LDK(KP146730474), T7R, VMUL(LDK(KP989176509), T7Y));
			 }
			 T80 = VADD(T7G, T7Z);
			 T8n = VSUB(T8e, T8d);
			 T8f = VADD(T8d, T8e);
			 T8j = VSUB(T7Z, T7G);
		    }
		    {
			 Var T8w, T8F, T8z, T8G;
			 {
			      Var T8u, T8v, T8x, T8y;
			      T8u = VSUB(T7q, T7x);
			      T8v = VSUB(T7E, T7B);
			      T8w = VFMA(LDK(KP803207531), T8u, VMUL(LDK(KP595699304), T8v));
			      T8F = VFNMS(LDK(KP595699304), T8u, VMUL(LDK(KP803207531), T8v));
			      T8x = VSUB(T7J, T7Q);
			      T8y = VSUB(T7X, T7U);
			      T8z = VFNMS(LDK(KP595699304), T8y, VMUL(LDK(KP803207531), T8x));
			      T8G = VFMA(LDK(KP595699304), T8x, VMUL(LDK(KP803207531), T8y));
			 }
			 T8A = VADD(T8w, T8z);
			 T8P = VSUB(T8G, T8F);
			 T8H = VADD(T8F, T8G);
			 T8L = VSUB(T8z, T8w);
		    }
		    {
			 Var T77, T8r, T88, T8C, T7m, T8D, T8b, T8s, T76, T87;
			 T76 = VADD(T72, T75);
			 T77 = VADD(T6Z, T76);
			 T8r = VSUB(T6Z, T76);
			 T87 = VADD(T85, T86);
			 T88 = VADD(T84, T87);
			 T8C = VSUB(T87, T84);
			 {
			      Var T7e, T7l, T89, T8a;
			      T7e = VFMA(LDK(KP956940335), T7a, VMUL(LDK(KP290284677), T7d));
			      T7l = VFNMS(LDK(KP290284677), T7k, VMUL(LDK(KP956940335), T7h));
			      T7m = VADD(T7e, T7l);
			      T8D = VSUB(T7l, T7e);
			      T89 = VFNMS(LDK(KP290284677), T7a, VMUL(LDK(KP956940335), T7d));
			      T8a = VFMA(LDK(KP290284677), T7h, VMUL(LDK(KP956940335), T7k));
			      T8b = VADD(T89, T8a);
			      T8s = VSUB(T8a, T89);
			 }
			 T7n = VADD(T77, T7m);
			 T8M = VSUB(T8D, T8C);
			 T8O = VSUB(T8r, T8s);
			 T8c = VADD(T88, T8b);
			 T8k = VSUB(T8b, T88);
			 T8t = VADD(T8r, T8s);
			 T8E = VADD(T8C, T8D);
			 T8m = VSUB(T77, T7m);
		    }
		    {
			 Var T81, T8g, T8N, T8Q;
			 T81 = VADD(T7n, T80);
			 T8g = VBYI(VADD(T8c, T8f));
			 ST(&(xo[WS(os, 125)]), VSUB(T81, T8g), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 3)]), VADD(T81, T8g), ovs, &(xo[WS(os, 1)]));
			 T8N = VBYI(VSUB(T8L, T8M));
			 T8Q = VSUB(T8O, T8P);
			 ST(&(xo[WS(os, 45)]), VADD(T8N, T8Q), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 83)]), VSUB(T8Q, T8N), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T8R, T8S, T8h, T8i;
			 T8R = VBYI(VADD(T8M, T8L));
			 T8S = VADD(T8O, T8P);
			 ST(&(xo[WS(os, 19)]), VADD(T8R, T8S), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 109)]), VSUB(T8S, T8R), ovs, &(xo[WS(os, 1)]));
			 T8h = VSUB(T7n, T80);
			 T8i = VBYI(VSUB(T8f, T8c));
			 ST(&(xo[WS(os, 67)]), VSUB(T8h, T8i), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 61)]), VADD(T8h, T8i), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T8l, T8o, T8B, T8I;
			 T8l = VBYI(VSUB(T8j, T8k));
			 T8o = VSUB(T8m, T8n);
			 ST(&(xo[WS(os, 35)]), VADD(T8l, T8o), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 93)]), VSUB(T8o, T8l), ovs, &(xo[WS(os, 1)]));
			 T8B = VADD(T8t, T8A);
			 T8I = VBYI(VADD(T8E, T8H));
			 ST(&(xo[WS(os, 115)]), VSUB(T8B, T8I), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 13)]), VADD(T8B, T8I), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T8J, T8K, T8p, T8q;
			 T8J = VSUB(T8t, T8A);
			 T8K = VBYI(VSUB(T8H, T8E));
			 ST(&(xo[WS(os, 77)]), VSUB(T8J, T8K), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 51)]), VADD(T8J, T8K), ovs, &(xo[WS(os, 1)]));
			 T8p = VBYI(VADD(T8k, T8j));
			 T8q = VADD(T8m, T8n);
			 ST(&(xo[WS(os, 29)]), VADD(T8p, T8q), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 99)]), VSUB(T8q, T8p), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T4K, T5d, T55, T59, T5q, T5F, T5x, T5B, T2f, T5C, T5E, T52, T5a, T5j, T5u;
		    Var T5c;
		    {
			 Var T3u, T53, T4J, T54;
			 {
			      Var T3g, T3t, T4v, T4I;
			      T3g = VADD(T2G, T3f);
			      T3t = VADD(T3p, T3s);
			      T3u = VFMA(LDK(KP998795456), T3g, VMUL(LDK(KP049067674), T3t));
			      T53 = VFNMS(LDK(KP049067674), T3g, VMUL(LDK(KP998795456), T3t));
			      T4v = VADD(T3V, T4u);
			      T4I = VADD(T4E, T4H);
			      T4J = VFNMS(LDK(KP049067674), T4I, VMUL(LDK(KP998795456), T4v));
			      T54 = VFMA(LDK(KP049067674), T4v, VMUL(LDK(KP998795456), T4I));
			 }
			 T4K = VADD(T3u, T4J);
			 T5d = VSUB(T54, T53);
			 T55 = VADD(T53, T54);
			 T59 = VSUB(T4J, T3u);
		    }
		    {
			 Var T5m, T5v, T5p, T5w;
			 {
			      Var T5k, T5l, T5n, T5o;
			      T5k = VSUB(T2G, T3f);
			      T5l = VSUB(T3s, T3p);
			      T5m = VFMA(LDK(KP740951125), T5k, VMUL(LDK(KP671558954), T5l));
			      T5v = VFNMS(LDK(KP671558954), T5k, VMUL(LDK(KP740951125), T5l));
			      T5n = VSUB(T3V, T4u);
			      T5o = VSUB(T4H, T4E);
			      T5p = VFNMS(LDK(KP671558954), T5o, VMUL(LDK(KP740951125), T5n));
			      T5w = VFMA(LDK(KP671558954), T5n, VMUL(LDK(KP740951125), T5o));
			 }
			 T5q = VADD(T5m, T5p);
			 T5F = VSUB(T5w, T5v);
			 T5x = VADD(T5v, T5w);
			 T5B = VSUB(T5p, T5m);
		    }
		    {
			 Var T11, T5h, T4Y, T5s, T2e, T5t, T51, T5i, T10, T4X;
			 T10 = VADD(TI, TZ);
			 T11 = VADD(Tr, T10);
			 T5h = VSUB(Tr, T10);
			 T4X = VADD(T4V, T4W);
			 T4Y = VADD(T4U, T4X);
			 T5s = VSUB(T4X, T4U);
			 {
			      Var T1C, T2d, T4Z, T50;
			      T1C = VFMA(LDK(KP098017140), T1s, VMUL(LDK(KP995184726), T1B));
			      T2d = VFNMS(LDK(KP098017140), T2c, VMUL(LDK(KP995184726), T23));
			      T2e = VADD(T1C, T2d);
			      T5t = VSUB(T2d, T1C);
			      T4Z = VFNMS(LDK(KP098017140), T1B, VMUL(LDK(KP995184726), T1s));
			      T50 = VFMA(LDK(KP995184726), T2c, VMUL(LDK(KP098017140), T23));
			      T51 = VADD(T4Z, T50);
			      T5i = VSUB(T50, T4Z);
			 }
			 T2f = VADD(T11, T2e);
			 T5C = VSUB(T5t, T5s);
			 T5E = VSUB(T5h, T5i);
			 T52 = VADD(T4Y, T51);
			 T5a = VSUB(T51, T4Y);
			 T5j = VADD(T5h, T5i);
			 T5u = VADD(T5s, T5t);
			 T5c = VSUB(T11, T2e);
		    }
		    {
			 Var T4L, T56, T5D, T5G;
			 T4L = VADD(T2f, T4K);
			 T56 = VBYI(VADD(T52, T55));
			 ST(&(xo[WS(os, 127)]), VSUB(T4L, T56), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 1)]), VADD(T4L, T56), ovs, &(xo[WS(os, 1)]));
			 T5D = VBYI(VSUB(T5B, T5C));
			 T5G = VSUB(T5E, T5F);
			 ST(&(xo[WS(os, 47)]), VADD(T5D, T5G), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 81)]), VSUB(T5G, T5D), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T5H, T5I, T57, T58;
			 T5H = VBYI(VADD(T5C, T5B));
			 T5I = VADD(T5E, T5F);
			 ST(&(xo[WS(os, 17)]), VADD(T5H, T5I), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 111)]), VSUB(T5I, T5H), ovs, &(xo[WS(os, 1)]));
			 T57 = VSUB(T2f, T4K);
			 T58 = VBYI(VSUB(T55, T52));
			 ST(&(xo[WS(os, 65)]), VSUB(T57, T58), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 63)]), VADD(T57, T58), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T5b, T5e, T5r, T5y;
			 T5b = VBYI(VSUB(T59, T5a));
			 T5e = VSUB(T5c, T5d);
			 ST(&(xo[WS(os, 33)]), VADD(T5b, T5e), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 95)]), VSUB(T5e, T5b), ovs, &(xo[WS(os, 1)]));
			 T5r = VADD(T5j, T5q);
			 T5y = VBYI(VADD(T5u, T5x));
			 ST(&(xo[WS(os, 113)]), VSUB(T5r, T5y), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 15)]), VADD(T5r, T5y), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T5z, T5A, T5f, T5g;
			 T5z = VSUB(T5j, T5q);
			 T5A = VBYI(VSUB(T5x, T5u));
			 ST(&(xo[WS(os, 79)]), VSUB(T5z, T5A), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 49)]), VADD(T5z, T5A), ovs, &(xo[WS(os, 1)]));
			 T5f = VBYI(VADD(T5a, T59));
			 T5g = VADD(T5c, T5d);
			 ST(&(xo[WS(os, 31)]), VADD(T5f, T5g), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 97)]), VSUB(T5g, T5f), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T9i, T9B, T9t, T9x, T9O, Ta3, T9V, T9Z, T93, Ta0, Ta2, T9q, T9y, T9H, T9S;
		    Var T9A;
		    {
			 Var T9a, T9r, T9h, T9s;
			 {
			      Var T96, T99, T9d, T9g;
			      T96 = VSUB(T94, T95);
			      T99 = VSUB(T97, T98);
			      T9a = VFMA(LDK(KP514102744), T96, VMUL(LDK(KP857728610), T99));
			      T9r = VFNMS(LDK(KP514102744), T99, VMUL(LDK(KP857728610), T96));
			      T9d = VSUB(T9b, T9c);
			      T9g = VSUB(T9e, T9f);
			      T9h = VFNMS(LDK(KP514102744), T9g, VMUL(LDK(KP857728610), T9d));
			      T9s = VFMA(LDK(KP857728610), T9g, VMUL(LDK(KP514102744), T9d));
			 }
			 T9i = VADD(T9a, T9h);
			 T9B = VSUB(T9s, T9r);
			 T9t = VADD(T9r, T9s);
			 T9x = VSUB(T9h, T9a);
		    }
		    {
			 Var T9K, T9T, T9N, T9U;
			 {
			      Var T9I, T9J, T9L, T9M;
			      T9I = VADD(T95, T94);
			      T9J = VADD(T97, T98);
			      T9K = VFMA(LDK(KP242980179), T9I, VMUL(LDK(KP970031253), T9J));
			      T9T = VFNMS(LDK(KP242980179), T9J, VMUL(LDK(KP970031253), T9I));
			      T9L = VADD(T9b, T9c);
			      T9M = VADD(T9f, T9e);
			      T9N = VFNMS(LDK(KP242980179), T9M, VMUL(LDK(KP970031253), T9L));
			      T9U = VFMA(LDK(KP970031253), T9M, VMUL(LDK(KP242980179), T9L));
			 }
			 T9O = VADD(T9K, T9N);
			 Ta3 = VSUB(T9U, T9T);
			 T9V = VADD(T9T, T9U);
			 T9Z = VSUB(T9N, T9K);
		    }
		    {
			 Var T8V, T9F, T9m, T9Q, T92, T9R, T9p, T9G, T8U, T9k;
			 T8U = VSUB(T86, T85);
			 T8V = VSUB(T8T, T8U);
			 T9F = VADD(T8T, T8U);
			 T9k = VSUB(T75, T72);
			 T9m = VSUB(T9k, T9l);
			 T9Q = VADD(T9l, T9k);
			 {
			      Var T8Y, T91, T9n, T9o;
			      T8Y = VFMA(LDK(KP471396736), T8W, VMUL(LDK(KP881921264), T8X));
			      T91 = VFNMS(LDK(KP471396736), T90, VMUL(LDK(KP881921264), T8Z));
			      T92 = VSUB(T8Y, T91);
			      T9R = VADD(T91, T8Y);
			      T9n = VFNMS(LDK(KP471396736), T8X, VMUL(LDK(KP881921264), T8W));
			      T9o = VFMA(LDK(KP881921264), T90, VMUL(LDK(KP471396736), T8Z));
			      T9p = VSUB(T9n, T9o);
			      T9G = VADD(T9o, T9n);
			 }
			 T93 = VADD(T8V, T92);
			 Ta0 = VSUB(T9R, T9Q);
			 Ta2 = VSUB(T9F, T9G);
			 T9q = VADD(T9m, T9p);
			 T9y = VSUB(T9p, T9m);
			 T9H = VADD(T9F, T9G);
			 T9S = VADD(T9Q, T9R);
			 T9A = VSUB(T8V, T92);
		    }
		    {
			 Var T9j, T9u, Ta1, Ta4;
			 T9j = VADD(T93, T9i);
			 T9u = VBYI(VADD(T9q, T9t));
			 ST(&(xo[WS(os, 117)]), VSUB(T9j, T9u), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 11)]), VADD(T9j, T9u), ovs, &(xo[WS(os, 1)]));
			 Ta1 = VBYI(VSUB(T9Z, Ta0));
			 Ta4 = VSUB(Ta2, Ta3);
			 ST(&(xo[WS(os, 37)]), VADD(Ta1, Ta4), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 91)]), VSUB(Ta4, Ta1), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var Ta5, Ta6, T9v, T9w;
			 Ta5 = VBYI(VADD(Ta0, T9Z));
			 Ta6 = VADD(Ta2, Ta3);
			 ST(&(xo[WS(os, 27)]), VADD(Ta5, Ta6), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 101)]), VSUB(Ta6, Ta5), ovs, &(xo[WS(os, 1)]));
			 T9v = VSUB(T93, T9i);
			 T9w = VBYI(VSUB(T9t, T9q));
			 ST(&(xo[WS(os, 75)]), VSUB(T9v, T9w), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 53)]), VADD(T9v, T9w), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T9z, T9C, T9P, T9W;
			 T9z = VBYI(VSUB(T9x, T9y));
			 T9C = VSUB(T9A, T9B);
			 ST(&(xo[WS(os, 43)]), VADD(T9z, T9C), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 85)]), VSUB(T9C, T9z), ovs, &(xo[WS(os, 1)]));
			 T9P = VADD(T9H, T9O);
			 T9W = VBYI(VADD(T9S, T9V));
			 ST(&(xo[WS(os, 123)]), VSUB(T9P, T9W), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 5)]), VADD(T9P, T9W), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T9X, T9Y, T9D, T9E;
			 T9X = VSUB(T9H, T9O);
			 T9Y = VBYI(VSUB(T9V, T9S));
			 ST(&(xo[WS(os, 69)]), VSUB(T9X, T9Y), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 59)]), VADD(T9X, T9Y), ovs, &(xo[WS(os, 1)]));
			 T9D = VBYI(VADD(T9y, T9x));
			 T9E = VADD(T9A, T9B);
			 ST(&(xo[WS(os, 21)]), VADD(T9D, T9E), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 107)]), VSUB(T9E, T9D), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	       {
		    Var T68, T6r, T6j, T6n, T6E, T6T, T6L, T6P, T5T, T6Q, T6S, T6g, T6o, T6x, T6I;
		    Var T6q;
		    {
			 Var T60, T6h, T67, T6i;
			 {
			      Var T5W, T5Z, T63, T66;
			      T5W = VSUB(T5U, T5V);
			      T5Z = VSUB(T5X, T5Y);
			      T60 = VFMA(LDK(KP427555093), T5W, VMUL(LDK(KP903989293), T5Z));
			      T6h = VFNMS(LDK(KP427555093), T5Z, VMUL(LDK(KP903989293), T5W));
			      T63 = VSUB(T61, T62);
			      T66 = VSUB(T64, T65);
			      T67 = VFNMS(LDK(KP427555093), T66, VMUL(LDK(KP903989293), T63));
			      T6i = VFMA(LDK(KP903989293), T66, VMUL(LDK(KP427555093), T63));
			 }
			 T68 = VADD(T60, T67);
			 T6r = VSUB(T6i, T6h);
			 T6j = VADD(T6h, T6i);
			 T6n = VSUB(T67, T60);
		    }
		    {
			 Var T6A, T6J, T6D, T6K;
			 {
			      Var T6y, T6z, T6B, T6C;
			      T6y = VADD(T5V, T5U);
			      T6z = VADD(T5X, T5Y);
			      T6A = VFMA(LDK(KP336889853), T6y, VMUL(LDK(KP941544065), T6z));
			      T6J = VFNMS(LDK(KP336889853), T6z, VMUL(LDK(KP941544065), T6y));
			      T6B = VADD(T61, T62);
			      T6C = VADD(T65, T64);
			      T6D = VFNMS(LDK(KP336889853), T6C, VMUL(LDK(KP941544065), T6B));
			      T6K = VFMA(LDK(KP941544065), T6C, VMUL(LDK(KP336889853), T6B));
			 }
			 T6E = VADD(T6A, T6D);
			 T6T = VSUB(T6K, T6J);
			 T6L = VADD(T6J, T6K);
			 T6P = VSUB(T6D, T6A);
		    }
		    {
			 Var T5L, T6v, T6c, T6G, T5S, T6H, T6f, T6w, T5K, T6a;
			 T5K = VSUB(T4W, T4V);
			 T5L = VSUB(T5J, T5K);
			 T6v = VADD(T5J, T5K);
			 T6a = VSUB(TZ, TI);
			 T6c = VSUB(T6a, T6b);
			 T6G = VADD(T6b, T6a);
			 {
			      Var T5O, T5R, T6d, T6e;
			      T5O = VFMA(LDK(KP773010453), T5M, VMUL(LDK(KP634393284), T5N));
			      T5R = VFNMS(LDK(KP634393284), T5Q, VMUL(LDK(KP773010453), T5P));
			      T5S = VSUB(T5O, T5R);
			      T6H = VADD(T5R, T5O);
			      T6d = VFNMS(LDK(KP634393284), T5M, VMUL(LDK(KP773010453), T5N));
			      T6e = VFMA(LDK(KP634393284), T5P, VMUL(LDK(KP773010453), T5Q));
			      T6f = VSUB(T6d, T6e);
			      T6w = VADD(T6e, T6d);
			 }
			 T5T = VADD(T5L, T5S);
			 T6Q = VSUB(T6H, T6G);
			 T6S = VSUB(T6v, T6w);
			 T6g = VADD(T6c, T6f);
			 T6o = VSUB(T6f, T6c);
			 T6x = VADD(T6v, T6w);
			 T6I = VADD(T6G, T6H);
			 T6q = VSUB(T5L, T5S);
		    }
		    {
			 Var T69, T6k, T6R, T6U;
			 T69 = VADD(T5T, T68);
			 T6k = VBYI(VADD(T6g, T6j));
			 ST(&(xo[WS(os, 119)]), VSUB(T69, T6k), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 9)]), VADD(T69, T6k), ovs, &(xo[WS(os, 1)]));
			 T6R = VBYI(VSUB(T6P, T6Q));
			 T6U = VSUB(T6S, T6T);
			 ST(&(xo[WS(os, 39)]), VADD(T6R, T6U), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 89)]), VSUB(T6U, T6R), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T6V, T6W, T6l, T6m;
			 T6V = VBYI(VADD(T6Q, T6P));
			 T6W = VADD(T6S, T6T);
			 ST(&(xo[WS(os, 25)]), VADD(T6V, T6W), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 103)]), VSUB(T6W, T6V), ovs, &(xo[WS(os, 1)]));
			 T6l = VSUB(T5T, T68);
			 T6m = VBYI(VSUB(T6j, T6g));
			 ST(&(xo[WS(os, 73)]), VSUB(T6l, T6m), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 55)]), VADD(T6l, T6m), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T6p, T6s, T6F, T6M;
			 T6p = VBYI(VSUB(T6n, T6o));
			 T6s = VSUB(T6q, T6r);
			 ST(&(xo[WS(os, 41)]), VADD(T6p, T6s), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 87)]), VSUB(T6s, T6p), ovs, &(xo[WS(os, 1)]));
			 T6F = VADD(T6x, T6E);
			 T6M = VBYI(VADD(T6I, T6L));
			 ST(&(xo[WS(os, 121)]), VSUB(T6F, T6M), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 7)]), VADD(T6F, T6M), ovs, &(xo[WS(os, 1)]));
		    }
		    {
			 Var T6N, T6O, T6t, T6u;
			 T6N = VSUB(T6x, T6E);
			 T6O = VBYI(VSUB(T6L, T6I));
			 ST(&(xo[WS(os, 71)]), VSUB(T6N, T6O), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 57)]), VADD(T6N, T6O), ovs, &(xo[WS(os, 1)]));
			 T6t = VBYI(VADD(T6o, T6n));
			 T6u = VADD(T6q, T6r);
			 ST(&(xo[WS(os, 23)]), VADD(T6t, T6u), ovs, &(xo[WS(os, 1)]));
			 ST(&(xo[WS(os, 105)]), VSUB(T6u, T6t), ovs, &(xo[WS(os, 1)]));
		    }
	       }
	  }
}
