# T0のR-primary部分を分ける：canonicalな二重DVR quotientと残るsocle問題

2026-11-04。これは有限source計算の追加ではなく、既存のactual mapを使うmodule論である。
R0は現在の共通九変数環、D0はactual零次数algebra、T0=ker(D0→D4)とする。

## 1. Saturationに使う二つの実際の函数

同R0上のN moduleから

```
F=v^2-nu*z1^2,   v=rho61-a(R8),   z1=b61(R8)
```

を取る。N projectionのsupportはV(F)。nilcurveもこのsupportに含まれるため、
Supp_R0(D4)=V(F)である。Rのactual common-cut image P_Rは八つのregular graph式で指定され、
R0/P_R=QQ[[U]]、U=rho112である。先の19+56 scalar-unit計算により

```
F|R=-81600 T^4+O(T^5),   U|R=-2T+O(T^2).
```

従ってFはR上の零函数ではなく、P_R+(F)はm_R0-primaryである。
全positive pieceがD4に等しいため、R generic pointではD0=T0が全Dそのものとなり、長さは正確に2。

## 2. D0の被約support

元Dの被約supportはNとRで尽くされる。非負torus gradingのdegree0抽出はexactであり、
Nil(D0)=Nil(D) intersect D0である。したがってD0_red=(D_red)_0。
Nのdegree0 subalgebraのR0 kernelはFであり、Rのdegree0部分は上の固定DVR curveである。
よって有限R0-module D0のsupportはV(F)とV(P_R)の合併である。
この議論はT0がcurve-supportedであることを仮定しない。

## 3. R-primary quotientは自動的にCM

Noetherian性により、次の増大列は有限段階で停止する。

```
S_F=0:_{D0}F^infinity=union_n(0:_{D0}F^n),
B_R=D0/S_F.
```

B_RにはFが非零因子として作用する。saturationのassociated-prime公式により、
Ass(B_R)はAss(D0)からFを含む素idealを除いたものとなる。
上のsupportのうちV(F)側は全て除かれる。V(P_R)上では、R0/P_RがDVRであり、
P_R以外の唯一のprimeは閉点で、そこはFを含む。よって

```
Ass_R0(B_R)={P_R}.
```

B_Rは非零な一次元CM moduleである。Uはそのparameterかつ非零因子であり、
B_R/U B_Rは有限長なのでB_RはQQ[[U]]上有限freeとなる。
generic R長2とR0/P_R=QQ[[U]]からそのrankは**2**である。

さらにS_FはD0のidealなのでB_Rはalgebraで、被約化はQQ[[U]]。
実際のscalar Uがこの被約写像のsectionを与える。nilradical L_RはQQ[[U]]上rank1のfree module。
L_R^2はgeneric ordinary-double ringで零となるが、B_RはU-torsionfreeなので全体で零。
従って抽象的なcomplete algebraとして

```
B_R isomorphic QQ[[U,epsilon]]/(epsilon^2).
```

epsilonの元physical座標での表示、R0作用の全embedding、停止するsaturation指数は未計算。
原点で元D0そのものがこの二重DVRだと主張したものではない。

## 4. T0とsocleへの正確な帰結

D4はFの外で零なので、T0_F=D0_Fである。よってT0のB_R内の像I_Rはfinite-colength idealであり、

```
T0/(0:_{T0}F^infinity)=I_R subset B_R.
```

このI_RもQQ[[U]]上free rank2。ただしその二生成式・colengthやT0全体のflatnessは未決である。

B_RとI_Rには有限socleが無いので、二つの完全列から

```
H0_m(D0)=H0_m(S_F),
H0_m(T0)=H0_m(T0 intersect S_F).
```

D0/T0はdepth1のD4のsubmoduleであるためH0_m(D0)=H0_m(T0)も既受理。
したがって全て同じ有限socleを測る。R-primary CM部分はcanonicalに除けたが、
残るS_FやT0 intersect S_Fが有限長・curve-supported・socle-freeとはまだ分からない。

## 5. 次に必要な有限入力

socleを閉じるには、actual D0の完全9/640表示において、F-torsion部分を抽出し、
そのm-primary torsionを判定する完全なmodule計算が必要である。
一次matrixのrank38やvisible gA*tau classがsocleでないことだけでは、この判定を代用できない。
少なくとも有限saturation指数、またはS_Fに対する十分なactual annihilator/presentationを証明しなければならない。
既存612 Q0 cofactorを保ったmapping coneは一つの正当な経路であるが、ここでは未実行。
新しいjet、point、prime、GBや大きいcomparison jobは行わない。

Code `code/ambient_B18_T0_R_primary_separation_20261104_v63.py`。
Report `reports/AMBIENT_B18_T0_R_PRIMARY_SEPARATION_20261104_v63/report.json`。
