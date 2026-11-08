# h5 residual P*：pure-even MPS源と四hyperplane-log quadratic jet の whole 排除

2026-11-08、PRIMARY E_multiple_hyperplane。純全源監査。
新しい点、operator、rank/minor、jet延長、GB又はCASを実行しない。
既印刷全sourceと、既強制点の唯一order二zeroに関する独立受理済み固定式だけを用いる。

## 1. 入力となるP*と既global marked-H条件

`UNIMODULAR_POLE_GLOBAL_MARKED_JET_AND_T36_T38_EXCLUSION.md` を入力にする。
P* は全 h5源の既未帰属一軌道で、D=h5、projective tangent24／orbit22。
唯一のPoisson hyperplaneはh=0、そこでrank[dQ|L]=5。
同noteのproper marked-h theoremにより、全sourceがあるPoisson Hのaffine H=1で
C+Q（L=0）となる既知componentは全閉包を含めP*を含まない。
この lemma は rank T_Q が落ちる補集合も弱い5角minorで保持し、
一般点でoriginを動かせることだけから特殊点を有界frameへ外挿していない。

## 2. Whole pure-even source の直接table

raw又はtrue completed sourceでHを下表の座標に選ぶ。
全H行はHで割れ、そのquotient X_i=P_iH/H はHを含まないspatial linear field。
P_H=Praw+E wedge Xで全H行を消す。
残る全entryはH次数0又は2だけであり、H=1でC+Q、L=0。
この判定は全period/edgeをsymbolicのまま全monomialに対して行う。

| source | H | Hを含む非toric coefficient（他はspatial二次） |
|---|---|---|
| T02 | x4 | t0H² |
| T03 | x2 | 無し |
| T04 | x3 | 無し |
| T06 | x3 | t1H² |
| T07 | x3 | 無し |
| T08 | x4 | t0H² |
| T09 | x3 | 無し |
| T13 | x3 | 無し |
| T15 | x3 | 無し |
| T16 | x4 | t0H²,t1H² |
| T17 | x4 | t0H²,t0t1H²/p1 |
| T20 | x4 | 無し |
| T21 | x3 | 無し |
| T22 | x4 | t0H²,t1H²,t2H² |
| T23 | x4 | t0H²,t1H²,-3t1t2H²/p0 |
| T25 | x4 | t0H²,-3t0t1t2H²/p0² |
| T28 | x4 | 9t1t2H²/p0,t2H² |
| T29 | x4 | t1H²,t2H² |
| T30 | x4 | t1H²,-t1t2H²/p0 |

T17/T20/T22以降は**true completion**を含め、denominatorを勝手に零へ代入しない。
必要 source denominator 非零の全generic sourceでこの恒等式が成立し、
proper marked-H theorem の閉性により同じfull componentの全境界へ延びる。
従って表の十九componentはwholeにP*の親でない。
sourceが一つのh5点を持つというだけの主張ではない。

## 3. 四linear log one-form閉包における triple-pole square lemma

四次積I=ell0*ell1*ell2*ell3をintegrating factorとし、

```
alpha/I=sum rho_i*dlog(ell_i), sum rho_i=0
```

となるcubic horizontal one-form全源のprojective閉包を考える。
residueは任意であり、一点の整数ratioを全sourceと呼ばない。

P*の唯一order二zero pでのcubic alpha*は既入力により非零で非定数rational first integralを持たず、
quartic integrating factor I*=x*h³ を持つ。
二つのintegrating factorの比はrational first integralとなるので、I*はscalarを除き唯一。
properな四hyperplane flagsをDVRへ延長すると、targetでは一つがx、残る三つがhへ行く。
同一のregular PGL変更で単純flagをx、一つのcluster flagをhへ正則に固定できる。
他の二つをh+a,h+bと書く。a,bはspatial linear formsでDVRのmaximal idealに属する。
原点h=1のformal spatial expansionでは、xのlog residueはtargetで有限なので除け、
clusterの一次・二次moment

```
M1=lambda*a+mu*b,
M2=lambda*a tensor a+mu*b tensor b
```

は共にcoefficientwise bounded。
これはclosed one-form familyが単元denominatorで正則specializeすることから従い、
residueが有界であると仮定していない。

任意の二spatial indexについて

```
a wedge M1=mu*(a wedge b),
b wedge M1=-lambda*(a wedge b)
```

はa,b→0とM1 boundedから零へ行く。
Cauchy–BinetによりM2の各2角minorは
`lambda*mu*(a wedge b)_ij*(a wedge b)_kl`、すなわち上の二つの零へ行く量の積。
よってtargetのM2はrank≤1。
clusterのpole order三が非零なら、そのhighest quadratic exact numeratorもrank≤1となる。
これには全不等valuation、residue発散、階層的cluster、単元frame変更を含む。

ところがP*ではh=1 chartで

```
alpha*/(x*h³)=dlog(x/h)-d(xy/h²+z/h),
```

quadratic exact numeratorはxyでrank二。
よってalpha*は四linear-log全source閉包の外である。
zero/degenerate alpha targetは使わず、unique integrating factorと非零rank二を保持する。

## 4. Whole T04–T15への適用

T04,…,T15の全true sourceはx0に依存しない。
q=[1:0:0:0:0]でprojective0/1jetは零、quadratic jetはW=(x1,x2,x3,x4)上の
純logcanonical Qij=cij*x_i*x_j。
全smoothing/correctionはvertex行P01/P02だけにあり、affine vertex chartではcubic以上を変えるだけ。
従ってquadratic jetのcubic1formは

```
alpha_Q/(x1*x2*x3*x4)=sum c_i*dlog(x_i)
```

であり、§3の全四linear-log classに入る。
独立四hyperplane上の源からproper flag/closureを取った全closed coneを使い、零alphaも含める。
order二zero qとこのcone条件のuniversal incidenceは閉、qはP4でproper。
よってT04–T15の全component closureに「あるqでそのjet条件」が成立する。

P*はorder二zeroがpだけで、alpha*は§3によりそのconeの外。
従って**全T04–T15はP*を含まない**。
§2と重複する七IDを二重に数えない。

## 5. 残す親と射程

§§2–4のunionは T02–T17（T18/T19を含まない）、T20/T21、
T22/T23/T25/T28/T29/T30 である。
T18/T19、T24/T26/T27/T31/T32/T33/T34/T35、T37はこのnoteでは未排除。
T01/T39/T40などdimension25以上は既P*接上界24の別定理、
T36/T38は既入力のwhole jet/zero-surface定理に依存する。
全MPS排除、P*成分性、actual parent帰属又は全E完了をまだ主張しない。
