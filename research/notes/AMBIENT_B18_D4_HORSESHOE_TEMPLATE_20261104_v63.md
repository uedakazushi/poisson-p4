# D4 horseshoeの有限templateと、relationを持ち上げるcone homotopy

2026-11-04。actual N cone8/36/84を用いた有限演算のtemplate。
このtemplate単独では36 residue係数をまだ代入していない。
T3のstructural零性証明とは独立で、今後の低次数module比較に使える。

## 1. 36 residueが与えられたときのliteral44/112

同じR0座標を(z1,...,z8,v)、nilcurve parameterをZ=z1とする。
actual N relation alphaのD4 liftへの値を正符号でc_alpha(Z)*Fとする。
curveのgraph idealは

```
P=(z2-psi2(Z),...,z8-psi8(Z),v-psiv(Z)).
```

D4 generatorの順をF,n0,...,n7とする。d1の最初の八列はP_i*F、次の36列は
`Nrelation_alpha-c_alpha(Z)*F`。

Nのd2全entryはz1,...,z8,vの一次式だけで、nuは現れない。
各N syzygy betaについて

```
h_i,beta=sum_alpha c_alpha(Z)*[coefficient of x_i in d2_N(alpha,beta)],
x=(z2,...,z8,v).
```

と定義する。Nrelation間のsyzygyをD4へ持ち上げれば
`(sum_alpha c_alpha*d2_N(alpha,beta))*F=0`。
FはR0/P上free basisなので、このsumのcurve restrictionは零。
d2_Nが一次式だから正確に

```
sum_alpha c_alpha*d2_N(alpha,beta)=sum_i P_i*h_i,beta.
```

従ってd2_Dは、nilcurveのKoszul28列と、上段h・下段d2_Nの84列、計112列である。
新しいHadamard係数もsource jetも要らない。
templateは全112個のd1*d2を整数多項式で照合し、84個のcurve cocycle以外の残差が無いことを示す。
実際のc,psiを単位DAGで与えれば、`build_horseshoe`がliteral matrixを返す。

正確性は0→L→D4→N→0のhorseshoeで従う。
下段がNの全関係を尽くし、残る上段はPの全関係であり、Pのfirst syzygyは全28 Koszul。
従ってここではsyzygyを数値rankから推測していない。

## 2. N coneの有限lifting homotopy

旧original free9のbasisをE0,...,E8とする。
v作用のchain lift T0は、厳密に

```
T0^3=q*T0, q=nu*z1^2,
(v-T0)*(v^2+v*T0+T0^2-q*I)=(v^3-q*v)*I
```

を満たす。全matrix係数を疎整数多項式で確認した。
Nへ零に写る任意のfree vector h(v)について、monic cubicでWeierstrass divisionし

```
h=(v^3-q*v)*k+a0+v*a1+v^2*a2,
e=a0+T0*a1+T0^2*a2,
b=(v^2+v*T0+T0^2-q*I)*k+a1+(v*I+T0)*a2
```

とすれば、h=e+(v-T0)*b。
eはR8係数のvectorで、even成分0、odd成分はsum zi*e_i=0を満たす。
regular z8のKoszul contractionでe=d_K*aと書けば、(a,b)がoriginal coneへのliftとなる。
unit pivot列を除けば、literal36-column N presentationへのliftを得る。

Koszul contractionも有限の規則で指定できる。
z1=0でtail8-1成分を再帰的に処理し、引いたtail成分をz1で割る。
分子のrestriction-zeroが証明されているので、その商は正則である。
これを八変数まで繰り返す。

## 3. Old2069を持ち上げる際の実装上の残り

B_fullでold honest9 relationを(F,N8)へ移す。
N成分を上のcone homotopyで36列へ持ち上げ、それを引いたscalar F成分を
nilcurveのregular8 graph ideal PでHadamard分解する。
これで44列へのliftが得られ、完全な112 syzygyと合わせればT3の有限presentationになる。

この手順で旧r/g coefficientは**元のR0の正則functionのまま**扱う。
monic divisionやregular-coordinate pullbackを、r/gの分子・分母を個別代入して
再除算する方法で実装してはいけない。零/零になることがある。
正則algebraic power seriesとしてのpointed母系、またはそのexact division operatorを保持する。
今回このWeierstrass/Koszul operatorを全2069列へ実行していない。
従って計算済みのrelation liftや追加の有限jet情報として主張しない。
monic cubic divisionの有限operator formulaは、未知SyzJinfを単に要求し直すcriterionとは異なるが、
そのcoefficient-DAG評価を完了したという意味でもない。

Code `code/ambient_B18_D4_horseshoe_template_20261104_v63.py`。
Report/folder `AMBIENT_B18_D4_HORSESHOE_TEMPLATE_20261104_v63`。
新しい点、source solve、jet、GB、primeは無い。
