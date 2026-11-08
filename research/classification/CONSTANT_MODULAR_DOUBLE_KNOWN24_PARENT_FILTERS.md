# constant modular DOUBLE canonical：既知24親の純構造filter

2026-11-08、PRIMARY E。新算術・点・operator・minor・jet無し。
DOUBLE canonicalのexact算術は全D101/M48、従って射影接24、actual orbit23。
そのproducer証明書を固定入力とする。ここで独立算術を実行したとは称さない。
新しく手で確認するのはMPS全sourceの因子／actual torusとproper closureへの矢印である。

## 1. 大きい親と小さい親の正確な扱い

接次元24は、この点を含む全成分の局所dimensionを24以下に抑える。
したがってCA/CB28及び全25成分は親にならない。
三つのlog/period表の次元上界や一点のgeneric表記からではなく、full原Dの接上界を使う。
一方orbit23なのでdimension22以下の既知成分には入らない。
dimension23の親があれば、そのcomponentはこの23軌道の閉包でなければならず、
その一般元のD／canonical cofactorと直接比較する必要がある。
以下は**既知24候補だけ**のfilterで、全prior23比較の代用ではない。

## 2. MPS T02/T03：全Dは3H+Qのproper閉条件内

正しい原reportで全projective Jacobi TRUEのperiod/smoothing源を用いる。
どちらも唯一の0行項はP01=t0*x4²、他は四変数のlog項Pij=kij*xi*xj。
元Euler-signed Pfを手展開すると

```
D=x2*x3*x4*[κ*x0*x1+t0*(k34-k24+k23)*x4²],
κ=k12*k34-k13*k24+k14*k23.
```

である。従って3H+Qのproper multiplication像に入る。
このPGL不変閉条件は、MPS局所全deformation定理が与えるSNC近傍のdense源から
各全component閉包へ延びる。chosen sliceの像だけで全componentを置き換えない。
targetのD=h²*C3はC3 primeでhを割らず、linear primeの総重複度は二だけなので、この閉条件の外。

## 3. T16/T17：全DはH+2Qのproper閉条件内

T16の同じ全Jacobi TRUE源でu=x0*x1、v=x2*x3、s=x4²と置くと

```
D=x4*(κ*u*v+a*u*s+b*v*s+c*s²),
κ=(p2-p1)*(p1+p2+2p0)/4,
a=(p2-p1)*t1/2,
b=t0*(p1+p2+2p0)/2, c=t0*t1.
```

κc=abなのでκ≠0開で二つのquadricへ因数分解する。
この開はperiodの非退化源に稠密であり、proper factor incidenceによって全境界へ延びる。
T17はnaiveでなく**true completed** reportを用いる。
その全項t0*t1*x4²/p1を保持すれば厳密に

```
D=x4*(p1*u+t0*s)*[(p0+p2)*v+t1*u+(t0*t1/p1)*s].
```

p1≠0はperiodの元unitであり、零境界も同じproper closed factor条件で保持する。
従って両全componentはH+Q1+Q2に入り、target h²*C3にはprime cubicがあるため入らない。

## 4. T21：二つのactual torusが全closureのorbitを22以下へ抑える

元全sourceは、0/1/2 blockの全三smoothingと全二periodを持ち、
3/4方向とのcrossはxi*x3又はxi*x4、P34=p0*x3*x4だけである。
したがってx3∂3、x4∂4は全period/smoothing係数でPoisson場となる。
両者はprojective linear場として独立でEulerの倍にはならない。
実Cstar² subgroupが全source tensorを固定するので、その全PGL orbitは24−2=22以下。
actual orbit tangent matrixのrank上界は閉条件なので、SNC近傍dense源から全component閉包にも延びる。
targetは保存されたM48 nonzero minorでorbit≥23であり、T21全closureの外である。
単なるPfの2H+cubic型一致からT21へ所属を推論しない。

## 5. 他の既知24と残す四候補

CPLUS_EXTは全D=H*Q1*Q2のproper閉条件、DOUBLE_RANK3_TWZEROはH*Q²、
DEPENDENT_3H_Qは三linear factor+quadricのproper閉条件である。
どれもtarget h²*C3（C3 prime degree三）を持たない。
前者のgeneric orbit22も整合するが、因子条件だけで十分である。

固定した現24 projectionの十二件から、以上八件を除くと
**T36、T37、T38、GA9_CUSPIDAL_PLANE_QUARTIC_24**だけが残る。
これはそれらへの所属ではない。別の新24親又はtarget自身の23成分性も未証明。
Ga9のa=1 witnessがD100であることからwhole component D≤100とは推論せず、Ga9候補を保持する。
同様にT37 genericにdoublezeroが無いというopen条件から、その全boundaryを除かない。

元SNC-local全sourceのdensity、MPS成分同定は既受理のpublished-local-deformation入力を用い、
ここで未完成naive T36/T37 sourceや全βcoverageを代用しない。
独立pure review／元exactD/Mの別作者再構成は統合側で別途記録する。
