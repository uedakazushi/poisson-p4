# Actual N tailの同じR0上の8・36・84解像度

2026-11-04。`AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63`の実source/同時R0結合を使い、
前packetで条件付きだったmodule coneをliteralな係数DAGへした。
このconstructorはN bindingそのものの独立再監査ではない。
新しいsource jet、点、GB、prime、source producerの再実行は無い。
D4のnil extensionとT3の非零性は未計算である。

## 1. 係数とgeneratorは全てactual physical functions

同じphysical R0九座標はrows

```
27,108,109,112,113,117,118,119,61.
```

R8の順序はrows108,109,112,113,117,118,119,27である。
既保存の56-unit Hensel系(det 298598400)とquadratic source basisを再利用し、
全150 source functionの正確なp-basis pairとtracefree pairを回復した。
これは新しいjetではなく、既に全次数で指定されたquadratic coefficient入力への
有限の有理演算である。元のR0九pair、selected八odd coefficient nodeが全て文字どおり一致する。

```
tau=p-S/2,  tau^2=nu,
z=(b61,b18,b19,b28,b29,b45,b46,b60),
v=rho61-a61(R8).
```

unit8 determinant -4/9によりzはR8のregular coordinates。
accepted bindingは同じR0で

```
R0=R8[[v]],   A_N=R8+m8*tau,   v acts as z1*tau
```

を与える。nu,a61,z_iは保存56系の具体的なnodeであり、独立な形式係数ではない。
元rho61はR0の九番目の独立変数として保持し、N上でのみv=z1*tauへ送る。

N moduleの八generatorを次の順に固定する。

```
n0=1,
n1,...,n7 = rho_j-a_j(R8), j=18,19,28,29,45,46,60.
```

Nのnormalizationではnj=z_(j+1)*tauである。
全a_jを含む全150 even/odd pairを保存するため、これらのgeneratorはD4側へも同じ
physical target functionとして持ち上げられる。CA/CBのtraceをNのtraceへ同一視しない。

## 2. 全36 relationと84 syzygyのliteral matrix

R8上で1,z1*tau,...,z8*tauは、even free summandとm8のKoszul解像度を持つ。
その自由加群F0のbasisをE0,E1,...,E8とし、v作用を

```
T0(E0)=E1,
T0(Ei)=nu*z1*zi*E0,
Tj=0 (j>=1)
```

とchain liftする。T0は全28 Koszul関係を正確に零へ送る。
R0へ係数拡大した解像度で`v-T`のconeを取るとrank9,37,84を得る。
row E1、shifted E0列のunit -1を消す。
残るeight free generatorsは上のn0,...,n7と一致する。

36列の順序も固定した。

* 最初の28列：z1,...,z8のlexicographic pairに対するKoszul relationでE1=v*n0を代入。
  z1を含む七列は `z1*n_j-z_(j+1)*v*n0`。
  残る21列は `z_i*n_j-z_j*n_i`の対応する添字変更。
* 最後の八列：元のz_i*tauへのv作用。
  最初は `(v^2-nu*z1^2)*n0`。
  残る七列は `v*n_j-nu*z1*z_(j+1)*n0`。

84列のd2は、最初の56列がKoszul3のlexicographic triple、最後の28列がcone側の
shifted Koszul2である。元37行からunit pivotの一行を消したものを保存する。
全d1*d2=0を、actual nodeへ代入する前の完全疎多項式恒等式で照合した。
36関係自体も、`v=z1*tau, tau^2=nu`への代入で全係数零となる。
全matrix entryは、この普遍多項式をactual z,nu,v nodeへ代入したliteral DAGである。

## 3. 正確性と最小性

元のR8解像度はregular sequence z1,...,z8のKoszul解像度から得る。
v作用は位相的nilpotentで、module上でT^2=nu*z1^2。
従ってmonic module divisionにより`v-T`のcokernelはactual A_Nであり、kernelは零。
coneは正確なR0-free resolutionになる。
unitを一つ消去した後の全微分係数はm_R0内にあるので、解像度は最小である。
全Betti ranksは

```
8,36,84,126,126,84,36,9,1.
```

特にactual同じR0上で

```
R0^84 --d2--> R0^36 --d1--> R0^8 --> A_N --> 0
```

は正確である。これは単にd1*d2=0からexactnessを推論したものではなく、
Koszul解像度とmonic coneの証明をactual bindingへ適用したものである。

## 4. D4 extensionへのインターフェース

literal tail packetは同じR0上で

```
0 -> QQ[[T]]*F_old -> D4 -> A_N -> 0
```

を与える。nilcurveはregular8 CIで、physical row117がuniformizer。
上のn0,...,n7は同じtarget functionとしてD4へ持ち上がる。
36個のd1 relationをこのliftへ適用した値は、一意なc_i(T)*F_oldである。
今後nil extractorがこの36個を計算し、curve idealのKoszulでd2との整合を持ち上げれば、
horseshoeでD4のfree ranks9,44,112を得る。

今回c_i(T)、F_oldのH_CA lift、horseshoeのoff-diagonal blockを計算したとは主張しない。
T3の完全presentationへA3 base changeをするのは、このextension/comparisonを完成した後である。
mu(D4)=9という既受理結果だけからextensionの分裂を推論しない。

## 保存

Code `code/ambient_B18_N_literal_module_cone_20261104_v63.py`。
Report/folder `AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63`。
`actual-N-literal-cone-DAG.json`には、元56-unit系、実R0写像、全150 p/trace pairs、
八generatorのphysical/even lift、全36/84 sparse matrix、普遍多項式、順序を保存する。
N binding primaryの独立reviewは別agentのgateであり、このconstructorが代用するものではない。
元B18全Dのdepth、T3の零性、nilradicalの全typeは未決定である。
