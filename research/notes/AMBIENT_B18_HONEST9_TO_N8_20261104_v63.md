# Honest9から同じR0上のN module8へのactual有限map

2026-11-04。nil generatorをCA二行basisへ戻すgateを使わず、honest9 generatorの
被約N上の像を実際の有限unit DAGで構成した。
nil補正行、9×9全transition、extension coefficients、T3の零性は今回未計算。

## 1. sourceとtargetのtraceを混ぜない96単位base change

honest9のgeneratorは、単位Schur後に残る

```
1, row18-a_CA18(R0), row19-a_CA19(R0), row28-a_CA28(R0),
row29-a_CA29(R0), row67-a_CA67(R0), row72-a_CA72(R0),
row46-a_CA46(R0), row95-a_CA95(R0).
```

ここでa_CAはCA自身のtrace-even関数であり、純粋なCA48単位系で定義されるR0 scalar関数。
これをNへ移すと、Nのquadratic basisではodd成分を持つ場合もあるので、
単にCA evenとN evenを等置しない。

Nの同R0結合は、actual R0九関数のN basis(1,p_N)での全pairを与える。
これをCA48式のR0変数へ代入し、CA48個の未知scalar関数それぞれを
Ai+Bi*p_Nと置いた。Nのp_N²=S_N*p_N+T_Nを使って偶奇係数を全て保持すると、
96未知数への96式となる。全分母は原点でunit。
Jacobian determinantは

```
7213895789838336/390625 = (J_CA48)^2 != 0.
```

この96系の全式とJacobianを保存した。unit Henselの一意性により、
これは元CA scalar関数をN normalizationへbase changeした実際の関数である。
CAのtauをNのtauと同一視していない。

各honest generatorから得たN pairをtracefree basis tau_N=p_N-S_N/2へ直し、
`h_j=a_j(R8)+b_j(R8)*tau_N`とした。全b_jの原点は0。
R8とtau_Nは既受理のactual N同R0結合のものである。

## 2. 実際のunit8座標でoddを分解する

z1,...,z8はphysical odd rows

```
61,18,19,28,29,45,46,60
```

のN odd係数であり、R8のunit座標である。
N56系、z8座標式、CA96系を合わせた160系の変数順を
N56未知数・元R8八変数・CA96未知数とした。
Jacobian determinantを実際に計算し、block積とも一致することを確認した。

```
det J160 = -38294359833110460235776/15625 != 0.
```

七つのprefix160系をliteralに構成した。
八つ全てのzを残すinverseは元の関数とunit一意性で一致するので、新しい八番目の系は不要。
各b_jについて

```
c_ij = [b_j(z1,...,zi,0,...,0)-b_j(z1,...,z_(i-1),0,...,0)]/zi
```

を64個のregular coordinate quotientとして指定した。
分子はzi=0で厳密に零となるため、全次数で正則かつalgebraicである。
原点値も160系のimplicit ADから求めて保存した。
telescopingにより、全次数の恒等式`b_j=sum zi*c_ij`が成立する。
有限jetからmembershipを推測する操作ではない。
旧full186のr/g商は一切composeしておらず、全関数はpure unit系で定義される。

## 3. B_Nの全8×9係数

同R0上のN moduleの順序は既存literal coneに合わせる。

```
n0=1,
n1,...,n7 = rows[18,19,28,29,45,46,60]-それぞれのN even関数。
```

v=rho61-a_N61(R8)=z1*tau_NはR0 scalarである。
従ってhonest generator h_jは正確に

```
h_j=(a_j+c_1j*v)*n0 + c_2j*n1 + ... + c_8j*n7
```

へ送られる。定数generator1のcolumnは(1,0,...,0)。
この全8×9 matrix B_Nをliteral DAGとして保存した。
pair恒等式は、N trace arithmeticと上のtelescopingのみによって全次数で証明される。
実際のphysical generatorから構成したmapなので、D3/D4から被約Nへの商mapに一致する。

原点matrixはrank8である。constant kernelは、honest generator順で

```
(0,0,0,0,0,0,0,1,1)
```

すなわち最後のphysical rows46+95方向。
これはN moduleへのmapのkernelの一次情報であり、T3の非零性ではない。
nil補正を加えた9×9 transitionの単位性も、このmatrixだけから計算したことにはしない。

## 4. 残るnil補正の扱い

各honest generatorから上のN式をactual D4で引いた差は既知のnil idealに入る。
その係数を決めるには別のnil derivation extractorが必要である。
新しいcoordinate quotientをnilcurveへcomposeする際、ziの像が恒等零なら
分子・分母を0/0として再除算してはならない。元の正則関数のpullback、または
そのpure unit prefix系を微分した有限unit入力を使う必要がある。

ただしtransitionの中心だけならHadamard係数の微分は不要である。
n0=1、他のn_iは中心で零、vも中心で零なので、中心のnil補正には
`delta(a_j)+c_1j(0)*delta(v)`と各generatorのdeltaだけが寄与する。
このscopeの違いを後続extractorへ引き継いだ。

## 保存・実測

Source `code/ambient_B18_honest9_to_N8_20261104_v63.py`。
`reports/AMBIENT_B18_HONEST9_TO_N8_20261104_v63/honest9-to-N8-DAG.json`に
CA96、全160 Jacobian、七prefix系、64商、全B_Nと原点matrixを保存。
全675735 scalar node、証明書13377046bytes。
一回2.551398秒、最大RSS366556KiB。60秒/2GB/32MB内。
jet、GB、新点、p逆元、old nonunit商のcompose、opaque branch pullbackは無い。
元B18原点のfull D depth/nil、T3 presentationや零性は主張しない。
