# B18四係数加群によるdepth判定：独立構造監査

2026-11-04。**構造的還元を受理する。actual depth19/20は未決。**
対象は同名author packetであり、actual係数の新展開、source/gauge/GB再実行は行わない。
受理済み相対weighted modelとその独立review、共通G19積、core gauge recipeの保存hashを照合する。
小さい有限代数で別のFraction線形代数を使い、係数ideal、w作用のkernel、最上段生成元を検算する。
この算術toyは一般証明やactual係数の代用ではない。

## actual模型への適用

既受理の相対模型は同じphysical19切断をCA/CBとそのnormalizerに適用する。
切断後の単射は共通equivariant群積から保たれ、任意base changeの単射を仮定していない。
共通even baseはraw rows
`27,108,109,112,113,117,118,119,140,61`、正weightは`w=u140`だけ。
その残り九変数の完成正則環をR0とする。これはfold61を除いた別の有限13模型のbaseではない。

CAのactual像は `A=R0[[w]]+(I0+wJ0)R0[[w]]τ`、`τ²=ν∈R0`。
八本のweight0 odd、二本のweight1 odd、τのweight0は全次数のequivariant theoremを入力とする。
したがって `C0=R0+I0τ`、`H=R0+(I0+J0)τ` とすると、
`A0=C0`、`An=w^nH (n≥1)` が正確に成立する。
ここでC0,HはR0²内の有限部分加群であり、freeである必要はない。

CA kernelで先に割り、CB kernel高々99本と**全87 graph差**をAへ送った閉idealをQとする。
これでactual `D=A/Q` を保つ。graph重みは4/53/24/6、Q生成元上界は186、
degree3はraw rows `9,39,54,84,99,144` の六本である。
旧abstract99模型同士の相対alignmentや低jetからの全次数推測は使わない。

## 完成、安定tail、torsion

degree d生成元の係数のC0-spanをLdとする。有限生成homogeneous idealのdegree n成分は
`Q[n]=Σ_(d<n) H Ld+Ln`（n≥1）、`Q[0]=L0`。
従って `J=Σ_(d=0..3)H Ld` とおくと、n≥4で `Q[n]=J`。

完成に関する抜けはない。AはNoetherian complete環R0[[w]]上有限で、Qは有限生成ゆえ閉である。
固定weight成分は有限和で決まり、逆に各degreeで選んだ有限個生成元の係数をw-adicに足せる。
R0加群も有限で閉じており、(m0,w)-adic topologyで同じ完成を与える。
従ってdegree≥4の無限tailに、成分検査から見えないw-torsionは残らない。

`T=H0_(w)(D)` は正確に

```
T0=(C0∩J)/L0,  Tn=J/Q[n] (n=1,2,3),  Tn=0 (n≥4)
```

となり、四有限R0加群の直和、`w^4T=0`である。degree0だけ係数のambientがC0である点も正しい。
有限商Dは同じ点でlocalなので `sqrt((m0,w)D)=mD`。
w⁴がTを消すことから
`H0_mD(D)=H0_m0(T)=⊕_(n=0..3)H0_m0(Tn)`。
非零な有限長部分加群は非零socleを持つため、四つのsocleのいずれか非零 iff depthD=0。
四つとも零ならdepthD≥1、受理済み上限1によりdepthD=1。
共通G19を戻すとcone depth19/20、scalarを除くとprojective depth18/19となる。
これは実際のどちらの値かを決めたものではない。

R0²内で各分母加群をm0でcolonし、分子と交叉して分母で割るauthorの四式は、
上のsocleそのものである。m0-colonは九つ全てのparameter条件であり、一つのparameterや
有限jetのvanishingに置き換えない。Tはnilradicalではなく、被約枝上でもwが消える場合を含む。

## 最上段12生成元と小例

H/C0はj1τ,j2τの高々二元でR0生成される。
C0のodd部分を掛けた積はeven部分R0へ入りC0に戻るため、この生成上界はC0作用を加えても増えない。
従って六つのcubic係数cαと掛けた12元で
`T3=(Vlow+H L3)/(Vlow+L3)` をR0生成できる。
`cα=aα+bατ` に対するpairは `(ν j_e bα,j_e aα)` である。
Vlowを捨てない。T3≠0だけではsocle非零を意味せず、T3=0だけでも四検査を完了しない。

別実装の有限環 `k[x,τ]/(x^N,τ²-x)` を用い、四つの異なるI0,J0とmixed係数idealについて、
degree0〜8のactual idealの積を直接生成し、w⁴kernelとtail、上の最上段R0生成元をFractionで照合する。
鋭さの例 `I0=(x),J0=(1),Q=(w³τ)` では同じclass `[τ²]` がdegree0〜3を通って残り、
`w³[τ²]≠0,w⁴[τ²]=0`。したがってdegree3を残す仮定だけでw³へ改善できない。

次に必要なのは同じactual CA trace基底でのν,I0,J0と全Q係数pairの全次数identityである。
core gaugeは23本のunit implicit systemとmonic quadraticによるalgebraic power-series recipeを与える。
係数は今回未評価で、ungauged rank7をrank8とみなす逆写像や有限jetのmembershipは受理しない。

再現コマンド:

```
python code/ambient_B18_four_module_depth_criterion_review_20261104_v63.py \
  --output /tmp/B18-four-module-independent-review
```
