# 非被約R・m三squarefree残二Ore源：全支持の別作者確認

2026-11-07、`/root/quartic_review`。**全source幾何と独立Fraction算術はPASS。**
対象はrank二Poisson plane、R=2L+M、normal位数三、generic normal leadingがsquarefree、
actual Dが被約既約quinticである枝。
既十九源reviewで残ったDIAGONAL_-2_-2及びNIL_SPLITを閉じる。
全m三squarefree非被約RのPGL像は次元二十四以下となるが、
c-dependent/g二の全点の既成分帰属、normal repeated-root、m二及び他Sは結論しない。

## 1. 不変の入力と今回の独立算術

producer public reportは `CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_20261107.json`
（c3ea0064…）。そこに埋め込まれた実行source60491daf…をhash確認し、
mutable current runner又はcurrent producer jobsを入力同一性に用いていない。
旧十九源はNEW_REPLAY bf8be1a3…と別作者receipt d0c2ff07…を利用し、再実行・改変しない。

新checker `classification_plane_double_line_m3_ore_review.py` はproducer/CASをimportせず、
既独立のFraction/Laurent sparse演算だけから式を再構成した。
有限job `jobs/plane_double_line_m3_ore_separate_review.json` のgate PASS後、
**0.779秒、44,208KiB**で次の全検査に成功した。

* Dの全P1二十＋全P2九、Nの十五＋九を、全十一又は六Q係数＋全W十八へ移すliteral Euler gauge。
  両完全係数行列の行列式は−1。旧完全K1との係数一致、全六P2 magnetic項も含め全九P2を保持。
* 任意全Q係数のstrict Jacobi、NのPfとsquarefree判別式、実normal shear、
  diagonal/scalar変換及び全central coboundaryの恒等式。
* 四canonical Qの全350×18 affine operatorと定数、全140×18 Ore moduleと定数。
  全行の符号付き対応を検査し、numerical rank比較でなく全affine系の同値を証明。
* 全rank十四minor、全十八行の完全kernel四方向とその単位free block、
  augmented minor及び完全kernel、全行に対するparticular solution。
* actual Euler-signed Pfをparticularと全kernel方向で再構成し、全affine familyのproper共通因子を確認。
* whole complex tau line。generic matrix/augmented minorsの積は4294967296 tau²。
  全factor-specialization分母は1又はtauであり、唯一wall tau零を元の全350/140から別確認。

新report `CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_SEPARATE_REVIEW_20261107.json`
にchecker全文、全input hashと各検査結果を保存した。
成功report及び旧reportのbytesは以後不変とする。

## 2. 全Ore表示とA零補集合

座標は(a,b,c,z,w)、y=(a,b,z,w)。旧完全sourceからcを含む項はac/bcの
cA(z,w)、cB(z,w)だけである。L=−A∂a−B∂bとしEuler_total∧Lを引くと

`Π=Q(y)+∂c∧W(y)`

を得る。元の全係数との行列式−1によりchosen sliceへの制限ではない。
全projective Jacobiをc次数で比較すると、c係数はstrict[Q,Q]であり、
これが零の後の残条件は正確に `E_y∧[Q,W]=0`。
projective Schoutenとstrict Jacobiを混同していない。

DではAを全binary linear、Xを全binary quadratic vectorとし
Qab=−bA+f2、Qan=−A E_n、Qbn=X、Qzw零。
strict Jacobiはz X(A)、w X(A)、零、零で、全f2は残る。
**A零をnormalizationで捨てない。** この場合Q=∂b∧(X−f2∂a)、L_∂bQ零。
一般補題として、Q=v∧Z、L_vQ零ならPfQ零なのでactual Dはc非依存。
module条件をvとwedgeすると、v∧[W,v]∧Zの項だけが残り、
`v(D)=vol(E_y∧Q∧∂vW)=0`を与える。
従ってここではDはbにも依存せずessential高々三。
被約既約quinticはbinaryにはなれないのでessential三であり、
既受理の八free-quintic全sourceに入り、全像≤23。
これは大きいrank二subfamilyを無理に次元計算する議論ではない。

## 3. 全actual normalizerとcentral classの補集合

DのA非零ではnormal GL2でA=z、X(A)零からX=R2∂w。
leading squarefree条件はR2 squarefreeかつz∤R2と正確に同値。
zを保つnormal Borelは残二根をw零、w=zへ送り、b-scaleが振幅を一に戻す。
全inverseはこのsquarefree openでactual linear変換である。

a,bへの全linear normal shearのcentral変化は−X(u)。そのimageはR2の一方向だけなので、
残class αz²+βzwを全て保持する。
producer記載のactual diagonal/scalarはcanonical anchorsとfirst symbolを保ち、
f2全体をr/s倍することを独立算術でも確認した。
従って全classは零、α零/β非零、α非零の三strataであり、
f2零、zw、z²+tau zwの全complex lineが必要十分。
tauをgeneric値へ置き換えたりfinite root permutationsで補集合を捨てない。

Nでは全strict Q Jacobiが恒等零、PfQ=−(v²+12ur)z⁴/6。
full leadingの判別式は正確に729 a⁸u²(v²+12ur)。
従ってu(v²+12ur)非零が対象SF枝の全openであり、零補集合は対象leadingを持たない。
actual w→w+hzでv→v+6hu、r→r−hv−3h²u、Δ=v²+12ur不変。
h=−v/(6u)でv零、r=Δ/(12u)非零とし、
記載diagonal/scalarの全式によりnonzero scalesと有限square-root choiceでu=r=1へ戻る。
これは単にconnectionを同値化するのでなく、P0も保つactual PGL/tensor操作である。

canonical Nで全a,b shearのcentral像は
`(gamma−beta)z²+(2delta+6alpha)zw+9beta w²`。
独立checkerは係数subminor54を確認したので全f2を除去できる。
Nに余分なcentral modulus又は分母零枝は残らない。

以上のfiltered linear操作はfull higher sourceにも同時に作用する。
変換後のWは全positive weight二十係数にあり、
正当なEuler gauge `W→W+ell(z,w)E_y` の二係数で
Wzのz²とWwのw²を除く。残りは使用した全W十八そのもの。
従ってunknown higher sourceをnormalizer途中で落としていない。

## 4. 四全moduleとwhole tauの最弱結論

| canonical Q | matrix/augmented rank | 全fiber | 全actual Dのproper因子 |
|---|---:|---:|---|
| D f2零 |14/14|4|z w(w−z)|
| D f2=zw |14/14|4|z w(w−z)|
| D f2=z²+tau zw、tau非零全体 |14/14|4|z w(w−z)|
| D tau零wall |14/14|4|z w(w−z)|
| N fixed |14/14|4|z|

actual DはWにaffine linear。
particularと全四kernel方向のDを割るproper factorは全solutionのDを割る。
未知magnetic係数の一部を零にした点だけのfactorizationではない。
正確なwhole tau!=0 identityと別のtau零certificateにより全complex tauを覆う。
全四normalform源はactual被約既約quinticを持たず、inverse actual変換後も同じ空性である。

A零のessential三枝だけを上記八23成分へ送り、残四源は空とするこの結論が、
全ambient completed ring又は追加jetより弱く十分である。

## 5. aggregateの適用範囲

既十九源別作者receiptは十六leading-empty、c-dependent/a²/g二の全像≤24、
及び本稿で閉じた二Ore源を尽くす。R=3Lの既受理全二モデルleading排除も合わせれば、
**非被約R・m三squarefree全枝のPGL像は二十四以下**。
有限個のconstructible像の閉包も同じ上界なので、未知dimension25–27の優越成分の一般元はこの枝にない。

c-dependent/a²/g二について全pointwise E帰属は依然未証明である。
本稿は全G、m三repeated leading、m二、reduced R、他のrank二surface又は新成分数を主張しない。
