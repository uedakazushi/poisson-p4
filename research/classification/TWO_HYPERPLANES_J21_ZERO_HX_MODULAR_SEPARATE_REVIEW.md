# Jordan21・(0,h+x)のn6非零全六源：modular実変形の独立確認

2026-11-08、PRIMARY E、`/root/log_manuscript`。**当該全六源のMPS40被覆PASS**。
producer `TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_CONFLUENCE.md` の全source、実曲線、
全半単純supportへの適用を別作者として確認した。新成分は零。
独立Fraction checkerのisolated `--j21-zero-hx` の最終成功実測は
0.138494508秒、RSS 25196 KiB。producer、S61又はCASをimportしていない。
新point、anchor、normal又はjetを使っていない。

## 1. 原19係数と全六sourceの実inverse

accepted nilpotent sourceのj21_0_hx全19係数を表示式から独立に組み立て、
全raw Jacobiを元の定義から展開した。
rawには十八の非零空間係数があるが、同じ式-2n18n6が二monomialsに現れるため、
旧reportが保存するdistinct residual equationsは十七である。
最初のreviewer実装はmultisetの長さ一致を要求してこの重複で停止した。
失敗時のsource及び契約全文は
`CLASSIFICATION_TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_SEPARATE_REVIEW_IMPLEMENTATION_FAILURE_20261108.json`
へ保存した。最終版では十八全係数を残し、十七distinct式との完全set一致と
唯一の重複-2n18n6を明示確認した。producerの数学的欠陥ではなく、失敗runを成功とは呼ばない。

S=h+xと置く。全source cubicのu係数は-n6 S²。
actual (u,v)尺度amplitudeによる全19sourceの変換を計算し、
その係数が-amplitude*n6 S²、cofactorがv∂u+S∂wのままになることを確認した。
n6非零chartに限りamplitude=1/n6でn6一へ戻せる。n6零は別補集合として保持する。

n6一のaccepted必要十分六source bindingを戻し、全strict10と元350を恒等零と確認した。
全自由parameterはn0,n2,n3,n4,n5,n16。
actual変換

```
U=u-n0h-n3x, V=v,
W=w+(n16+n2)h+(n16+n4)x+n5v
```

を独立に組み立て、det一、全逆、両方向全150、cofactor exact保存を確認した。
この六source正規化ではsource parameterを一つも反転しない。
k=n16、K=-Su+vwとして結果は全域で

```
xu=xv, xv=0, xw=xS,
uv=K-kSv, uw=-K+kSw+2kxv, vw=xv+kS²,
D=hx[-S²u+Svw+xv²].
```

三次のu係数S²と残部v(Sw+xv)は互いに素で、全kについて絶対既約。
五つの余分なsource係数を指定点で落としたsliceではない。

## 2. 全parameterのQeと正確なmodular twist

T=S+e w、Ke=-Tu+vwとして、h行零のquadratic tensorを

```
Qxu=xv, Qxv=e xv, Qxw=xT,
Quv=Ke, Quw=-Ke, Qvw=xv
```

と定義した。全eでstrict10・元350を直接展開して零と確認し、
actualDも hx[xv²+(T+ev)Ke] と一致させた。
row x/xからHamiltonian fieldを、実affine volumeのdivergenceからmodular fieldを独立に作ると

```
Ve=(0,0,v,e v,T),
Ze=(0,-2e x,2v-w-eu,ev-T,2h+4x+2e w)
```

である。両Lie derivative L_Ve Qe、L_Ze Qeは零、[Ve,Ze]=0、div Ve=2eも全係数で確認した。
従ってPe,k=Qe+k Ve wedge Zeは真のPoissonであり、
これを全e,kのliteral tensorから改めてstrict10／元350でも検査した。
全150係数はQQ[e,k]のpolynomialで、kを反転しない。
actualDは厳密に(1+2ke)D(Qe)、e零端点は§1の全k模型へ一致する。

## 3. 半単純正常形の全逆と必要な開

actual coordinate change

```
R=u-v/e, A=v, B=w+(h+x)/e
```

のdet一と両向きinverseを確認した。
全150 tensorを移送して逆にも戻し、Veが
e(A∂A+B∂B)、R方向translation零へ変わることを確認した。
移送tensorのx行はe(1+2ke)x(0,A,B)である。
従ってe(1+2ke)非零での非零tensor scaleによりexact cofactorはdiag(0,1,1)となる。

producerの正規化後全150 rational entriesも、別のFraction AST parserでnum/denへ読み、
各entryを元の移送tensorとpolynomial cross multiplicationで比較した。
全分母は非零scalarを除きe及び1+2keの冪だけであることも確認した。
新たな未知source係数、k又は選んだpointを反転していない。

適用する `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §§4.2–4.3は、
diag(0,1,1)の全十五係数・全binary rank・零S2を含む既受理の全支持定理である。
rank二S2のSNC極限、rank一の実deformation、S2零の残二項形式deformationを全て含むので、
今回のtensorに追加のbinary-rank又は一般係数仮定は要らない。

各固定kについて、e零のcubicは上の絶対既約F0である。
reducible cubicはproper multiplicationの閉像だから、一般eでもcubicは既約。
e(1+2ke)はどの固定kでも恒等零にならず、e零では後者因子は一である。
従って各固定kの曲線の非空稠密開に全半単純定理を適用でき、その全点がMPS40有限閉和に入る。
Pe,kはe零に正則なので全kの端点も同じ閉和内。
source inverseとPGL安定性・閉性によりn6非零の全六source及びその全飽和閉包が被覆される。
e零又は1+2ke零を正規化(7)へ直接代入する推論は使わない。

## 4. 保存と残余

最終reportは `CLASSIFICATION_TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_SEPARATE_REVIEW_20261108.json`、
jobは `two_hyperplanes_j21_zero_hx_modular_separate_review.json`、
入口は `inputs/two_hyperplanes_j21_zero_hx_modular_separate_review_manifest_20261108.json`。
旧checker全FunctionDef ASTと四currentjob bytesは
`inputs/two_quadrics_snc_pre_j21_zero_hx_review_source_contracts_20261108.json`に保存し、不変と確認した。
旧producer・旧report・main・Git・台帳はこのreviewで変更していない。

n6零全枝、j21(h+x,h) p零、三plane1+1+3残源、全E/G及び分類主定理は別義務である。
本結果からambient局所一意性、新成分性、他branchへの陽性伝播を結論しない。
