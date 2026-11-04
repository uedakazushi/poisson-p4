# 正weightの full Bianchi 作用素から来るweight4関係の包含

2026-11-04。追加のcoefficient searchやrank計算を行わず、凍結済みsourceを
加群として結び直す。結論は **旧22関係に、今回の作用素から通常の正weight
多項式倍で得るweight4関係を加えた加群が、rank37の既知加群に等しい** こと。
full kernel80、分母付き・飽和・高weightからの新しい関係は扱わない。

## 1. 実際の63 cubic sourceのwhole kernel

`R0=QQ[a,e,d,1/(aFBC)]` とする。元63 sourceは42個のquadratic×lightと
21個のgenuine weighted cubic、標的は91個のweight3 monomialsである。
T34/Tn座標と各行のunitは受理済みadapterで一致している。

既受理のsix exact QQ row relationsを `r1,...,r6` と書く。deleted rowsは

```
48,55,59,57,61,14
```

で、順次のpivot coefficientsは

```
1,1,1,a,−a,12a/5.
```

後の関係の係数に前のdeleted rowは現れない。triangular determinantは
`−12a³/5` でunitである。従って

`R0^63 / <r1,...,r6> ≅ R0^57`。

実際の91列行列への誘導写像は既証generic rank57を持つ。R0はdomainで、
free57に入る加群kernelはgenericに零、従ってtorsionfree性から零である。
よってこの6本は **R0-linear kernel全体のbasis** である。特殊点での周囲rank57
やcokernel flatnessを仮定していない。このwhole-kernel結論そのものも
`S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW` に受理済み。

任意のkernel vectorに対し、上のdeleted coefficientを対応pivotで割って
順次riを引くと、残る57係数は零となる。従ってここでいう包含はfraction field上
だけではなく、R0上で成立する。未知因子の逆元は不要である。

## 2. 6本のweight3 piecesの実source

full175+tG6作用素の全weight decompositionでは、weight3 pieceは6本、
weight4 pieceは15本、weight3未満は無い。weight3は全てspatial degree5から来て、
保存operatorのpiece indicesは `35,56,71,81,87,90` である。

選択34式をimplicit graph上で零とすると、各weight3 pieceに残るのは

* weight2の6実quadrics ×local weight1 coordinate;
* weight3の21実genuine cubics ×constant

だけである。local weight0 nuisance19/23はT/Uとも恒等的に零。weight1 coordinate
にはlinear graph `T*D*v` だけが寄与し、正weightゆえhigher correctionは無い。
従って各pieceは、同じ63 cubic sourceのR0-linear kernel vectorを与える。
全135 source identityとTn/D/row-unit alignmentは既受理packetで固定されている。

§1より、これら6本それぞれはr1,...,r6のR0-linear combinationである。
その具体的な係数を新たに展開する必要はない。全pieceのsource supportとweightを
今回の小さいdata-binding checkerで照合した。

## 3. light倍と354→319削除

weight4に到達する通常の正weight polynomial multipliersは、weight4 pieceに
対してconstant、weight3 pieceに対して7 light variablesのlinear combination
しかない。Hはweight2なので、この場合のmultiplierに入らない。

r1,...,r5はgenuine rowsの5個のunit deletion relationsである。各々を7 lightで
掛けた35関係は、凍結済み `S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63` の
35 deleted source vectorsそのものであり、354→319 source quotientで零になる。
normalization factorsも同certificateで固定されている。

r6は **product rows0..41だけ** を持つ。従ってそのlight倍は、6quadricsの
quadratic-cofactor関係となる。保存された形ではcofactorはlight²のみで、H係数は
既に零である。仮にHを含む表記で見ても、H係数は六quadricsの定数関係であり、
FULLU上の定数独立性から全係数が零になる。

旧22本はFULLU=`D(a*d*F*B*C*J)`上で、このpure quadratic-cofactor mapの
**whole module kernel** のbasisである。従ってr6の全light倍は旧22内に入る。
特殊fiberでのpointwise rank146を使う推論ではない。

以上より、6個のweight3 piecesの全light倍を319行へ移したものは旧22内に入る。

## 4. 正確な37方向と限定

E22を旧22が張るkernel submodule、D15を今回の直接weight4 pieces15本が張る
submoduleとする。B4を、今回の全作用素piecesから通常の正weight倍により作る
**weight4 source関係**の像とする。§§2–3から

`D15 ⊂ B4 ⊂ E22 + D15`、従って `E22+B4=E22+D15`。

受理済み11quartic pivotsと4cubic pivotsはD15の15方向をaFBC全域で独立にする。
その15列のblockはtriangularで対角I11/I4である。E22はquartic/cubic columnsが
零なので、FULLU上で交わりが零。よって

`E22+B4 = E22 ⊕ D15 ≅ R^37`。

**Bianchi単独のB4がrank37だとは主張しない。** 正確なのは旧22を加えた和が
37であること。これにより、この通常の正weight作用素生成法に旧22を加えても
80方向のgateには **43本不足** することが確定する。

これは与えた175 Schouten rowsと追加tG6 identityのweight≤4から通常に生成する
weight4関係についての限定である。tG6はweight6であり、負weightの除算はしない。
高weightからの別の変換、変数で割った関係、飽和、他のsyzygies、full80 kernel、
actualC122のrank64、全点quartic239、affine exhaustionはこの包含の対象外である。

次の方法へ移る根拠にはなるが、別の80方向構成が不可能という意味ではない。
今回のcheckerは凍結済み係数のlabels/pivots/support/hashを結ぶだけであり、
追加のcoefficient expansion、field、点、GB、rank searchを実行しない。
