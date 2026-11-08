# rank一long・非定数gの25次元成分：同一点全ambient独立検証

2026-11-08、quartic_review、PRIMARY E_multiple_hyperplane。**PASS**。
producerのSymPy/general ambient engineをimportせず、既独立の標準ライブラリFraction疎多項式／行列演算で再構成した。
一つの同じ点、同じ全150／350順序、既選択minorだけを使い、新point/profile/minor探索をしていない。

## 点はsourceのg,U,Kから再構成した

座標(h,x,y,z,w)、g=z、U=2z∂z+(3+5z)∂w、
K=1+2z+3w+5z²+7zw+11w²を使う。
divU=2、U(g)=(divU)gを確認し、
F=(U+divU+1)Kを微分と乗算から作った。
F=12+46z+75w+70z²+145zw+33w²と一致した。

入口の元affine六括弧にr=y−x²/2、H=gK、Fを代入し、各degree二以下の式を斉次化して二倍した。
保存150vectorを点の定義として使ったのではない。
元sourceから作った全150係数はanchor JSONとproducer reportの両方へ全係数で一致した。
全十strict Jacobiと全350projective Jacobiは零。
全modularは(0,2h,2x,0,0)、全五Pfとsigned quintic全126係数も一致した。

```
D=−4h²zQ,
Q=hy−x²/2+h²+2hz+3hw+5z²+7zw+11w².
```

QのHessian/2から作る対称行列の行列式は171/32≠0で、rank五である。
一般位数、proper因子又はsource次元を一点の数値から外挿していない。

## 全ambient行列と既minorを独立に再構成した

全150個の元degree二lift方向を一つずつ加え、
J(P+B)−J(P)−J(B)から全350×150微分を作った。
全Euler25列、実GL25列、tensor scalar列も元tensorから独立に作った。
producerの全entryと比較して、D350×150／Euler150×25／actual150×51が全て一致した。

既保存の100／24／47minorの行列式をFractionで再計算し、全て保存値と一致して非零。
各complete kernelを元全行へ代入し、全消滅とfree-coordinate単位行列を確認した。
従って新minor探索なしにrankの上下界が一致し、正確に

```
rank D=100, rank Euler=24, rank actual(Euler+GL+scalar)=47.
```

全350行についてD*actual51=0も確認した。
これは選択submatrixだけのrank推測ではない。
後二rankは同じ保存証明書の安価な追加照合であり、
成分性の弱い矢印に必要なのはDの既100minorと全source下界25である。
actual rank47からのorbit22は**この点だけ**の値で、generic orbit次元を結論しない。

## 全source下界と成分性の接続

全long非定数g sourceの飽和閉包Zが既約・正確に25次元であり、
全prior129と異なるという純証明は、別作者
`RANK_ONE_LONG_NONCONSTANT_ENTRANCE_SEPARATE_REVIEW.md` を固定入力とする。
今回の一点からその全parameter次元を推測していない。
この点は元sourceそのものなのでZに属する。

全150 lift係数空間で、Zのaffine coneのEuler24逆像は既約50次元。
一方、今回の全微分rank100によりZariski tangent次元は50。
局所次元はこの下界と接上界の間で50に等しく、局所環はregular、全空間成分は一意。
Euler24とscalar一を除くとprojective次元25。
閉既約Zも25次元なので、その全成分に等しい。
従って全long非定数g sourceの閉包は、別prior129分離と合わせて新しい25次元全Poisson成分となる。

全sourceの各点のsmooth性、generic orbit、全nonreduced構造は主張しない。
long g零／定数g、short源への別curveや、全分類の完成はこのambient receiptの射程外。
成分表／主稿／Gitへの採録はrootが別に行う。

## 有限実行と保存

fresh job `rank_one_long_nonconstant_ambient_separate_review.json` を実行前gate PASS。
予算180秒／2GiB、実測**0.875117秒、最大RSS17,844KiB、ALL PASS**。
`CLASSIFICATION_RANK_ONE_LONG_NONCONSTANT_AMBIENT_SEPARATE_REVIEW_20261108.json` に
全input hash、inline program全文、正確な既minor determinantとkernel検査結果を保存した。
producer code／report、旧source証明書、全体台帳、Gitは変更していない。
