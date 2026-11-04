# Actual E4の全二次initial：cached source2と正確なlinear Schurのみ

2026-11-04。対象は前packetで構成した二つのactual normalizer tensor E4である。
元D0へnormalizer base changeが単射だとは仮定しない。
二次initialからE4の全完成環、socle、depthを結論しない。

## 1. 共通physical R0へ座標を戻す

FULL9 source2の各族57未知数には、1+9+45=55個の実際の係数が保存されている。
このcacheを使って両族の全57式を55係数全てで照合する。
sourceを追加して解くのではなく、既存解でnuと全十二physical outputのeven/odd係数を評価する。

cacheの自由変数は族ごとに異なるactual F9なので、そのまま二族を引かない。
各cacheの九physical R0 readback x(Y)を取り、unit9の固定inverseと既存二次係数から逆map Y(x)を作る。
九つのcomposition identityを全55係数で検査し、両族のnu,a_j,b_jを同じx9へ合成する。
全一次・二次inverse補正を保存する。

## 2. 十一変数と全十四式

typed変数はx9,tauA,tauBで、tauは各族の既存tracefree基底である。

```
tauA²-nuA(x), tauB²-nuB(x),
(a_Aj-a_Bj)(x)+b_Aj(x)tauA-b_Bj(x)tauB, j=全十二physical座標。
```

これら十四式の全次数2係数を構成する。bの一次項だけで二次部分は正確だが、証明書にはbの二次係数も残す。
二つのnuのlinear partを含む全14x11 matrixをQQでSchur消去する。
選んだlinear pivot式に対し、消える変数の二次correctionも同時に計算するので、
残りのequationからそのlinear contributionを落とさない。

全linear式を消した後の残りgeneratorは全てorder2以上である。
したがって、その二次成分のQQ spanは実際のidealの二次initial全体である。
高次cancellationがさらに三次以上を生むことはあっても、未考慮の二次式は生まれない。

## 3. 結果とall-order gateの区別

CA/CBのphysical R0 by F9 determinantはそれぞれ3/128、96。
linear rankは2で、monic二式から(x27,x61)を消せる。そのminorは-20/9である。
残る変数を

```
x108,x109,x112,x113,x117,x118,x119,tauA,tauB
```

とすると、全二次initial spanのrankは**5**、正確に

```
L*(x108,x109,x113,x118,x119),
L=5*x112-20*x117+24*tauA+6*tauB
```

である。五独立式の係数は順に1/12,1/12,1/12,1/12,1/24というunitだけ異なる。
全raw十四式・Schur補正後の十四式・rational factorを証明書に保存する。
二次の有理因数分解だけを行い、GBや高次source jetは使わない。

仮に二次式がL times linear-formsの形でも、これだけで全idealがL times同じlinear-idealとは言えない。
必要な次のgateは、actual N normalizerのpair liftを通じた全次数hypersurface Lと、
そのLで割った残りidealの完全生成・高次relationの同定である。
余った低次元curve成分や閉点embedded部分を二次initialだけで排除しない。

そのactual N pair liftがsmoothなcodimension1埋込としてbindされれば、正式な函数L_formalを選べ、
全十二式はL_formalで割れる。上の五独立linear quotientをunit IFTで消した後、
残る問題は**四変数(x112,x117,tauA,tauB)上の高々七函数のideal**となる。
これらの全函数と残余curveのdepthは未計算であり、単に五quadricsを零にした四平面をactual branchとは呼ばない。
このfactor/quotientのall-order bindingが、次に必要な具体的な有限入力である。

Code `code/ambient_B18_zero_normalizer_quadratic_20261104_v63.py`。
Report/folder `AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63`。
新しいsource jet、48-variable refit、old r/g composition、point、prime、GBは無い。
