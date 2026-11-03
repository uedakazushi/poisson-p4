# 独立review：全未知残余のdegree43・degree42 scalar表示

2026年11月3日。二つの凍結packet `S61_AFFINE_RESIDUAL43_20261103_v63`、
`S61_AFFINE_RESIDUAL42_20261103_v63` を独立検査し、両方 **ACCEPTED** とする。

## 結論と範囲

betaだけの非空開で、元全134式のideal Iについて、スキームとして正確に

    I+(R42)=I+(R43)=I+(D0,D1,D2,D3,D4)=I+(e)

となる。eは既存の大域known-five splitter、R42とR43は指定finite recipeで正確にdegree42、43である。
未知q、cyclic determinant、quintic Hによる局所化を追加しない。
未知残余の全非被約構造を保つが、その空性・有限性・一般五点尽くしを証明するものではない。

## 新しい独立exact計算

両producerをfresh replayした後、producerのfactor helperやdeterminant helperを使わず、次を再計算した。

* 元134式のQQ係数を、既知平面の45×3 embeddingと42 normal単位列によるpullbackへ変換した。
  quadratic monomialの変換をcacheした独立の1081係数vector engineで全134行を作り、
  選択五平面式のspanを引いた全129行を保存129×42 affine-linear matrixの再乗算と比較した。
  **全129行・全QQ係数が一致**し、純平面項の見落としもない。
* 保存42-row selectorで、既知五点における全42次determinantを別のmodular Gaussian engineで計算した。
  値は1680,23338,12913,13967,26226で一致し、全て非零である。
  保存最高次方向のdegree42係数も31743で一致した。
* 五点のh-Vandermonde inverseからellのsplit weightsを独立に作った。
  全46関数(1と45q座標)のellによる消失と、指定4×4非零minorを確認した。
* 元134式から別のHessian-arrayを作り、degree45/44/43のdeterminant係数をtrace公式で再計算した。
  ellによる三係数は(0,0,15015)、既知因子上のR43 normは17591で一致した。
  raw trace idempotentのdegree45係数は10048で、次数低下が単なるtraceの自動消失でないことも確認した。

測定6.87秒は二つのfresh producerと全独立再構成を含む。有限体値を一般残余の空性とは解釈しない。

## degree43証明の監査

既知finite etale algebra EのE1=span(1,p0,...,p44)は、既存の全16beta普遍三平面定理によりdim≤4である。
中心の非零4-minorと合わせ一般にはdim4であり、ell4=1の四元連立でellをbeta-only有限inverseとして定義できる。
midpoint Jacobian determinantのq-degree45部分はpに依存せず、degree44部分はpにaffine linearである。
ellがE1を殺すのでR43=ell(Dp)の次数は43以下、非零最高次係数で正確に43となる。

分裂後のD_p(p_i)=d_i L_iから、R43の既知因子での値はd_i ell(L_i)となる。
その全値がunitであり、Dp=eDpなのでR43は未知因子で0、既知因子でunitである。
従って(R43)=(e)Aで、ideal等号はschemeとして成立する。
同じeDp=Dpの係数比較とtrace式により、五つのD係数idealも(e)Aに等しい。

## degree42証明の監査

普遍三平面の**実scheme交叉がknown Eそのもの**である既存結果を使う。
全129残差がA(q)nへexactに因数化され、選択42行へのadjugateからR42*n_j∈Iを得る。
normal ideal N=(n)Aの商はEなので、大域分裂A=E×Aresの下でN=(1−e)Aとなる。
したがってR42(1−e)=0。一方R42のEでの全値はunitなので(R42)=(e)Aとなる。

一般betaのmatrixは普遍三平面枠、五つの独立なrestriction、first-normal factor ruleによる有限inverse recipeである。
凍結QQ行列をそのまま全betaの定数matrixと誤解しない。
普遍平面はbeta関数体上に定義され、rational affine originを選べる。
平面内のorigin変更はnを変えず、因数化とknown-point評価を保存するため、凍結toric中心の計算はgeneric recipeの非空性を正しく証明する。

## 凍結review

* `code/s61_affine_residual42_43_review_20261103_v63.py`
* `reports/S61_AFFINE_RESIDUAL42_43_20261103_v63_REVIEW/report.json`

両新packet、全evidence、普遍三平面と既知五点scheme交叉への入力hashを検査した。
次の実問題はI+(R42)のunitであり、その証明は今回のreviewに含まれない。
