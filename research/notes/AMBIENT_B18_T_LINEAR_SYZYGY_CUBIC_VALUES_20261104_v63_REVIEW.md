# 独立監査：実 cubic witness 71

## 受理範囲

実119一次係数 syzygy のうち71番目について、二つの cubic quotient class
の値を独立に計算し、

```
(-251099375/1051776, 0)
```

を確認する。この一つの非零実関係だけで、既知の二次関係の変数倍を除いた
二 channel のうち `e4*z5*r²` class が消え、真の `dim gr³(T)<=1` が従う。
119全値の rank計算を再実行したという主張はしない。残る pure r³ class の
非零性、真の HF3、全長T、nil三乗はこの結果だけでは決まらない。

## 独立な quotient と adjoint

producer の cubic rewrite/rankから dual を採用せず、別途独立監査した
39 homogeneous generator の Gröbner 基底へ、全600 cubic monomial を
逆順 divisor で簡約する。標準項 `e4*z5*r²,e4*r³` への係数として Pi2 を
得る。この Pi2 は、元の37一次式、178実二次式、それらの全変数倍を零に
する。37一次式だけを Gröbner 基底と仮定していない。

凍結 W0の84単位列から U0の逆を固定消去で求め、`Ystar1*U0^-1` の全84列
を作り、逆に U0 を掛けて元 Ystar1 に戻ることを確認する。それから三次
入力 monomial を forward に微分し、正しい output Euler 分母を使って
contraction を scatter する。m0 と実 pi を最後に合成して、全56 source
成分の cubic adjoint を独立に得る。producer の29+1 nonzero配列とは最後に
だけ比較する。

## source71 と pivot 補正

保存された119本の syzygy の71番目と、選択37本の実元768列係数から、八個
の変数に対する元 source column 結合を再構築する。旧640列側と native128列
側について、各変数の定数像が **別々に全56成分で零**であることを確認する。
これにより effective native `(L0+L1)S_B` をこの結合へ使える。

実 a2 と Xp1 から

```
C3=sum l_j X_j2-Xp1*a2
```

を読む。その負号と ordinary coefficient convention を独立に実装し、
mixed monomial は異なる一次因子ごとに一度だけ寄与させる。Hessianの2倍
規約を混ぜない。得られる source functional と既知補正値を元の配列に
照合する。必要な二次 pair は既存 B11 内だけにあることも確認する。

## source値と未知微分の扱い

直接C0値は、既に独立監査済みの108 composite source値を、元nodeIDと
actual CA-F9/B11 frameを固定して再利用する。この依存は full57 inverse
のコピー、同じ CB57source2、実105根からの B11輸送、全単位inverseと opaque
微分消去を独立に検証済みである。今回新しい source求解や値を作らない。

旧 B の typed roots は別の forward B11 algebra で再評価する。各 atom の
定数・一次は保存79原子の実値を使う。108 source cache に入らない atom の
二次係数は五つの独立な形式変数として残し、零という実値を与えない。
sum/product の真の係数規則で運び、source71の最終 contractionでそれらが
全て消えると確認する。これにより producer の oldB adjoint値を単に
読み上げる方法とは異なる検査になる。

native値は独立監査済みの actual G8/monic の B11値を使う。元の L2の項を
各列で零としたという主張ではなく、前述の zero-origin source結合の後に
だけ effective値を適用する。旧 r/g 正則函数を新しい0/0式として再除算
しない。typed rational定数の分母は有理数単位で、他の source単位は上記
受理済みsource依存に固定される。

## cubic class と全次数の論理

X2-only定理では、sourceを二次 pivotで補正すると次数三から始まり、m1,
m2,X3が不要となる。Euler contraction の非R線形差は、一次Koszul複体
`K(u1,...,u7,0)` の H2=0 により d3_linear の一次倍であり、Pi2が零にする
既知二次関係×変数に入る。この projected theorem の範囲で今回の値は
実三次初期関係である。単なる数値配列の非零判定ではない。

残余第一channelが非零、第二channelが零なので、既知三次spanを差し引き
定数で正規化すれば canonical な `e4*z5*r²` が真の initial module に入る。
さらにどの関係があるかは主張しない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_linear_syzygy_cubic_values_review_20261104_v63.py
```

作者scriptをimport/実行せず、独自 Fraction算術で検証する。新source jet、
IFT求解、rank探索、348 scalar行列はない。60秒・2GB・32MiB・150000係数
参照の guard内である。
