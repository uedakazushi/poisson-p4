# B-free 模型の重い変数の平方係数

2026-11-02。

## 1. 係数環と結論

係数環を K[d], K=QQ(a,e) とする。元135式と既証明の B-free frame
から生成された raw implicit jet を入力とし、ordinary degree2 の
残余135式の H² 係数をすべて抽出する。H は八つの normal 変数のうち
index5、重み2の変数である。

非零係数は21個ある。そのうち source77 と source131 の二つだけで
K[d] における ideal が単位 ideal となる。producer は拡張 Euclid 法で
d 一次以下の二つの係数を求め、正確な Bezout 恒等式1を確認する。

これは K[d] 上の係数の主張である。特に以下に現れる 6e+1 は K の
単元として扱うだけであり、以前の全 core 開上で単元とは主張しない。
全四次関係空間への適用には別証明の238方向全射を用いる必要がある。
長さ65、全 d における二次・三次階数、周囲環や層の飽和除去は未主張。

## 2. raw 係数と実際の残余係数

frame の右逆を U、d に依存しない共通分母を u とし、native には uU
を入力する。ordinary degree2 で native が出力するのは

\[
W_2^{\rm raw}=-uUQ_2,\qquad
A_2^{\rm raw}=uQ_2+JW_2^{\rm raw}.
\]

したがって実際の残余の二次部分は A2_raw/u である。u は K の非零元。
ordinary degree3 以上の項は H² に寄与しないので、raw の最初の
45個の W2 を飛ばした後の135個の A2 だけで十分である。

この抽出は元の B を反転した frame からの移送ではない。
fresh B-free raw の全135係数を読み、後述する二つの因数分解も
得られた係数との正確な比較で検証する。

## 3. 二つの係数と Bezout

\[
\begin{aligned}
F&=500a+90e^2+3e-24,\\
B&=50ad+18de^2-9de-2d+60e-40,\\
C&=100a+18e^2-9e-2,\\
D_0&=50a+9e^2+3e-6,\\
P&=125a^2d-450ae-200a-54e^3-27e^2+33e+6,\\
Q&=15ade-10ad+200a+18e^2+6e-12.
\end{aligned}
\]

raw の H² 係数を h_i とすると、

\[
h_{77}=-\frac{3a^2D_0CFP^2}{3125000000000},\qquad
h_{131}=\frac{3a^2D_0CFBQ}{62500000000}.
\]

report に \(b_{77},b_{131}\in K[d]\) を明示保存し、

\[
b_{77}h_{77}+b_{131}h_{131}=1
\]

を全係数検証する。この二つの b の分母は d に依存せず、因子は
a,F,D0,C,6e+1 のみである。B を分母に持たない。
実際の残余係数に対しては ub_77,ub_131 が Bezout 係数となり、
こちらの恒等式も別に検査する。

p=32003,(a,e)=(1,2) では表示された係数の全分母が可逆であることも
確認する。この指定点の性質を全 core 開へ拡張する主張ではない。

## 4. 再利用 API と検証範囲

```
from s61_binary_hyperplane_Bfree_heavy_square_20261102 import run_heavy_square
certificate = run_heavy_square(wd, raw)
# Fresh source の coeff ファイルも比較する場合:
certificate = run_heavy_square(wd, raw, denominator_input)
```

分母は frozen BFREE_FRAME report の U_common_denominator を使う。
第三引数を渡した場合は、その最初の疎 QQ[a,e,d] 多項式が同じ分母である
ことも検査する。raw は外側の元135式 producer が生成した fresh 入力を
受け取る。入力生成自体をこの小 helper が再実装したとは主張しない。

全135係数、21個の非零位置、二つの raw/actual Bezout 係数を保存する。
正式 source/report/note の hash と生成 raw の hash は分離する。
CLI 出力先は必須で、既存出力は上書きしない。
