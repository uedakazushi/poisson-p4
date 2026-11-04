# 一般五変数源の constrained83 raw-polar 行列

2026-11-04。K=QQ(H,L,s0,q0,t0)。既受理の元134式、q45 の中心 C、
d=Q52−C、N45×20、L51×20 をそのまま用いる。全45未知座標を保ち、
新しい点・素数・GB・Gaussian・Gram・残余の所属計算を行わない。

## 再計算しなかった旧5列

当初の提案 [P F(C),J(C)W0,J(C)W1,polar(W0,d),polar(W1,d)] は、保存済み
five-column-matrix の列と文字通り同一だった。項数は順に
116846,146249,139502,148818,141538、合計692953である。これは旧670個の
source-binding gate で既に確認済みなので再構築しない。既存の小さい核の
共通分母は P290項、分子合計13493項で、旧 Cramer Δ886項・34673項の核とは
異なる。核の列 content-gcd を新たに計算したわけでもない。

代わりに核 T を代入せず、補助変数 alpha20,beta20 を残す。
未知線形係数を (c,alpha20,beta20,n42) の順とし、実行列を

```
[ F(C)   J(C)N   polar_Q(d,N)   polar_Q(d,e_normal42) ]  134行
[   0      L          0                  0           ]   51行
[   0      0          L                  0           ]   51行
```

と定める。法線は現行 pivot q行[3,6,13]以外の42標準座標である。
旧 constant_unit_prefix は L51×20 のみを消去したもので、この236×83行列を
既に消去していたという扱いはしない。

## 新しい正確な係数

source/support preflight の入力は29028項、出力上界53315項、合計82343項。
この結果に従い、許可された一回の QQ 五変数多項式計算を実行した。

| block | 実際の項数 |
|---|---:|
| F(C) | 1483 |
| J(C)N | 17606 |
| polar(d,N) | 16715 |
| raw d-polar42 | 12703 |
| L 二組 | 532 |
| 合計 | 49039 |

9431非零成分、最大一多項式42項、明示的に所有する live 項数の最大77788であった。
独立な二つの有限組立てによって10988個の QQ 多項式恒等式を照合した。
J(C)N は q=C+N alpha の一次係数、polar(d,N) は
q=C+d sigma+N beta の sigma*beta 係数、raw42は sigma*n の係数とも一致する。
q_i² の対角は両方の寄与を残す。追加した L 二組も元入力の全成分と一致する。
native 計算は0.340秒、行列 JSON は約1.31MB。内部 GCD や分数体演算はない。

## 同じ sigma と rank の意味

実 plane は q=C+N*T(a,b)+d sigma であり、Q52 の sigma は1、Q51 は−1/2、
trace中心を持つ toric 三角形は0である。sigma を無正規化の E-pivot 座標と
同一視しない。保存済み W=P*N*T から、制約 L alpha=L beta=0 の核へ戻した
最初の五列は、旧5列全体を非零 parameter P で割ったものになる。
canonical trace-zero 三角形の五 ideal generators は

```
1+sigma−2sigma²+Q0, a+Q1, b+Q2, a*sigma, b*sigma
```

であり、旧五汎関数の単射性と同じ source にある。

受理済み rank L=18 から下102行の rank は36。制約核上の134行は47列である。
五次元 plane image で割った raw d-polar42 は、元 RAW134 residual cross の
d方向の42列と同じである。受理済み RAW134 cross126単射と
d=(H*q0*t0/144)*E0 の実非零 parameter 正規化から単射42が従う。
従って全 constrained 行列の rank は83となる。これは既存の rank 証明からの
表示の同値性であり、新しい83 minor・逆行列・rank 検査を計算したものではない。

## 次の逆計算に渡す範囲

literal QQ 定数成分は187、parameter 単項式成分は1919ある。
source 行・列位置と二部 graph の matching を preflight として保存する。
matching は係数の相殺や消去後の fill を判定しないため、rank 証明や実行可能な
unit-pivot 列の証拠とはしない。定数・単項式の消去もこの packet では未実行。

この実際の小さい行列の parameter 左逆を構成できれば、核 T と大きい5列を
代入せず、plane image を除いた最初の42方向の係数逆へ戻せる。現時点では
その逆も残り84方向も未計算である。一般 fibre の尽くし、normal ideal 所属、
飽和や残余の空性は依然未証明。
