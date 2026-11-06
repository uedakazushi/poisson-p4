# 冪零(3,1)型の25次元成分：全支持と独立数値検証

2026-11-07。主課題 E_multiple_hyperplane。
modular 場 `N=x2∂x1+x3∂x2`、c=0 の全 affine 二次支持のうち、
以下で定める n7≠0 の全源の成分性を受理する。
全90係数の線形還元、残る全 Jacobi の当該開上の消去、全周囲接空間を別実装で確認した。
n7=0 の補集合を被覆済みと数えない。

## 全90係数から全10源まで

六つの bracket に定数・一次・二次の全90係数を置き、
`div P=N,[N,P]=0,Q∧N=0` の全係数を独立に作った。
この affine 線形系は QQ 上 rank75、従って解空間は15次元である。
入力 `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` case0 の
全15源は全式を満たし、その係数微分の rank は15なので、線形解空間そのものである。
行の重複や零係数の保存方式による元 producer との行数の相違は、この同値を変えない。

この全15源の十 raw Jacobi を疎な単項式演算で独立再構成した。
次の六式は、実際の Jacobi 係数としてそれぞれ非零有理数倍で現れる。

\[
\begin{gathered}
 (2n_{10}+n_{14}+1)n_7+2n_{11}n_6,\quad
 n_{11}n_7,\quad n_{11}n_5+n_7n_9,\\
 (6n_{10}+2)n_6+3n_{13}n_7,\quad
 (2n_{10}+1)n_5+n_{12}n_7+2n_6n_9,\quad n_5n_9.
\end{gathered}
\]

n7≠0 なら、順に n11=0、n9=0 を得て、残る三式から

\[
 n_{14}=-1-2n_{10},\qquad
 n_{13}=-\frac{2(3n_{10}+1)n_6}{3n_7},\qquad
 n_{12}=-\frac{(2n_{10}+1)n_5}{n_7}
\]

を得る。逆にこの代入後は **全** raw Jacobi が恒等零となる。
従って必要式だけを解いて十分性を省略せず、n7開の全解を
`Gm_n7×A9_(n0,n1,n2,n3,n4,n8,n10,n5,n6)` という既約源で尽くす。

全源の affine Pfaffian は

\[
 -\tfrac12 x_3^2\left(n_5+\tfrac23n_6x_3+n_7x_4\right).
\]

n7≠0 だから全点で一般 rank四である。
斉次化し bivector を6倍した際の Pfaffian も全パラメータで独立確認した。
実退化因子は

\[
 D=2\{h=0\}+2\{x_3=0\}
   +\{n_5h+\tfrac23n_6x_3+n_7x_4=0\}.
\]

三つの一次式は全源で独立である。

## 全周囲の25次元成分性

全10自由パラメータを1とする実点を取り、全150個の二次 lift 係数に対する
全350射影 Jacobi 方程式の微分 D を独立再生成した。
元 producer の Jacobi、軌道、rank helper は import していない。
別作者の疎整数係数ループ `classification_rank3_mixed_component_review.py` を再利用した。

QQ 上の全行還元と非零 minor は

| 行列 | rank |
|---|---:|
| 全350×150 D | 100 |
| Euler鉛直方向 | 24 |
| 鉛直・PGL・scale | 47 |
| 上記と実全10源方向 | 50 |

を与える。全61本の実方向を D が零化することも厳密に確認した。
従って射影全接次元は `150−100−24−1=25`、実既約族の次元下界も `50−24−1=25`。
この点は25次元の滑らかな点で、全空間成分は唯一である。
全10源の PGL 飽和閉包はその成分そのものとなり、n7開の全支持が被覆される。

証明書は `reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_REVIEW_20261107.json`。
全線形系・全源・全接空間を合わせて実測2.982秒、90秒・1536MiBの有限契約内で完了した。
元 producer の有限体非零 minor だけによる下界と区別し、今回の rank は全て QQ 上の等号である。

## 既知成分との相違

25以外の次元の全既知成分とは次元で異なる。
MPS の同次元成分とは、全新源閉包の rank四部分に成立する閉条件「非被約五次因子」と、
MPS が被約 SNC 因子の点を持つことから異なる。
以前からある25次元成分 NEW25 は全閉包に三乗一次因子を持つ。
本源の全点は上記独立な三平面の重複度 (2,2,1) で、三乗一次因子を持たない。
従ってその成分とも一致しない。
この旧 NEW25 の閉条件は、既存の `RANK3_TWZERO_NEW_COMPONENT_EXISTENCE_REVIEW.md` の
全三変数源 h³Q に対する依存監査を入力とし、旧証明書全体を再実行したとはしない。

## 再現

```sh
python3 scripts/classification.py gate research/classification/jobs/rank2_nil_component_review.json
.venv/bin/python code/classification_rank2_nil_component_review.py --output /tmp/rank2-nil31-review.json
```

コード・完全な入力源・producer結果・独立疎演算エンジンの SHA256、
全行列 hash、非零 minor、実10源の全微分係数を公開証明書へ保存した。
