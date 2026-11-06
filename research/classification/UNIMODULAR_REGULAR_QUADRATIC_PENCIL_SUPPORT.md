# 全二次pencil最高部のunimodular延長は一般一次像を持たない

2026-11-07。G。`α=J+β+γ` を係数次数0,1,2の閉二形式とし、
Pf(α) が非零定数で、最高部が `γ=df∧dg`、f,gは二次斉次式とする。
第1–5節は pencil span(f,g) が非退化な二次式を含む全場合を扱い、
第6–7節は全memberが退化する補集合を扱う。
**全低次延長のprimitive一次像はW16を優越しない。**
対角化可能性、相異なる固有値、一般pencilを仮定せず全Segre型を保持する。
Jacobian表示を持たない最高部は別に残る。

## 1. 全有限Segre型と線形項の完全性

非退化memberをfとし、対称行列をF,G、自己随伴作用素をT=F⁻¹Gと置く。
複素数体上で、各Jordan block Jm(λ)には反対角対称行列Smを取り、
組(F,G)をblock `(Sm,Sm Jm(λ))` の直和へ正規化できる。
これは非退化対称形式に自己随伴な作用素のJordan分解であり、
一般の対角化可能性を仮定していない。
pencilのGL2変更は固有値へのMöbius作用を与えるので、相異なる固有値が三個以下なら
0,1,−1の部分集合に正規化できる。非退化memberを失わないよう全値を有限に取る。
全体の非零scalarは以後のkernelと評価階数を変えない。

四固有値が相異なる場合は `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第10節で全延長が空と証明済み。
単一固有値でblockが(1,1,1,1)ならγ=0、(2,1,1)ならg−λfはrank一の平方でγに線形因子がある。
相異なる二固有値の重複度(3,1)で三重側が半単純なら、やはりpencil内にrank一の平方がありγに因子がある。
これらは同ノート第15節の全因子付き排除に入る。
残る十型は以下で尽くす。角括弧は同じ固有値のJordan blockをまとめたもの。

| 型 | 非退化Ωの可能性 |
|---|---|
| [4] | なし |
| [3,1] | なし |
| [2,2] | あり、§4 |
| [3];[1] | なし |
| [2,1];[1] | なし |
| [2];[2] | なし |
| [2];[1,1] | なし |
| [1,1];[1,1] | あり、§3 |
| [2];[1];[1] | なし |
| [1,1];[1];[1] | なし |

最高部γを固定し、全24係数の一次二形式βに対して

\[
 d\beta=0,\qquad \beta\wedge\gamma=0
 \tag{1}
\]

を課す。十型全てで全係数行列の階数は16で、解空間は八次元。
その全解は定数一形式u,vにより

\[
 \beta=df\wedge v+u\wedge dg
 \tag{2}
\]

と書かれる。全8列の独立性と(1)への代入をQQ上で検算した。
従ってF1=f+u·x、G1=g+v·x、Ω=J−u∧vと置けば、元の全形式は厳密に

\[
 \alpha=\Omega+dF_1\wedge dG_1.                         \tag{3}
\]

未知の低次係数を部分族へ制限していない。
定数Pfaffianの二次・一次条件は順に

\[
 \Omega\wedge df\wedge dg=0,\qquad
 \Omega\wedge(df\wedge v+u\wedge dg)=0.                  \tag{4}
\]

表の最終列は、最初の条件の全六係数kernelを求め、そのPfaffianを全パラメータで評価した結果である。
八型ではこのPfaffianが恒等零であり、残る二型だけが非零を許す。

## 2. Ωが退化する場合には非退化な全延長が存在しない

この補題はfの非退化性だけを使い、Jordan型を問わない。
Ω=0なら(3)のrankは二以下なので除外される。
rank Ω=2なら線形座標でΩ=dx1∧dx2と置き、U=kerΩの二変数をy1,y2とする。
最初の(4)は

\[
 (\partial_{y_1}f)(\partial_{y_2}g)
 -(\partial_{y_2}f)(\partial_{y_1}g)=0.                  \tag{5}
\]

f非退化より二つの一次式∂y1f,∂y2fは独立であり、互いに素。
UFDにおける次数一の比較から、ある定数tで
`∂_U g=t ∂_U f` が成り立つ。
h=g−tfはUの変数に依存しない。
第二の(4)をU方向へ制限すると

\[
 (df)|_U\wedge(v-tu)|_U=0.
\]

同じ独立性から `(v-tu)|U=0`。
従ってΩ∧u∧v=Ω∧u∧(v−tu)=0であり、
`J=Ω+u∧v` のPfaffianも零となる。これはJの非退化性に矛盾する。
よって全延長ではΩ自身が非退化でなければならない。

表でΩのPfaffianが恒等零だった八型は、全て空である。
以下の二型では、全ての残る低次延長の一次評価が非単射となることを示す。

## 3. 二つの半単純二次元固有空間

pencilの取り直しにより

\[
 f=(x_1^2+x_2^2)/2,\qquad g=(x_3^2+x_4^2)/2
\]

とできる。(4)の二次条件の全解は
`Ω=A dx1∧dx2+B dx3∧dx4`。
Ω非退化だからAB≠0。
一次条件はu3=u4=v1=v2=0を与える。
従ってβの全像は

\[
 \operatorname{span}\{dx_1\wedge v,dx_2\wedge v,
                     u\wedge dx_3,u\wedge dx_4\}
\]

に入る。uとvが共に非零なら二つの二次元空間はu∧vを共有し、次元は高々三。
一方が零ならさらに小さい。従って全ての低次項でrank β≤3である。

## 4. 単一固有値の二つのJordan2 block

固有値を零へ移し、

\[
 f=x_1x_2+x_3x_4,\qquad g=(x_2^2+x_4^2)/2
\]

とする。(4)の二次条件の全解は

\[
 \Omega=t(dx_1\wedge dx_4+dx_2\wedge dx_3)+s\,dx_2\wedge dx_4,
 \qquad\operatorname{Pf}\Omega=t^2.
\]

従ってt≠0。一次条件のx1,x3,x2,x4係数を順に取ると

\[
 v_3=v_1=0,\qquad v_4=u_3,\qquad v_2=u_1.
\]

これを(2)へ代入した全係数は

\[
\begin{aligned}
 \beta_{12}&=2u_1x_2,&\beta_{13}&=0,\\
 \beta_{14}&=u_3x_2+u_1x_4,&\beta_{23}&=-u_3x_2-u_1x_4,\\
 \beta_{24}&=u_3x_1-u_1x_3+u_2x_4-u_4x_2,
 &\beta_{34}&=2u_3x_4.
\end{aligned}
\]

全てx2,x4,u3x1−u1x3の三一次式だけに依存するため、rank β≤3。
u1,u3等を除算せず、全退化境界を含む。

## 5. 分類義務への帰結と検証

primitive一次項は体積形式・Darboux frameによる可逆な線形同型を除いて上のβである。
第3・4節の非単射条件は真の閉条件なので、これらの全源は一般W16を優越しない。
八つの空型、四固有値の空型、全factor型を合わせ、
**非退化memberを持つ全二次pencil最高部の低次延長を処理した。**

公開係数証明書は `CLASSIFICATION_UNIMODULAR_QUADRATIC_PENCIL_20261107.json`、
producerは `code/classification_unimodular_quadratic_pencil.py`、
有限契約は `jobs/unimodular_quadratic_pencil.json`。
全十型の24係数と定数六係数のkernel、期待する全八源、PfaffianをQQ上で検算した（0.290秒）。
統合担当が全Segre表とMöbius規格化、退化Ω補題、二残型の全境界を独立に読解確認した。
独立Fraction証明書 `CLASSIFICATION_UNIMODULAR_QUADRATIC_PENCIL_REVIEW_20261107.json` は
全十型の24係数rank16、全8源、定数六係数kernelとPfaffianを別実装で再構成した（0.031秒）。

## 6. 共通kernelを持つ全pencil

fとgが同じ非零方向を消す場合、座標を取り両者をx1,x2,x3だけの式とする。
γに非定数共通因子がある場合は既受理の因子付き排除に入るので、以下ではgcdγ=1とする。
三変数の二行三列gradient行列の最大minor idealは高さ二である。
Hilbert–Burch完全性により、その全linear syzygyはdf,dgの定数線形結合であり、
非零constant syzygyは存在しない。行列の全成分は一次式なので、生成元の次数を下げてはいない。

一次二形式の全体を

\[
 \beta=\beta_3(x_1,x_2,x_3)+x_4\Omega_3
       +dx_4\wedge\{\eta_1(x_1,x_2,x_3)+x_4\eta_0\}
\]

と分ける。Ω3は定数二形式、η0は定数一形式、η1は一次一形式。
`β∧γ=0`からη1は上のlinear syzygy、η0はconstant syzygyとなる。
従って `η1=a df+b dg`、η0=0であり、dη1=0。
閉性dβ=0からΩ3=0も従う。
よってβはx4に依存せず、全一次評価が非単射である。
この証明は共通kernelの次元や三変数pencilのJordan型によらず、全支持を扱う。

## 7. 共通kernelを持たない全退化pencil

まずこの四変数対称pencilが一つの型で尽くされることを確認する。
generic rankが二以下なら、generic memberを `F=diag(Ir,0)` と置く。
r=2の場合、Gのblock表示を(A,B;Bᵀ,D)とすると、全pencilのrank≤2から
Schur complementが恒等零であり、D=0とBᵀB=0を得る。
二次元非退化対称空間の等方部分空間の次元は一以下なのでrank B≤1。
従ってker B内の非零ベクトルはF,Gの共通kernelとなる。
r≤1でも二つの独立rank一対称形式の和はrank二となるため、共通kernelが存在する。
これは仮定に反するのでgeneric rankは三である。

primitive polynomial kernelベクトルv(t)を取ると、対称性とrank三から

\[
 \operatorname{adj}(F+tG)=r(t)v(t)v(t)^{\mathsf T}
\]

となる。primitive性とUFDによりr(t)も多項式である。
左辺の次数は高々三なのでdeg v≤1。
定数kernelがないから `v(t)=u+t v`、u,vは独立であり、
`Fu=Gv=0,Fv=−Gu≠0`。
対称性によりspan(u,v)は両形式について完全等方。
`Fv` が1となる第三basisを選び、残る第四basisをker(Fv)内に取ると、
第三・第四basisへu,vを加える二つのshearで全余剰cross項と第三対角項を消せる。
従ってpencilは対称Kronecker三blockと一つのscalar blockの直和となる。
scalar blockが零なら共通kernelが生じるので非零。
pencilのGL2変更とscalarを用いると、唯一の形

\[
 f=x_1x_2,\qquad g=x_1x_3+x_4^2/2                         \tag{6}
\]

となる。

(6)についても、全24係数の`dβ=0,βγ=0`のkernelは正確に八次元で、(2)の全源である。
`Ω∧γ=0`の全kernelは `Ω=t dx1∧dx2` のみ。
t=0ならJ=u∧vが退化する。
t≠0なら、一次条件Ω∧β=0のx4,x1係数からu3=u4=0が従う。
従ってΩ∧u∧v=0、Pf J=0となり、やはり非退化性に反する。
よってこの全型の延長は空である。

全係数証明書は `CLASSIFICATION_UNIMODULAR_SINGULAR_PENCIL_20261107.json`、
同じrunnerの `--singular` モードで再現できる（0.072秒）。
契約は `jobs/unimodular_singular_pencil.json`。
旧十型の証明書とproducer bytesは変更せず保持した。
統合担当は第6節のHilbert–Burch議論、第7節の全退化pencilの標準形還元を独立に読解確認した。
独立Fraction証明書 `CLASSIFICATION_UNIMODULAR_SINGULAR_PENCIL_REVIEW_20261107.json` は
全24係数rank16、全八源、定数kernel一次元とPfaffian零を別実装で再構成した（0.005秒）。

これで共通kernelの有無、非退化memberの有無、全因子付き・零最高部の補集合を全て扱い、
**最高部が二つの二次斉次式のJacobianとなる全unimodular支持のG非優越が従う。**
残る義務は、この表示も一次Casimirも持たない最高部である。
